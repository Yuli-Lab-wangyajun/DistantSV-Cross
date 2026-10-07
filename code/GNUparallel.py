#!/usr/bin/env python3

import argparse
import re
import stat
import sys
from collections import defaultdict
from pathlib import Path

MAX_CORES = 60

def safe_name(name):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return name.strip("_") or "unknown_species"

def parse_lastz_script(script_path, species_names):
    suffix = "_lastz.sh"
    filename = script_path.name
    if not filename.endswith(suffix):
        return None
    stem = filename[:-len(suffix)]
    for species in sorted(species_names, key=len, reverse=True):
        species_suffix = f"_{species}"
        if not stem.endswith(species_suffix):
            continue
        chromosome = stem[:-len(species_suffix)]
        if chromosome and species:
            return chromosome, species
    return None

def get_requested_cores(script_path):
    patterns = (
        re.compile(r"^\s*#SBATCH\s+--ntasks-per-node=(\d+)\s*$"),
        re.compile(r"^\s*#SBATCH\s+--ntasks-per-node\s+(\d+)\s*$"),
    )
    with script_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            for pattern in patterns:
                match = pattern.match(line)
                if match:
                    cores = int(match.group(1))
                    if cores < 1:
                        raise ValueError(
                            f"Invalid core count in {script_path}: {cores}"
                        )
                    if cores > MAX_CORES:
                        raise ValueError(
                            f"Single task exceeds {MAX_CORES} cores: "
                            f"{script_path} ({cores} cores)"
                        )
                    return cores
    raise ValueError(
        f"Cannot find #SBATCH --ntasks-per-node in: {script_path}"
    )

def scan_lastz_scripts(target_dir, reference_species, species_names):
    tasks_by_species = defaultdict(list)
    skipped = []
    for script_path in sorted(target_dir.glob("*_lastz.sh")):
        parsed = parse_lastz_script(script_path, species_names)
        if parsed is None:
            skipped.append(
                f"Skip unrecognized filename: {script_path.name}"
            )
            continue
        chromosome, species = parsed
        if reference_species and species == reference_species:
            skipped.append(
                f"Skip reference species script: {script_path.name}"
            )
            continue
        cores = get_requested_cores(script_path)
        tasks_by_species[species].append(
            {
                "chromosome": chromosome,
                "species": species,
                "cores": cores,
                "script_path": script_path.resolve(),
            }
        )
    return tasks_by_species, skipped

def split_by_cores(tasks, max_cores):
    ordered = sorted(
        tasks,
        key=lambda task: (
            -task["cores"],
            task["chromosome"],
            str(task["script_path"]),
        ),
    )
    groups = []
    used_cores = []
    for task in ordered:
        selected = None
        for index, current_cores in enumerate(used_cores):
            if current_cores + task["cores"] <= max_cores:
                selected = index
                break
        if selected is None:
            groups.append([task])
            used_cores.append(task["cores"])
        else:
            groups[selected].append(task)
            used_cores[selected] += task["cores"]

    for group in groups:
        group.sort(
            key=lambda task: (
                task["chromosome"],
                str(task["script_path"]),
            )
        )
    return groups

def split_aux_scripts(script_paths, max_scripts):
    script_paths = sorted(script_paths, key=lambda path: path.name)
    if not script_paths:
        return []
    num_scripts = len(script_paths)
    num_groups = (num_scripts + max_scripts - 1) // max_scripts
    base_size = num_scripts // num_groups
    remainder = num_scripts % num_groups
    groups = []
    start = 0
    for group_index in range(num_groups):
        group_size = base_size + (1 if group_index < remainder else 0)
        groups.append(script_paths[start:start + group_size])
        start += group_size
    return groups

def get_sbatch_lines(script_path, num_cores, job_name):
    sbatch_lines = []
    has_cores = False
    with script_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#SBATCH"):
                continue
            stripped = line.strip()
            if re.match(r"#SBATCH\s+-J\s+", stripped):
                sbatch_lines.append(f"#SBATCH -J {job_name}\n")
            elif re.match(
                r"#SBATCH\s+--ntasks-per-node(?:=|\s+)",
                stripped,
            ):
                sbatch_lines.append(
                    f"#SBATCH --ntasks-per-node={num_cores}\n"
                )
                has_cores = True

            else:
                sbatch_lines.append(line)
    if not has_cores:
        sbatch_lines.append(f"#SBATCH --ntasks-per-node={num_cores}\n")

    return sbatch_lines

def write_combined_script(target_dir, species, tasks, part_no, total_parts):
    species_name = safe_name(species)
    if total_parts == 1:
        job_name = f"combined_lastz_{species_name}"
    else:
        job_name = f"combined_lastz_{species_name}_part{part_no}"

    output_path = target_dir / f"{job_name}_job.sh"
    total_cores = sum(task["cores"] for task in tasks)
    if total_cores > MAX_CORES:
        raise RuntimeError(
            f"Internal error: {output_path.name} requests "
            f"{total_cores} cores"
        )
    sbatch_lines = get_sbatch_lines(
        tasks[0]["script_path"],
        total_cores,
        job_name,
    )
    with output_path.open("w", encoding="utf-8") as out:
        out.write("#!/bin/bash\n")
        out.writelines(sbatch_lines)
        out.write("\n")
        out.write("set -euo pipefail\n\n")
        out.write(f"# Species: {species}\n")
        out.write(f"# Chromosome scripts: {len(tasks)}\n")
        out.write(
            f"# Total cores: {total_cores}; "
            f"maximum allowed: {MAX_CORES}\n\n"
        )

        out.write("PIDS=()\n")
        out.write("PID_CORES=()\n")
        out.write("RUNNING_CORES=0\n")
        out.write("FAILED=0\n\n")

        out.write("wait_for_first() {\n")
        out.write("    local pid=${PIDS[0]}\n")
        out.write("    local cores=${PID_CORES[0]}\n")
        out.write('    if ! wait "$pid"; then\n')
        out.write("        FAILED=1\n")
        out.write("    fi\n")
        out.write('    PIDS=("${PIDS[@]:1}")\n')
        out.write('    PID_CORES=("${PID_CORES[@]:1}")\n')
        out.write("    RUNNING_CORES=$((RUNNING_CORES - cores))\n")
        out.write("}\n\n")

        for task in tasks:
            task_cores = task["cores"]
            script_path = task["script_path"]

            out.write(
                f'echo "Run {task["chromosome"]}: '
                f'{task_cores} cores"\n'
            )
            out.write(
                f"while (( RUNNING_CORES + {task_cores} "
                f"> {total_cores} )); do\n"
            )
            out.write("    wait_for_first\n")
            out.write("done\n")
            out.write(f'bash "{script_path}" &\n')
            out.write('PIDS+=("$!")\n')
            out.write(f'PID_CORES+=("{task_cores}")\n')
            out.write(
                f"RUNNING_CORES=$((RUNNING_CORES + "
                f"{task_cores}))\n\n"
            )

        out.write("while (( ${#PIDS[@]} > 0 )); do\n")
        out.write("    wait_for_first\n")
        out.write("done\n\n")
        out.write("if (( FAILED != 0 )); then\n")
        out.write(
            '    echo "At least one LASTZ task failed." >&2\n'
        )
        out.write("    exit 1\n")
        out.write("fi\n")

    mode = output_path.stat().st_mode
    output_path.chmod(
        mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
    )

    return output_path

def write_combined_aux_script(
    target_dir,
    script_type,
    script_paths,
    part_no,
    total_parts,
):
    if total_parts == 1:
        job_name = f"combined_{script_type}"
    else:
        job_name = f"combined_{script_type}_part{part_no}"

    output_path = target_dir / f"{job_name}_job.sh"
    total_cores = len(script_paths)

    if total_cores > MAX_CORES:
        raise RuntimeError(
            f"Internal error: {output_path.name} requests "
            f"{total_cores} cores"
        )

    sbatch_lines = get_sbatch_lines(
        script_paths[0],
        total_cores,
        job_name,
    )

    with output_path.open("w", encoding="utf-8") as out:
        out.write("#!/bin/bash\n")
        out.writelines(sbatch_lines)
        out.write("\n")
        out.write("set -euo pipefail\n\n")
        out.write(f"# Script type: {script_type}\n")
        out.write(f"# Child scripts: {len(script_paths)}\n")
        out.write(
            f"# Total cores: {total_cores}; "
            f"maximum allowed: {MAX_CORES}\n\n"
        )

        out.write("PIDS=()\n")
        out.write("FAILED=0\n\n")

        for script_path in script_paths:
            out.write(f'echo "Run: {script_path.name}"\n')
            out.write(f'bash "{script_path.resolve()}" &\n')
            out.write('PIDS+=("$!")\n\n')

        out.write('for pid in "${PIDS[@]}"; do\n')
        out.write('    if ! wait "$pid"; then\n')
        out.write("        FAILED=1\n")
        out.write("    fi\n")
        out.write("done\n\n")

        out.write("if (( FAILED != 0 )); then\n")
        out.write(
            f'    echo "At least one {script_type} task failed." >&2\n'
        )
        out.write("    exit 1\n")
        out.write("fi\n")

    mode = output_path.stat().st_mode
    output_path.chmod(
        mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
    )

    return output_path

def write_batch_submit(target_dir, combined_scripts):
    submit_path = target_dir / "batch_lastz_submit.sh"
    with submit_path.open("w", encoding="utf-8") as out:
        out.write("#!/bin/bash\n")
        out.write("set -euo pipefail\n\n")
        out.write('echo "Start submitting combined jobs..."\n\n')
        for script_path in combined_scripts:
            out.write(f'echo "Submitting: {script_path.name}"\n')
            out.write(f'sbatch "{script_path}"\n')
            out.write("sleep 1\n\n")
        out.write('echo "All combined jobs submitted."\n')
    mode = submit_path.stat().st_mode
    submit_path.chmod(
        mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
    )
    return submit_path

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Scan LASTZ, multiz and SV scripts, then create combined "
            "species-level and type-level SLURM scripts."
        )
    )
    parser.add_argument(
        "target_dir",
        help="Directory containing *_lastz.sh, *_multiz.sh and *_SV.sh",
    )
    parser.add_argument(
        "-r",
        "--reference-species",
        required=True,
        help="Reference species name to exclude, for example: mu_mu",
    )
    parser.add_argument(
        "-s",
        "--species",
        required=True,
        help=(
            "Comma-separated list of all species names. Species names may "
            "contain underscores, for example: rano,ra_no,mu_mu"
        ),
    )
    parser.add_argument(
        "-g",
        "--group-size",
        type=int,
        default=MAX_CORES,
        help="Maximum cores/scripts per combined script; maximum is 60",
    )
    parser.add_argument(
        "-c",
        "--cores",
        type=int,
        default=None,
        help="Kept for compatibility only and ignored",
    )
    args = parser.parse_args()
    species_names = sorted(
        {
            species.strip()
            for species in args.species.split(",")
            if species.strip()
        },
        key=len,
        reverse=True,
    )
    if not species_names:
        parser.error(
            "-s/--species must contain at least one species name"
        )
    if args.reference_species not in species_names:
        parser.error(
            f"Reference species '{args.reference_species}' "
            "must be included in -s/--species"
        )
    if args.group_size < 1:
        parser.error("-g/--group-size must be at least 1")
    if args.group_size > MAX_CORES:
        print(
            f"Warning: -g {args.group_size} exceeds {MAX_CORES}; "
            f"using {MAX_CORES}.",
            file=sys.stderr,
        )
        args.group_size = MAX_CORES
    if args.cores is not None:
        print(
            "Warning: -c/--cores is ignored. LASTZ core counts are read "
            "from each LASTZ script; multiz/SV use one core per script.",
            file=sys.stderr,
        )
    target_dir = Path(args.target_dir).expanduser().resolve()
    if not target_dir.is_dir():
        parser.error(f"Directory does not exist: {target_dir}")
    tasks_by_species, skipped = scan_lastz_scripts(
        target_dir,
        args.reference_species,
        species_names,
    )
    multiz_scripts = sorted(target_dir.glob("*_multiz.sh"))
    sv_scripts = sorted(target_dir.glob("*_SV.sh"))
    for message in skipped:
        print(f"Warning: {message}", file=sys.stderr)
    if not tasks_by_species and not multiz_scripts and not sv_scripts:
        print(
            "Error: no non-reference *_lastz.sh, *_multiz.sh, "
            "or *_SV.sh scripts were found.",
            file=sys.stderr,
        )
        sys.exit(1)
    old_script_patterns = (
        "combined_lastz_*_job.sh",
        "combined_multiz_job.sh",
        "combined_multiz_part*_job.sh",
        "combined_SV_job.sh",
        "combined_SV_part*_job.sh",
    )
    for pattern in old_script_patterns:
        for old_script in target_dir.glob(pattern):
            old_script.unlink()
    combined_scripts = []
    print(f"Scan directory: {target_dir}")
    print("=" * 60)

    for species in sorted(tasks_by_species):
        tasks = tasks_by_species[species]
        groups = split_by_cores(tasks, args.group_size)
        print(
            f"Species: {species}; chromosome scripts: {len(tasks)}; "
            f"combined scripts: {len(groups)}"
        )
        for part_no, group in enumerate(groups, start=1):
            group_cores = sum(task["cores"] for task in group)

            script_path = write_combined_script(
                target_dir,
                species,
                group,
                part_no,
                len(groups),
            )
            combined_scripts.append(script_path)

            print(
                f"  part {part_no}/{len(groups)}: "
                f"{len(group)} scripts, {group_cores} cores -> "
                f"{script_path.name}"
            )
    if multiz_scripts:
        multiz_groups = split_aux_scripts(
            multiz_scripts,
            args.group_size,
        )

        print(
            f"Multiz scripts: {len(multiz_scripts)}; "
            f"combined scripts: {len(multiz_groups)}"
        )

        for part_no, group in enumerate(multiz_groups, start=1):
            script_path = write_combined_aux_script(
                target_dir,
                "multiz",
                group,
                part_no,
                len(multiz_groups),
            )

            combined_scripts.append(script_path)

            print(
                f"  multiz part {part_no}/{len(multiz_groups)}: "
                f"{len(group)} scripts, {len(group)} cores -> "
                f"{script_path.name}"
            )
    if sv_scripts:
        sv_groups = split_aux_scripts(
            sv_scripts,
            args.group_size,
        )

        print(
            f"SV scripts: {len(sv_scripts)}; "
            f"combined scripts: {len(sv_groups)}"
        )

        for part_no, group in enumerate(sv_groups, start=1):
            script_path = write_combined_aux_script(
                target_dir,
                "SV",
                group,
                part_no,
                len(sv_groups),
            )

            combined_scripts.append(script_path)

            print(
                f"  SV part {part_no}/{len(sv_groups)}: "
                f"{len(group)} scripts, {len(group)} cores -> "
                f"{script_path.name}"
            )

    submit_path = write_batch_submit(target_dir, combined_scripts)
    print("=" * 60)
    print(f"Generated combined scripts: {len(combined_scripts)}")
    print(f"Submission script: {submit_path}")

if __name__ == "__main__":
    main()