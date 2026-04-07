#!/usr/bin/env python3

import os
import sys
import glob

SUFFIXES = ['_lastz.sh', '_multiz.sh', '_SV.sh']
MAX_CORES = 60
ENV_SETUP = [
    'source ~/miniconda3/bin/activate',
    'conda activate test'
]

def get_sbatch_lines(script_path, num_cores, job_suffix):
    """Extract #SBATCH lines from script, modify core count and job name"""
    sbatch_lines = []
    with open(script_path, 'r') as f:
        for line in f:
            if line.startswith('#SBATCH'):
                if '--ntasks-per-node=' in line:
                    sbatch_lines.append(f'#SBATCH --ntasks-per-node={num_cores}\n')
                elif '-J ' in line:
                    sbatch_lines.append(f'#SBATCH -J combined_{job_suffix}\n')
                else:
                    sbatch_lines.append(line)
    return sbatch_lines

def generate_combined_script(suffix, script_files, target_dir):
    """Generate combined SLURM script for scripts with specified suffix"""
    num_scripts = len(script_files)
    if num_scripts == 0:
        return

    num_cores = min(num_scripts, MAX_CORES)
    print(f"Processing {suffix}: Found {num_scripts} files, using {num_cores} cores")

    # Get SBATCH header from the first script
    first_script = script_files[0]
    sbatch_lines = get_sbatch_lines(first_script, num_cores, suffix[1:-3])

    # Build output file path
    output_name = f'combined_{suffix[1:-3]}_job.sh'
    output_path = os.path.join(target_dir, output_name)

    with open(output_path, 'w') as out:
        # Write shebang
        out.write('#!/bin/bash\n')
        # Write SBATCH lines
        out.writelines(sbatch_lines)
        out.write('\n')
        # Write environment setup
        for line in ENV_SETUP:
            out.write(line + '\n')
        out.write('\n')
        # Write parallel command
        out.write(f'# Use GNU Parallel to execute all {suffix} scripts in parallel\n')
        out.write(f'parallel -j {num_cores} \'bash {{}}\' ::: \\\n')
        # Write absolute paths of all scripts
        abs_paths = [os.path.abspath(sf) for sf in script_files]
        for i, path in enumerate(abs_paths):
            if i < len(abs_paths) - 1:
                out.write(f'    "{path}" \\\n')
            else:
                out.write(f'    "{path}"\n')
        out.write('\n# Wait for all tasks to complete (parallel waits by default)\n')
        out.write('wait\n')

    # Add executable permission
    os.chmod(output_path, 0o755)
    print(f"Generated: {output_path}")

def main():
    # Determine target directory
    target_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

    # Process each suffix separately
    for suffix in SUFFIXES:
        pattern = os.path.join(target_dir, f'*{suffix}')
        script_files = glob.glob(pattern)
        script_files.sort()
        generate_combined_script(suffix, script_files, target_dir)

    print("All combined scripts have been generated.")

if __name__ == '__main__':
    main()