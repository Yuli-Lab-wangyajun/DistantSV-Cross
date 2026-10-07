import sys

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <fileA> <fileB> <output_file>")
        sys.exit(1)
    
    fileA_path = sys.argv[1]
    fileB_path = sys.argv[2]
    output_file_path = sys.argv[3]
    
    # Read file A and store as a dictionary
    a_dict = {}
    with open(fileA_path, 'r') as fileA:
        for line in fileA:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                filename = parts[0]
                label = parts[1]
                # Store the third column if present; otherwise store empty string
                third_col = parts[2] if len(parts) >= 3 else ""
                a_dict[filename] = (label, third_col)
    
    # Process file B and write to output
    with open(fileB_path, 'r') as fileB, open(output_file_path, 'w') as out_file:
        # Skip the header line of file B
        header = next(fileB)
        # Write the header line for the output file
        out_file.write("Chromosome\tNewRegionStart\tNewRegionEnd\tLabel\tThirdColumnFromA\n")
        
        for line in fileB:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            filename_b = parts[0]
            if filename_b in a_dict:
                col4 = parts[3]   # Column 4 (Chromosome)
                col7 = parts[6]   # Column 7 (NewRegionStart)
                col8 = parts[7]   # Column 8 (NewRegionEnd)
                label, third_col_a = a_dict[filename_b]
                out_file.write(f"{col4}\t{col7}\t{col8}\t{label}\t{third_col_a}\n")

if __name__ == "__main__":
    main()