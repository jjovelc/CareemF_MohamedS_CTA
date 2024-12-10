import sys

# File paths from command line arguments
all_reads_file = sys.argv[1]
na_reads_file = sys.argv[2]
output_file = sys.argv[3]

# Dictionary to store the headers from all_reads_file
all_reads_dict = {}

print("Hashing long IDs")
# Read the all_reads_file and store only headers in a dictionary
with open(all_reads_file, 'r') as all_reads:
    for line in all_reads:
        line = line.strip()
        if line.startswith('>'):
            # Store the header, using the ID up to the first space as the key
            short_id = line.split(' ')[0]
            all_reads_dict[short_id] = line

# Open na_reads and output files
with open(na_reads_file, 'r') as na_reads, open(output_file, 'w') as output:
    for line in na_reads:
        line = line.strip()
        if line.startswith('>'):
            # Write the corresponding header from all_reads if it exists
            output.write(f'{all_reads_dict.get(line, line)}\n')
        else:
            # Write the sequence line as it is
            output.write(f'{line}\n')
