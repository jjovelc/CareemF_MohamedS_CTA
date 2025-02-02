import re
import sys

def extract_numbers(line):
    """Extract numbers from BUSCO result line"""
    # Try to extract percentage
    percent_match = re.search(r'(\d+\.\d+)%', line)
    percent = float(percent_match.group(1)) if percent_match else None
    
    # Try to extract count
    count_match = re.search(r'^\s*(\d+)\s+', line)
    count = int(count_match.group(1)) if count_match else None
    
    return count, percent

def parse_busco_results(input_text):
    # Dictionary to store results for all samples
    all_results = {}
    
    # Split the input text into individual BUSCO result blocks
    busco_blocks = input_text.strip().split("_______________________________")
    
    for block in busco_blocks:
        if not block.strip():
            continue
            
        # Extract sample name
        sample_match = re.search(r'/BUSCO/(.+?)_transcripts', block)
        if not sample_match:
            continue
        sample_name = sample_match.group(1)
        
        # Initialize results dictionary for this sample
        results = {
            'Complete': {'count': 0, 'percent': 0},
            'Single-copy': {'count': 0, 'percent': 0},
            'Duplicated': {'count': 0, 'percent': 0},
            'Fragmented': {'count': 0, 'percent': 0},
            'Missing': {'count': 0, 'percent': 0}
        }
        
        # Extract numbers for each category
        lines = block.split('\n')
        for line in lines:
            if 'Complete BUSCOs (C)' in line:
                results['Complete']['count'], _ = extract_numbers(line)
            elif 'Complete and single-copy BUSCOs (S)' in line:
                results['Single-copy']['count'], _ = extract_numbers(line)
            elif 'Complete and duplicated BUSCOs (D)' in line:
                results['Duplicated']['count'], _ = extract_numbers(line)
            elif 'Fragmented BUSCOs (F)' in line:
                results['Fragmented']['count'], _ = extract_numbers(line)
            elif 'Missing BUSCOs (M)' in line:
                results['Missing']['count'], _ = extract_numbers(line)
            elif 'C:' in line:
                # Extract percentages from summary line
                percents = re.findall(r'(\d+\.\d+)%', line)
                if len(percents) >= 3:
                    results['Complete']['percent'] = float(percents[0])
                    results['Single-copy']['percent'] = float(percents[1])
                    results['Duplicated']['percent'] = float(percents[2])
        
        all_results[sample_name] = results
    
    # Write consolidated results
    with open("consolidated_busco_results.txt", "w") as f:
        # Write header
        samples = sorted(all_results.keys())
        f.write("Category\t" + "\t".join(samples) + "\n")
        
        # Write counts
        f.write("Counts\n")
        for category in ['Complete', 'Single-copy', 'Duplicated', 'Fragmented', 'Missing']:
            f.write(f"{category}\t")
            counts = [str(all_results[sample][category]['count']) for sample in samples]
            f.write("\t".join(counts) + "\n")
        
        # Write percentages
        f.write("\nPercentages\n")
        for category in ['Complete', 'Single-copy', 'Duplicated', 'Fragmented', 'Missing']:
            f.write(f"{category}\t")
            percentages = [f"{all_results[sample][category]['percent']:.1f}" for sample in samples]
            f.write("\t".join(percentages) + "\n")

def main():
    if len(sys.argv) != 2:
        print("Usage: python parse_busco_reports.py <results_file>")
        sys.exit(1)
    
    results_file = sys.argv[1]
    
    with open(results_file, "r") as results:
        content = results.read()
    
    parse_busco_results(content)

if __name__ == "__main__":
    main()
