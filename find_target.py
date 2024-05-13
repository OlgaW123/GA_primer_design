#from Bio import SeqIO
import subprocess
import os

def run_blat(query_seq):
    # Write the query sequence to a temporary file
    with open("temp_query.fa", "w") as file:
        file.write(">query\n" + query_seq)

    # Run BLAT
    blat_command = ["blat", "hg38.2bit", "temp_query.fa", "blat_output.psl"]
    subprocess.run(blat_command, check=True)

    # Parse the BLAT output to get the best hit coordinates
    best_hit = None
    with open("blat_output.psl", "r") as file:
        lines = file.readlines()
        if len(lines) > 5:  # Skipping header lines in PSL file
            best_hit = lines[5].strip().split()

    # Clean up temporary files
    #os.remove("temp_query.fa")
    #os.remove("blat_output.psl")
    for index,field in enumerate(best_hit):
        print(f"Index {index}: {field}")
    t_starts = list(map(int, best_hit[20].strip(',').split(',')))
    q_sizes = list(map(int, best_hit[18].strip(',').split(',')))
    return best_hit[13], t_starts[0], t_starts[-1] + q_sizes[-1]  # Start and end coordinates of the best hit

def fasta_to_string(file_name):
    try:
        with open(file_name,'r') as file:
            content = file.read()
        return content
    except FileNotFoundError:
        print("File not found!")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def extract_sequence(genome, chrom, start, end):
    # Adjust coordinates to include 1000 nucleotides on each side
    start = max(0, start - 1000)
    print(start)
    end = end + 1000
    print(end)
    # Load the genome using twoBitToFa (if not already in a parseable format)
    subprocess.run(["twoBitToFa", f"-seq={chrom}", f"-start={start}", f"-end={end}", "hg38.2bit", "hg38_cut.fa"], check=True)
    # Extract the sequence using Biopython
    extended_seq = fasta_to_string("hg38_cut.fa")
    return extended_seq

def main():
    query_sequence = "ACTCTTCTGGTCCCCACAGACTCAGAGAGAACCCACCATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGTGAGGCTCCCTCCCCTGCTCCGACCCGGGCTCCTCGCCCGCCCGGACCCACAGGCCACCCTCAACCGTCCTGGCCCCGGACCCAAACCCCACCCCTCACTCTGCTTCTCCCCGCAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGGTGAGCGGCGGGCCGGGAGCGATCTGGGTCGAGGGGCGAGATGGCGCCTTCCTCTCAGGGCAGAGGATCACGCGGGTTGCGGGAGGTGTAGCGCAGGCGGCGGCTGCGGGCCTGGGCCGCACTGACCCTCTTCTCTGCACAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAAGCTGGAGCCTCGGTAGCCGTTCCTCCTGCCCGCTGGGCCTCCCAACGGGCCCTCCTCCCCTCCTTGCACCGGCCCTTCCTGGTCTTTGAATAAAGTCTGAGTGGGCAGCA"  # Replace with your sequence
    chrom, start, end = run_blat(query_sequence)
    extended_sequence = extract_sequence("hg38.fa", chrom, start, end)

    # Write the extended sequence to a new file
    with open("output_sequence.fa", "w") as file:
        file.write(">extended_region\n" + extended_sequence)

if __name__ == "__main__":
    main()
