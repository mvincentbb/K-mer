from collections import Counter, defaultdict
import argparse



def count_kmers(sequence, k):
    kmer_counts = Counter([sequence[i:i+k] for i in range(len(sequence) - k + 1)])
    print(kmer_counts)
    return kmer_counts


def reverse_complement(sequence):
    """
    Compute the reverse complement of a DNA sequence.

Args:
        sequence (str): DNA sequence.

    Returns:
        str: Reverse complement of the input sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in sequence[::-1])

def find_local_extrema(frequency_counts):
    """
    Find the first local minimum (i0) and the first local maximum after i0 (imax) in the frequency count distribution.

    Args:
        frequency_counts (dict): A dictionary mapping occurrence frequencies to their counts.

    Returns:
        tuple: (i0, imax), where i0 is the first local minimum and imax is the first local maximum after i0.
    """
    sorted_frequencies = sorted(frequency_counts.items())
    i0 = sorted_frequencies[0][0]  # Default value for i0 (minimum occurrence)
    imax = sorted_frequencies[-1][0]  # Default value for imax (maximum occurrence)
    for i in range(1, len(sorted_frequencies) - 1):
        prev_count, curr_count, next_count = (
            sorted_frequencies[i - 1][1],
            sorted_frequencies[i][1],
            sorted_frequencies[i + 1][1],
        )
        if prev_count > curr_count < next_count:
            i0 = sorted_frequencies[i][0]
        if i0 is not None and prev_count < curr_count > next_count:
            imax = sorted_frequencies[i][0]
            break
    return i0, imax


def correct_read(read, kmer_counts, i0, quality_scores=None):
    """
    Correct errors in a read by modifying base calls that create rare k-mers.

    Args:
        read (str): DNA sequence of the read.
        kmer_counts (dict): A dictionary mapping k-mers to their occurrence counts.
        i0 (int): Index of the first local minimum in the frequency count distribution.
        quality_scores (list, optional): List of quality scores for the read.

    Returns:
        tuple: (corrected_sequence, num_corrections)
    """
    k = len(next(iter(kmer_counts)))
    corrected_read = list(read)
    num_corrections = 0

    # Pre-calculate min quality scores for each position to avoid repetitive computation
    min_quality_scores = [min(quality_scores[i:i+k]) if quality_scores else None for i in range(len(read) - k + 1)]

    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        current_kmer_count = kmer_counts.get(kmer, 0)
        
        if current_kmer_count <= i0:
            # Rare k-mer found, attempt to correct
            corrected = False
            for j in range(k):
                if corrected:
                    break
                for base in "ACGT":
                    if base == kmer[j]:
                        continue
                    new_kmer = kmer[:j] + base + kmer[j+1:]
                    if kmer_counts.get(new_kmer, 0) > i0:
                        # Consider quality scores for correction
                        if quality_scores and quality_scores[i+j] <= min_quality_scores[i]:
                            corrected_read[i+j] = base
                            num_corrections += 1
                            corrected = True
                            break
                        elif not quality_scores:
                            corrected_read[i+j] = base
                            num_corrections += 1
                            corrected = True
                            break
            if corrected:
                break  # Remove if correction of all occurrences is required

    return "".join(corrected_read), num_corrections


def correct_read(read, kmer_counts, i0, quality_scores=None):
    """
    Correct errors in a read by modifying base calls that create rare k-mers.

    Args:
        read (str): DNA sequence of the read.
        kmer_counts (dict): A dictionary mapping k-mers to their occurrence counts.
        i0 (int): Index of the first local minimum in the frequency count distribution.
        quality_scores (list, optional): List of quality scores for the read.

    Returns:
        tuple: (corrected_sequence, num_corrections)
    """
    k = len(next(iter(kmer_counts)))
    corrected_read = list(read)
    num_corrections = 0

    # Pre-calculate min quality scores for each position to avoid repetitive computation
    min_quality_scores = [min(quality_scores[i:i+k]) if quality_scores else None for i in range(len(read) - k + 1)]

    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        current_kmer_count = kmer_counts.get(kmer, 0)
        
        if current_kmer_count <= i0:
            # Rare k-mer found, attempt to correct
            corrected = False
            for j in range(k):
                if corrected:
                    break
                for base in "ACGT":
                    if base == kmer[j]:
                        continue
                    new_kmer = kmer[:j] + base + kmer[j+1:]
                    if kmer_counts.get(new_kmer, 0) > i0:
                        # Consider quality scores for correction
                        if quality_scores and quality_scores[i+j] <= min_quality_scores[i]:
                            corrected_read[i+j] = base
                            num_corrections += 1
                            corrected = True
                            break
                        elif not quality_scores:
                            corrected_read[i+j] = base
                            num_corrections += 1
                            corrected = True
                            break
            if corrected:
                break  # Remove if correction of all occurrences is required

    return "".join(corrected_read), num_corrections







def process_fastq(input_file, k, output_file_prefix, i0=None):
    """
    Process the FASTQ file, perform error correction, and write the corrected reads to separate output files for each k value.

    Args:
        input_file (str): Path to the input FASTQ file.
        k (int): Length of k-mers.
        output_file_prefix (str): Prefix for the output FASTQ file names.
        i0 (int, optional): Index of the first local minimum in the frequency count distribution. If not provided, it will be calculated.

    Returns:
        tuple: (total_corrections, original_total_kmers, corrected_total_kmers)
    """
    kmer_counts = defaultdict(int)
    total_corrections = 0
    original_total_kmers = 0
    corrected_total_kmers = 0

    with open(input_file, 'r') as fastq_file:
        sequence = ''
        quality_scores = []
        for i, line in enumerate(fastq_file, start=1):
            if i == 3000:
                break
            if i % 4 == 1:
                header = line.strip()
            elif i % 4 == 2:
                sequence = line.strip()
            elif i % 4 == 0:
                quality_scores = [ord(char) - 33 for char in line.strip()]
                for kmer in count_kmers(sequence, k).keys():
                    kmer_counts[kmer] += 1
                original_total_kmers += sum(kmer_counts.values())
                if i0 is None:
                    i0, _ = find_local_extrema(defaultdict(int, {v: k for k, v in kmer_counts.items()}))
                corrected_sequence, num_corrections = correct_read(sequence, kmer_counts, i0, quality_scores)
                total_corrections += num_corrections
                for kmer in count_kmers(corrected_sequence, k).keys():
                    kmer_counts[kmer] += 1
                corrected_total_kmers += sum(kmer_counts.values())
                output_file = f"{output_file_prefix}_k{k}.fastq"
                with open(output_file, 'a') as corrected_file:
                    corrected_file.write(header + '\n')
                    corrected_file.write(corrected_sequence + '\n')
                    corrected_file.write('+\n')
                    corrected_file.write(''.join(chr(score + 33) for score in quality_scores) + '\n')

    return total_corrections, original_total_kmers, corrected_total_kmers

def run_error_correction(input_file, output_file_prefix, k_values):
    """
    Run the error correction process for different k values and save the corrected reads to separate output files.

    Args:
        input_file (str): Path to the input FASTQ file.
        output_file_prefix (str): Prefix for the output FASTQ file names.
        k_values (list): List of k values to test.
    """
    for k in k_values:
        process_fastq(input_file, k, output_file_prefix)

def generate_report(output_file_prefix, k_values):
    """
    Generate a report summarizing the results of the error correction process for different k values.

    Args:
        output_file_prefix (str): Prefix for the output FASTQ file names.
        k_values (list): List of k values tested.
    """
    with open(f"{output_file_prefix}_report.txt", "w") as report_file:
        report_file.write("Error Correction by k-mer Frequency Report\n\n")
        report_file.write("k\tCorrections\tOriginal k-mers\tCorrected k-mers\tAssembly Quality\n")
        for k in k_values:
            output_file = f"{output_file_prefix}_k{k}.fastq"
            # Read the assembly quality results for this k value
            # (This part will depend on how you evaluate the assembly quality)
            assembly_quality = "N50: 123456, Contigs: 789"
            report_file.write(f"{k}\t...\t...\t...\t{assembly_quality}\n")

        report_file.write("\nConclusions and observations:\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct errors in sequencing reads using k-mer frequency analysis.")
    parser.add_argument("input_file", help="Path to the input FASTQ file.")
    parser.add_argument("k_values", type=int, nargs="+", help="List of k values to test.")
    parser.add_argument("output_file_prefix", help="Prefix for the output FASTQ file names.")
    args = parser.parse_args()

    run_error_correction(args.input_file, args.output_file_prefix, args.k_values)
    generate_report(args.output_file_prefix, args.k_values)
