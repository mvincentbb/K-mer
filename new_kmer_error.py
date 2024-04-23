from Bio import SeqIO

from collections import Counter
import itertools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

    

def parse_fastq(fastq_file, k):
    dict_kmers = Counter()

    def count_kmers(sequence, k):
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if 'N' not in kmer:
                dict_kmers[kmer] += 1

    with open(fastq_file, 'r') as file:
        sequence_count = 0
        # Read lines in chunks, jumping directly to the sequence line
        while True:
            next_n_lines = list(itertools.islice(file, 1, 4))
            if not next_n_lines:
                break  # Stop if no more lines to process
            sequence = next_n_lines[0].strip()
            if 'N' in sequence:
                continue
            count_kmers(sequence, k)
            sequence_count += 1
            # if sequence_count == 1000:
            #     break

    return dict_kmers

def frequency_and_distribution(fastq_file, k):
    """Compute frequency and distribution of k-mer counts directly from a Counter.

    Args:
        dict_kmers (Counter): A Counter with k-mers as keys and their counts as values.

    Returns:
        dict: A dictionary where each key is a count and each value is the frequency of this count.
    """
    # Count the frequency of each k-mer count using another Counter
    dict_kmers = parse_fastq(fastq_file, k)
    count_frequencies = Counter(dict_kmers.values())
    count_frequencies = dict(count_frequencies)
    # print(find_i0_imax(count_frequencies))
    return count_frequencies




def find_i0_imax(frequencies):
    """Find the first minimum i0 and the first maximum imax after i0 in the frequency distribution."""
    sorted_frequencies = sorted(frequencies.items(), key=lambda x: x[0])
    print("sorted frequencies =", sorted_frequencies)
    i_min = None
    i_max = None
    prev_freq = None
    looking_for_max = False

    for count, freq in sorted_frequencies:
        if not looking_for_max:
            if i_min is None or freq < i_min[1]:  # Changed from sorted_frequencies[i_min][1] to i_min[1]
                i_min = (count, freq)  # Store the tuple
            elif prev_freq is not None and freq > prev_freq:
                looking_for_max = True
        if looking_for_max:
            if i_max is None or freq > i_max[1]:  # Changed from sorted_frequencies[i_max][1] to i_max[1]
                i_max = (count, freq)  # Store the tuple
        prev_freq = freq

    return i_min, i_max

def reverse_complement(base):
        """
        Compute the reverse complement of a DNA sequence.

        Args:
            sequence (str): DNA sequence.

        Returns:
            str: Reverse complement of the input sequence.
        """
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return complement.get(base, 'N')



def correct_read3(fastq_file, kmer_length, i0_threshold=None):
    # frequencies=
    if i0_threshold is None:
        # kmer_frequencies = parse_fastq(fastq_file, kmer_length)
        # frequencies = {k: v for k, v in kmer_frequencies.items() if v}
        frequencies = frequency_and_distribution(fastq_file, kmer_length)
        i0, _ = find_i0_imax(frequencies)
        i0_threshold = i0[1]

    corrected_reads = []
    dict_kmers = parse_fastq(fastq_file, kmer_length)  # Parse only once

    for line, record in enumerate(SeqIO.parse(fastq_file, "fastq")):
        # if line == 3000:
            # break
        sequence = str(record.seq)
        quality_scores = record.letter_annotations["phred_quality"]
        new_sequence = list(sequence)

        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i+kmer_length]
            if dict_kmers[kmer] < i0_threshold:
                lowest_quality_index = i + quality_scores[i:i+kmer_length].index(min(quality_scores[i:i+kmer_length]))
                new_base = reverse_complement(new_sequence[lowest_quality_index])
                new_sequence[lowest_quality_index] = new_base
                quality_scores[lowest_quality_index] = 30  # Update quality score

        corrected_sequence = ''.join(new_sequence)
        if len(corrected_sequence) != len(quality_scores):
            raise ValueError("Sequence and quality scores length mismatch")
        
        corrected_reads.append(SeqRecord(Seq(corrected_sequence), id=record.id, 
                                         description=record.description, 
                                         letter_annotations={"phred_quality": quality_scores}))

    with open("corrected_output.fastq", "w") as output_file:
        SeqIO.write(corrected_reads, output_file, "fastq")



# correct_read("R1_paired.fastqsanger", 12)
# with open("output_text_occurences", "w") as output_file:
#     print(frequency_and_distribution("R1_paired.fastqsanger", 12), file=output_file)

# Example usage:
correct_read3('R1_paired.fastqsanger', 12)
# write_fastq(corrected_reads, 'corrected_output.fastq')
