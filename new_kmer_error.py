from Bio import SeqIO
from collections import Counter
import itertools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# import sys

# def encodage(sequence):
#     code = {"A": 0, "C": 1, "G": 2, "T": 3}
#     encoded_sequence = 0
#     for i, nucleotide in enumerate(sequence):
#         encoded_sequence += code[nucleotide] * 4 ** (len(sequence) - i - 1)
#     return encoded_sequence



# def parse_fastq2(fastq_file,k):
#     dict_kmers = Counter()
#     with open("fastq_file", "r") as fastq_file:
#         for record in SeqIO.parse(fastq_file, "fastq"):
#             sequence = str(record.seq)
#             dict_kmers.update(count_kmers(sequence, k))


# def parse_fastq3(fastq_file, k):
#     dict_kmers = Counter()
#     with open(fastq_file) as r1:
#         for i, line in enumerate(r1):
#             if i == 3000:
#                 break
#             if i % 4 == 1:  # Adjusted to correctly select sequence lines in a FASTQ file
#                 sequence = line.strip()
#                 if 'N' in sequence:
#                     continue
#                 else:
#                     dict_kmers.update(count_kmers(sequence, k))  # k needs to be defined earlier

#     return dict_kmers

def parse_fastq_test(fastq_file, k):
    dict_kmers = Counter()

    def count_kmers(sequence, k):
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if 'N' not in kmer:
                dict_kmers[kmer] += 1
    
    

    with open(fastq_file, 'r') as file: 
        sequence_count = 0
        for record in SeqIO.parse(fastq_file, "fastq"):
            sequence = str(record.seq)

            if 'N' in sequence:
                continue
            count_kmers(sequence, k)
            sequence_count += 1
            if sequence_count == 100000:
                break
            
        # Read lines in chunks, jumping directly to the sequence line
        # while True:
        #     next_n_lines = list(itertools.islice(file, 1, 4))
        #     if not next_n_lines:
        #         break  # Stop if no more lines to process
        #     sequence = next_n_lines[0].strip()
        #     if 'N' in sequence:
        #         continue
        #     count_kmers(sequence, k)
        #     sequence_count += 1
        #     # if sequence_count == 100000:
        #     #     break

    return dict_kmers



    
    
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
            if sequence_count == 1000:
                break

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


def find_i0_imax2(frequencies):
    """Find the first minimum i0 and the first maximum imax after i0 in the frequency distribution."""
    # Sort the frequencies by count
    sorted_frequencies = sorted(frequencies.items(), key=lambda x: x[0])
    print("sorted frequencies =",sorted_frequencies)
    
    i_min = None
    i_max = None
    prev_freq = None
    
    for count, freq in sorted_frequencies:
        if prev_freq is None:
            prev_freq = freq
        else :
            if freq < prev_freq:
                i_min = (count, freq)
            elif freq > prev_freq:
                i_max = (count, freq)
            prev_freq = freq
       
    
    return i_min, i_max if i_max else i_min  # Return i_max if found, otherwise return i_min

def find_i0_imax2(frequencies):
    """Find the first minimum i0 and the first maximum imax after i0 in the frequency distribution."""
    sorted_frequencies = sorted(frequencies.items(), key=lambda x: x[0])
    print("sorted frequencies =",sorted_frequencies)
    i_min = None
    i_max = None
    prev_freq = None
    looking_for_max = False

    for count, freq in sorted_frequencies:
        if not looking_for_max:
            if i_min is None or freq < sorted_frequencies[i_min][1]:
                i_min = count
            elif prev_freq is not None and freq > prev_freq:
                looking_for_max = True
        if looking_for_max:
            if i_max is None or freq > sorted_frequencies[i_max][1]:
                i_max = count
        prev_freq = freq

    return i_min, i_max
# with open("output_text", "w") as output_file:
#     print(parse_fastq("R1_unpaired.fastqsanger", 12), file=output_file)


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



def correct_read4(fastq_file, kmer_length, i0 = None, quality_scores=None):
    corrected_reads = []
    
    if i0 is None:
        i0, _ =find_i0_imax(frequency_and_distribution(fastq_file, kmer_length))

       
    line = 0
    for record in SeqIO.parse(fastq_file, "fastq"):
        if line == 3000:
            break
        sequence = str(record.seq)
        quality_scores = record.letter_annotations["phred_quality"]
        new_sequence = list(sequence)
        line += 1
        


        dict_kmers = parse_fastq(fastq_file, kmer_length)
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i+kmer_length]
            if dict_kmers[kmer] < i0[1]:
                # Identify the lowest quality base in the k-mer
                lowest_quality_index = i + quality_scores[i:i+kmer_length].index(min(quality_scores[i:i+kmer_length]))
                lowest_quality_index = i + quality_scores[i:i+kmer_length].index(min(quality_scores[i:i+kmer_length]))
                # Replace it with the most probable alternative base
                new_base = reverse_complement(new_sequence[lowest_quality_index])
                new_sequence[lowest_quality_index] = new_base
                # Update quality score to an average high value to denote correction
                quality_scores[lowest_quality_index] = 30

        corrected_reads.append(SeqRecord(Seq(''.join(new_sequence)), id=record.id, description=record.description, letter_annotations={"phred_quality": quality_scores}))

    # Write corrected reads to new FASTQ file
    with open("corrected_output.fastq", "w") as output_file:
        SeqIO.write(corrected_reads, output_file, "fastq")

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
        if line == 3000:
            break
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

def parse_fastq_new(fastq_file):
    with open(fastq_file, 'r') as file:
        line = 0
        while True:
            if line == 3000:
                break
            identifier = file.readline().strip()
            if not identifier:
                break
            sequence = file.readline().strip()
            file.readline()  # skip the '+'
            quality = file.readline().strip()
            yield identifier, sequence, quality
            line += 1

def get_kmer_frequencies(sequences, k):
    kmer_counts = Counter()
    for sequence in sequences:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if 'N' not in kmer:
                kmer_counts[kmer] += 1
    return kmer_counts


def correct_reads(sequence, k, threshold, quality):
    # corrected_reads = []
    # sequences = [seq for _, seq, _ in parse_fastq_new(fastq_file)]
    # kmer_counts = get_kmer_frequencies(sequences, k)
    
    # for identifier, sequence, quality in parse_fastq_new(fastq_file):
    sequence = list(sequence)
    quality_scores = [ord(char) - 33 for char in quality]  # Convert Phred+33 ASCII to quality scores
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer_counts[kmer] < threshold:
            lowest_quality_index = i + quality_scores[i:i+k].index(min(quality_scores[i:i+k]))
            # Dummy logic for base correction: flip to 'A' if not already 'A'
            sequence[lowest_quality_index] = 'A' if sequence[lowest_quality_index] != 'A' else 'G'
            quality_scores[lowest_quality_index] = 30  # Update the quality score arbitrarily to 30
    corrected_sequence = ''.join(sequence)
    corrected_quality = ''.join(chr(q + 33) for q in quality_scores)
    corrected_reads.append((identifier, corrected_sequence, corrected_quality))
    
    return corrected_reads


def write_fastq(corrected_reads, output_file):
    with open(output_file, 'w') as f:
        for identifier, sequence, quality in corrected_reads:
            f.write(f"{identifier}\n{sequence}\n+\n{quality}\n")



# correct_read("R1_paired.fastqsanger", 12)
# with open("output_text_occurences", "w") as output_file:
#     print(frequency_and_distribution("R1_paired.fastqsanger", 12), file=output_file)

# Example usage:
correct_read3('R1_paired.fastqsanger', 12)
# write_fastq(corrected_reads, 'corrected_output.fastq')
