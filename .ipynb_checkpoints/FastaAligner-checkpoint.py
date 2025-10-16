import sys

# FASTA importer function

# Define fasta_importer that takes a file path as a string and returns a dictionary of the FASTA data.
def fasta_importer(path: str) -> dict:

    # Create empty dictionary fasta_dict in which to store the file.
    fasta_dict = {}

    # Initiate the variable first_line to True, which will be set to False once the first line is read.
    first_line = True

    # Define the set valid_chars containing the permitted characters.
    valid_chars = ['A', 'C', 'G', 'T', '-']

    # Open the FASTA file.
    with open(path, 'r') as f:

        # Initiate while loop to read lines one at a time.
        while True:

            # Strip whitespace characters from the start and end of the line and save the string to the variable
            # line.
            line = f.readline().strip()

            # Check if the line starts with '>', in which case it is a header.
            if line.startswith('>'):

                # Add the previous header and complete sequence read before it to the dictionary fasta_dict.
                if not first_line:
                    fasta_dict[head] = seq

                # If the first line of the file is being read, do not add anything to the dictionary since
                # a sequence has not been read yet. Set first_line to False.
                else:
                    first_line = False
            
                # Assign the line to the variable head and set header_read to True.
                head = line
                header_read = True

            # If the line does not start with '>' it is a sequence.
            else:
                
                # Check that the sequence contains only valid characters.
                if sum(line.upper().count(c) for c in valid_chars) != len(line):
                    raise TypeError(f'Invalid characters detected in sequence "{head.split()[0][1:]}". '
                                    f'Program interrupted, please re-align sequences.')

                # Check if the previously read line was a header. If so, assign the current line to the
                # variable seq.
                if header_read:
                    seq = line.upper()

                # If the previously read line was a sequence, add the current line to the previous line to fix
                # sequences split by newline characters.
                else:
                    seq += line.upper()

                # Assign the variable header_read to False if the line read did not start with a '>'.
                header_read = False
        
            # Check if there are no more lines to read. Add the final head and seq pair to fasta_dict and break
            # the while loop.
            if not line:
                fasta_dict[head] = seq
                break

    return fasta_dict


# ----------------------
# Load data into memory.
# ----------------------

# Assign file path to path. 
path = 'exampledata/score.example.fna'
output_path = 'exampledata/output_fasta.txt'

fasta_dict = fasta_importer(path)

# --------------------------           
# Generate alignment scores.
# --------------------------

# Assign scores to variables.
match_score = 1
transition = -1 # a<->g or c<->t
transversion = -2 # other mutations
gap_penalty = -1
gap_match = 0 # this script assumes that if two sequences both have a gap in the same position, that
              # position of the parallel alignment adds 0 to the alignment score.

# Assign an empty list to the variable score_summary.
score_summary = ''

# Loop through the keys in the dictionary, calculating each possible pairwise score.
for i, key1 in enumerate(fasta_dict.keys()):
    for j, key2 in enumerate(fasta_dict.keys()):

        # Only calculate scores when i < j to ensure the same sequence does not get scored with itself or that
        # duplicate scores are reported.
        if i < j:

            # Assign sequences to variables and normalize case. Ensure that they are the same length, otherwise
            # interrupt the program.
            seq_a, seq_b = fasta_dict[key1], fasta_dict[key2]
            if len(seq_a) != len(seq_b):
                raise ValueError('sequences are not of the same length. program interrupted')

            # Initialize alignment scoring variables.
            identity = 0
            gaps = 0
            score = 0
            alignment_len = len(seq_a)

            # Zip the sequences to compare alignment at each position and loop through nucleotide by nucleotide.
            for a, b in zip(seq_a, seq_b):

                # Check if characters are identical.
                if a == b:

                    # If they represent matching nucleotides, update score and identity count.
                    if '-' not in a:
                        score += match_score
                        identity += 1
                        #print(a, b, match_score)

                        # IMPORTANT: matched gaps ('-') are ignored
                
                # If characters are not identical:
                else:

                    # Check if the unidentical characters represent a nucleotide transition and update the score.
                    if len({a, b} & {'A', 'G'}) == 2 or len({a, b} & {'C', 'T'}) == 2:
                        score += transition
                        #print(a, b, 'transition')

                    # The unidentical characters do not represent a nucleotide transition.
                    else:

                        # Check if the unidentical characters represent a gap and update the score and gap count.
                        if '-' in (a, b):
                            score += gap_penalty
                            gaps += 1
                            #print(a, b, 'gap')

                        # The last possible combinatinon of characters represents a transversion. Update the score.
                        else:
                            score += transversion
                            #print(a, b, 'transversion')

            # Add the score summary to the string score_summary.
            score_summary += (f'{key1[1:]}-{key2[1:]}: '
                              f'Identity: {identity}/{alignment_len} ({round(identity/alignment_len*100, 1)}%), '
                              f'Gaps: {gaps}/{alignment_len} ({round(gaps/alignment_len*100, 1)}%), '
                              f'Score={score}\n')


# ------------------
# Write output file.
# ------------------
with open(output_path, 'w') as f:
    f.write(score_summary)