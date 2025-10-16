##### ==================== ###################################################
##### 1: Library importing ###################################################
##### ==================== ###################################################
import sys
import os
import matplotlib.pyplot as plt
import numpy as np


##### ============== #########################################################
##### 2: User inputs #########################################################
##### ============== #########################################################

# Check the number of arguments passed to the program.
# If nothing is passed after the program name, throw an error.
if len(sys.argv) <= 1:
    raise IndexError('Too few arguments passed. You must pass the path to '
                     'the input file to the program.')
    
# If three or more, raise an error.
elif len(sys.argv) > 2:
    raise IndexError('Too many arguments passed. Please pass only one file.')

# Check that the first file passed exists.
elif not os.path.exists(sys.argv[1]):
    raise ValueError('The input file does not exist. Ensure that the input '
                     'file is passed as the first argument!')
    
# If one argument is passed after the program name, assign it to input_path.
else:
    input_path = sys.argv[1]
    

###### ========================= #############################################
###### 2: User-defined functions #############################################
###### ========================= #############################################

### 2.1: Define FASTA importer function.
### ------------------------------------

# Define fasta_importer that takes a file path as a string and returns a 
# dictionary of the FASTA data.
def fasta_importer(path: str) -> dict:

    # Create empty dictionary fasta_dict in which to store the file.
    fasta_dict = {}

    # Initiate the variable first_line to True, which will be set to False
    # once the first line is read.
    first_line = True
    
    # Create a dummy variables head and seq for style consistency.
    head = ''
    seq = ''

    # Define the set valid_chars containing the permitted characters.
    valid_chars = ['A', 'C', 'G', 'T', '-']

    # Open the FASTA file.
    with open(path, 'r') as f:

        # Initiate while loop to read lines one at a time.
        while True:

            # Strip whitespace characters from the start and end of the line
            # and save the string to the variable line.
            line = f.readline().strip()

            # Check if the line starts with '>', in which case it is a header.
            if line.startswith('>'):

                # Add the previous header and complete sequence read before it
                # to the dictionary fasta_dict.
                if not first_line:
                    fasta_dict[head] = seq

                # If the first line of the file is being read, do not add
                # anything to the dictionary since a sequence has not been 
                # read yet. Set first_line to False.
                else:
                    first_line = False
            
                # Assign the line to the variable head and set header_read to
                # True.
                head = line
                header_read = True

            # If the line does not start with '>' it is a sequence.
            else:
                
                # Convert line to uppercase.
                line_upper = line.upper()
                
                # Check that the sequence contains only valid characters. If
                # it does not, replace invalid characters with N.
                if sum(line_upper.count(c) for c in valid_chars) != len(line):
                    
                    # Store the characters from line_upper in a list to 
                    # iterate through.
                    line_list = list(line_upper)
                    for i, l in enumerate(line_list):
                        
                        # Update the character at position i to N if it is not
                        # valid.
                        if l not in valid_chars:
                            line_list[i] = 'N'
                    
                    # Convert line_list back to a string.
                    line_upper = str(line_list)
                
                # Check if the previously read line was a header. If so, 
                # assign the current line to the variable seq.
                if header_read:
                    seq = line_upper

                # If the previously read line was a sequence, add the current 
                # line to the previous line to fix sequences split by newline 
                # characters.
                else:
                    seq += line_upper

                # Assign the variable header_read to False if the line read 
                # did not start with a '>'.
                header_read = False
        
            # Check if there are no more lines to read. Add the final head and 
            # seq pair to fasta_dict and break the while loop.
            if not line:
                fasta_dict[head] = seq
                break

    # Return the dictionary fasta_dict.
    return fasta_dict


##### ================ #######################################################
##### 3: Program logic #######################################################
##### ================ #######################################################

### 3.1: Import data.
### -----------------

# Run the fasta_importer function on the file passed to the program by the 
# user.
fasta_dict = fasta_importer(input_path)


### 3.1: Construct and save plots.
### ------------------------------

# Loop through dictionary entries, plotting each sequence against every other sequence in the dictionary once.
for i, item1 in enumerate(fasta_dict.items()):
    for j, item2 in enumerate(fasta_dict.items()):

        # Assign variables head and seq for each dictionary to improve
        # readability.
        head1 = item1[0]
        head2 = item2[0]
        seq1 = item1[1]
        seq2 = item2[1]

        # Ensure that a sequence does not get plotted against itself or that
        # two sequences get plotted twice.
        if i < j:

            # Create the list im_list to store the data to be plotted.
            im_list = []

            # Loop through the bases of both sequences to evaluate if they 
            # represent on or off axis matches. Start with the sequence that
            # will be plotted on the y axis. Since this program uses
            # plt.imshow to show the data, the data must be constructed top-
            # left down, because the y axis on image data increases as it goes
            # down.
            for l, base2 in enumerate(seq2[::-1]):

                # Assign a blank list to the variable row.
                row = []

                # Loop through the bases of the sequence plotted on the x 
                # axis, checking for matches. Diagonal matches are encoded as
                # 2, other matches are encoded as 1, non-matches are encoded
                # as 0.
                for k, base1 in enumerate(seq1):

                    # Check if the bases are a match (matched gaps are not
                    # considered matches).
                    if base2 == base1 and base2 in 'ACGT':

                        # If the two indecies and 1 sum to the length of the
                        # sequence, the match is a diagonal match, since the
                        # sequence represented on the y axis is being read top
                        # down. Append the appropriate value to row.
                        if k + l + 1 == len(seq1):
                            row.append(2)

                        else:
                            row.append(1)
                            
                    # In the case of a non-match, append 0 to row. 
                    else:
                        row.append(0)

                # Append row to im_list once the entire row has been read.
                im_list.append(row)

            # Convert im_list to a numpy array so that it can be handled by 
            # plt.imshow.
            im = np.array(im_list)

            # Assign variables xticks and yticks to lists of each sequence, 
            # with sequence 2 reversed.
            xticks, yticks = list(seq1), list(seq2[::-1])

            # Assign fig and ax using a matplotlib figure and axes.
            fig, ax = plt.subplots()

            # Check the length of the sequences. If they are shorter than 30 
            # bases, show each base and tick.
            # Otherwise, hide bases and ticks for readability.
            if len(seq1) < 30:
                plt.xticks(range(len(seq1)), xticks)
                plt.yticks(range(len(seq2)), yticks)
                
            else:
                plt.xticks([])
                plt.yticks([])

            # Set x and y axis labels with the protein ID from the header of 
            # the FASTA file, ensuring that only the ID is read by splitting 
            # it into a list and printing only the first element without the
            # '>'.
            id1 = head1.split()[0][1:]
            id2 = head2.split()[0][1:]
            
            ax.set_xlabel(f'Sequence 1 ({id1})')
            ax.set_ylabel(f'Sequence 2 ({id2})')

            # Set the plot title.
            plt.title('Dot plot (main diagonal matches in black)')
            
            plt.imshow(im, cmap='Greys')
            
            # Save the plot in the folder dotplots.
            if not os.path.exists('dotplots'):
                os.mkdir('dotplots')
            plt.savefig(f'dotplots/{id1}_{id2}.png')