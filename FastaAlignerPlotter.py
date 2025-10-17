#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
FastaAlignerPlotter.py

Description: This program will output a dotplot of each pairwise alignment 
in a FASTA file with multiple aligned sequences.

User-defined functions: 
    fasta_importer: takes a file path as an argument and returns a dictionary,
where FASTA headers are keys and sequences are values.

Non-standard modules: matplotlib.
    Use pip install matplotlib to install. Numpy is included in the
    installation of matplotlib.

Procedure:
    1: Import libraries
    2: Gather user inputs from the command line
    3: Define user-defined functions
    4: Run program, iterating through sequences and plotting alignments

Input: input file

Usage: ./FastaAligner.py input_file

Version 1.0
Date: 2025-10-17
Name: Oliver Todreas
'''


##### ==================== ###################################################
##### 1: Library importing ###################################################
##### ==================== ###################################################

# The sys module is used to access arguments from the command line.
import sys

# The os module is used to check if files exist and make directories.
import os

# The matplotlib module is used to make plots.
import matplotlib.pyplot as plt

# The numpy module is used to create numpy arrays used by matplotlib.
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
    
# If one argument is passed after the program name, assign it to input_path.
else:
    input_path = sys.argv[1]
    

###### ========================= #############################################
###### 3: User-defined functions #############################################
###### ========================= #############################################

### 3.1: Define FASTA importer function.
### ------------------------------------

# Define fasta_importer that takes a file path as a string and returns a 
# dictionary of the FASTA data.
def fasta_importer(path: str) -> dict:
    
    '''
    This function takes a FASTA file's location as a string, reads it, and 
    returns a dictionary, where keys are headers and values are sequences. It
    ensures the following conditions are met: 
        
        The file exists
        The file is of FASTA format
        The first character in the file is ">" 
        The file contains more than one sequence
        Sequences of equal length
            The sequences must be of equal length because the FASTA file should
            contain globally aligned sequences.
        
    A warning is also printed if invalid characters are detected. Since the
    program does not reject disallowed characters but replaces them with "N",
    any string can theoretically become a sequence.
    '''
    
    # Import libraries used in the function.
    import os
    import numpy as np
    
    # Check that the path passed exists.
    if not os.path.exists(path):
        raise ValueError('The input file does not exist.')

    # Check that the FASTA file is of the correct filetype.
    fasta_exts = ('.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', 
                  '.frn')
    if not path.endswith(fasta_exts):
        raise ValueError('The input file must be a FASTA file.')

    # Create empty dictionary fasta_dict in which to store the file.
    fasta_dict = {}

    # Initialize the variable line_count which will be updated at every line
    # read. Initialize the variable file_startswith_GT to False. This will be
    # changed to True if the input files start with '>'. Set header_read to
    # False and update at every line depending on whether the line is a header.
    first_line = True
    file_startswith_GT = False
    header_read = False
    
    # Initialize the variable invalid_characters_found to False.
    invalid_characters_found = False
    
    # Create a dummy variables head and seq for style consistency.
    head = ''
    seq = ''

    # Define the set valid_chars containing the permitted characters.
    valid_chars = ['A', 'C', 'G', 'T', '-', 'N']

    # Open the FASTA file.
    with open(path, 'r') as f:

        # Initiate while loop to read lines one at a time.
        while True:

            # Strip whitespace characters from the start and end of the line
            # and save the string to the variable line.
            line = f.readline().strip()

            # Check if the line starts with '>', in which case it is a header.
            if line.startswith('>'):

                # If the line starting with '>' is the first line of the file,
                # update file_startswith_GT.
                if first_line:
                    file_startswith_GT = True

                # Add the previous header and complete sequence read before it
                # to the dictionary fasta_dict.
                else:
                    fasta_dict[head] = seq
            
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
                            invalid_characters_found = True
                    
                    # Convert line_list back to a string.
                    line_upper = ''.join(line_list)
                
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
                
            # Check if the file starts with '>'.
            if first_line and not file_startswith_GT:
                raise ValueError('FASTA file corrupted. The FASTA file must '
                                 'start with ">".')
                
            # Update first_line to False.
            first_line = False
                
            # Check if there are no more lines to read. Add the final head and 
            # seq pair to fasta_dict and break the while loop.
            if not line:
                
                # Add the final head and seq pair to fasta_dict.
                fasta_dict[head] = seq
                
                # Warn user if invalid characters were found.
                if invalid_characters_found:
                    print('Warning: invalid characters found in the input '
                          'file. These have been converted to "N".')
                # Break the while loop.
                break

    # Check that all sequences are of the same length. If not, return an error.
    # Loop through the values of the dictionary appending them to seq_lengths.
    seq_lengths = []
    for i in fasta_dict.values():
        seq_lengths.append(len(i))
        
    # Check that more than one sequence was read.
    if len(seq_lengths) < 2:
        raise ValueError('The FASTA file does not contain more than 1 '
                         'sequence.')
        
    # Ensure that the first sequence's is the same as all the other sequences.
    if len(np.unique(seq_lengths)) > 1:
        raise ValueError('Not all sequences in the FASTA file are of the same '
                         'length.')
        
    # Return the dictionary fasta_dict.
    return fasta_dict


##### ================ #######################################################
##### 4: Program logic #######################################################
##### ================ #######################################################

### 4.1: Import data.
### -----------------

# Run the fasta_importer function on the file passed to the program by the 
# user.
fasta_dict = fasta_importer(input_path)


### 4.2: Construct and save plots.
### ------------------------------

# Loop through dictionary entries, plotting each sequence against every other
# sequence in the dictionary once.
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

                    # Check if the bases are a match (matched gaps and matched
                    # N nucletodes are not considered matches).
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
            # bases, show each base, tick, and gridlines.
            if len(seq1) < 30:
                plt.xticks(range(len(seq1)), xticks)
                plt.yticks(range(len(seq2)), yticks)
                plt.grid(alpha=0.5)
                
            else:
                plt.xticks([])
                plt.yticks([])

            # Set x and y axis labels with the protein ID from the header of 
            # the FASTA file, ensuring that only the ID is read by splitting 
            # it into a list and printing only the first element without the
            # '>'.
            id1 = head1.split()[0][1:]
            id2 = head2.split()[0][1:]
            
            ax.set_xlabel(f'Sequence A ({id1})')
            ax.set_ylabel(f'Sequence B ({id2})')

            # Set the plot title.
            plt.title('Dot plot (main diagonal matches in black)')
            
            # Plot the image to the axes. Use color map Greys to ensure that
            # non-matches are white.
            plt.imshow(im, cmap='Greys')
            
            # Save the plot in the folder dotplots.
            if not os.path.exists('dotplots'):
                os.mkdir('dotplots')
            plt.savefig(f'dotplots/{id1}_{id2}.png')