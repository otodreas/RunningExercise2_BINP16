#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
FastaAlignerPlotter.py

Description: This program will output a dotplot of each pairwise alignment 
in a FASTA file with multiple aligned sequences. This function produces a
dotplot using matplotlib's "imshow" function, which reads a 2D numpy array as
an image. The advantages to this are that the function is a neat solution to 
plotting neighboring squares with different colors without having to draw
polygons manually. The disadvantage is that the image needs bitmap coordinates,
in which the y-axis is inverted. Therefore, the sequence plotted on the y-axis
is read backwards.

User-defined functions: 
    fasta_importer: takes a file path as an argument and returns a dictionary,
where FASTA headers are keys and sequences are values.
    Note: separate documentation is provided inside the user-defined functions.

Non-standard modules: matplotlib.
    Use pip install matplotlib to install. Numpy is included in the
    installation of matplotlib.

Procedure:
    1: Import libraries
    2: Gather user inputs from the command line
    3: Define user-defined functions
    4: Run program, iterating through sequences and plotting alignments

The program addresses the following potential errors:
    Exactly 1 filepath must be passed to the function. It must exist and be of
the appropriate filetype
    Other potential errors are handled inside user-defined functions. Separate
documentation is provided for them

Input: input file

Usage: ./FastaAligner.py input_file

Version 1.0
Date: 2025-10-17
Name: Oliver Todreas
'''


# ---------------
# Library imports
# ---------------

# The sys module is used to access arguments from the command line.
import sys

# The os module is used to check if files exist and make directories.
import os

# The matplotlib module is used to make plots.
import matplotlib.pyplot as plt

# The numpy module is used to create numpy arrays used by matplotlib.
import numpy as np

# Import custom library.
from importers import fasta_importer


# -----------
# User inputs
# -----------

# Check the number of arguments passed to the program.
if len(sys.argv) <= 1:
    raise IndexError('Too few arguments passed. You must pass the path to '
                     'the input file to the program.')
    
# If three or more, raise an error.
elif len(sys.argv) > 2:
    raise IndexError('Too many arguments passed. Please pass only one file.')
    
# If one argument is passed after the program name, assign it to input_path.
else:
    input_path = sys.argv[1]


# -------------
# Program logic
# -------------

# Run the fasta_importer function.
fasta_dict = fasta_importer(input_path)

# Plot each nucleotide against every other nucleotide in a sequence pair.
for i, item1 in enumerate(fasta_dict.items()):
    for j, item2 in enumerate(fasta_dict.items()):

        # Assign variables head and seq for each dictionary.
        head1 = item1[0]
        head2 = item2[0]
        seq1 = item1[1]
        seq2 = item2[1]

        # Ensure only desired sequence pairs get plotted.
        if i < j:

            # Create the list im_list to store the data to be plotted.
            im_list = []

            # Check match status
            for l, base2 in enumerate(seq2[::-1]):

                # Assign a blank list to the variable row.
                row = []

                # Assign position values.
                for k, base1 in enumerate(seq1):

                    # Check if the bases are a match.
                    if base2 == base1 and base2 in 'ACGT':

                        # Diagonal/non-diagonal match.
                        if k + l + 1 == len(seq1):
                            row.append(2)

                        else:
                            row.append(1)
                            
                    # Non-match. 
                    else:
                        row.append(0)

                # Append row to im_list once the entire row has been read.
                im_list.append(row)

            # Convert im_list to a numpy array after reading.
            im = np.array(im_list)

            # Assign variables xticks and yticks to lists of each sequence.
            xticks, yticks = list(seq1), list(seq2[::-1])

            # Assign fig and ax using a matplotlib figure and axes.
            fig, ax = plt.subplots()

            # Show each base, tick, and gridlines.
            if len(seq1) < 30:
                plt.xticks(range(len(seq1)), xticks)
                plt.yticks(range(len(seq2)), yticks)
                plt.grid(alpha=0.5)
                
            else:
                plt.xticks([])
                plt.yticks([])

            # Set x and y axis labels.
            id1 = head1.split()[0][1:]
            id2 = head2.split()[0][1:]
            ax.set_xlabel(f'Sequence A ({id1})')
            ax.set_ylabel(f'Sequence B ({id2})')

            # Set the plot title.
            plt.title('Dot plot (main diagonal matches in black)')
            
            # Plot the image to the axes.
            plt.imshow(im, cmap='Greys')
            
            # Save the plot in the folder dotplots.
            if not os.path.exists('dotplots'):
                os.mkdir('dotplots')
                
            # Save the figure.
            plt.savefig(f'dotplots/{id1}_{id2}.png', dpi=fig.get_dpi())