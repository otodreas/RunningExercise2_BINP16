#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
FastaAligner.py

Description: This program will output each pairwise alignment score in a FASTA
file with multiple aligned sequences. It uses predifined scoring rules that can
be altered by the user.

User-defined functions: 
    fasta_importer: takes a file path as an argument and returns a dictionary,
where FASTA headers are keys and sequences are values.
    parameter_importer: takes a file path as an argument and returns a
dictionary where score types are keys and scores are values.
    Note: separate documentation is provided inside the user-defined functions.
    
Non-standard modules: None.

Procedure:
    1: Import libraries
    2: Gather user inputs from the command line
    3: Define user-defined functions
    4: Run program, iterating through sequences and calculating their alignment
scores.

The program addresses the following potential errors:
    Exactly 1, 2, or 3 filepaths must be passed to the function. They must
exist and be of the appropriate filetype
    Other potential errors are handled inside user-defined functions. Separate
documentation is provided for them

Input: input file, parameters [optional], output file [optional]

Usage: ./FastaAligner.py input_file parameters output_file

Version 1.0
Date: 2025-10-17
Name: Oliver Todreas
'''


# ---------------
# Library imports
# ---------------

# The sys module is used to access arguments from the command line.
import sys

# The os module is used to check if files exist.
import os

# Import custom functions.
from importers import fasta_importer, parameter_importer


# -----------
# User inputs
# -----------

# Initialize parameters_path and output_path.
parameters_path = None
output_path = 'output_fasta.txt'

# Check the number of arguments passed to the program.
# If nothing is passed after the program name, throw an error.
if len(sys.argv) <= 1:
    raise IndexError('Too few arguments passed. You must pass the path to '
                     'the input file to the program.')
    
# If four or more, raise an error.
elif len(sys.argv) > 4:
    raise IndexError('Too many arguments passed. Please pass no more than '
                     'three arguments to the program.')
    
# If one argument is passed after the program name, assign it to input_path.
elif len(sys.argv) == 2:
    input_path = sys.argv[1]

# If two, assign the first to input_path and ask the user which type of file
# the second argument is referring to.
elif len(sys.argv) == 3:
    input_path = sys.argv[1]
    argtype = input('The identity of the last argument is ambiguous. Are you '
                    'providing a path to an output file? ("yes"/"y" for '
                    'output file, "no"/"n" for parameter file): ')
    
    if argtype == 'yes' or argtype == 'y':
        output_path = sys.argv[2]
        
    elif argtype == 'no' or argtype == 'n':
        parameters_path = sys.argv[2]
        
    else:
        raise ValueError('Invalid answer. Please start over.')
    
# If three, assign them to input_path, parameters_path, and output_path.
else: # the only remaining possibility is that 4 arguments were passed.
    input_path, parameters_path, output_path = sys.argv[1:4]
    
    
# Check that the output file is a .txt file.
if not output_path.endswith('.txt'):
    raise ValueError('The output file must be a .txt file.')
    
# Check if the output path already exists and warn the user of overwriting
# risks unless the output path is left at its default.
if os.path.exists(output_path) and output_path != 'output_fasta.txt':
    overwrite = input(f'The output path "{output_path}" already contains a '
                      f'file. Do you want to overwrite it ("yes"/"y" or '
                      f'"no"/"n")? ')
        
    if overwrite == 'yes' or overwrite == 'y':
        pass
        
    elif overwrite == 'no' or overwrite == 'n':
        raise Exception('Program interrupted.')
        
    else:
        raise ValueError('Invalid answer. Please start over.')


# -------------
# Program logic
# -------------

# Run the fasta_importer function on the file passed to the program by the 
# user.
fasta_dict = fasta_importer(input_path)
parameters_dict = parameter_importer(parameters_path)

# Assign an empty list to the variable score_summary.
score_summary = ''

# Loop through the keys in the dictionary, calculating each possible pairwise
# score.
for i, key1 in enumerate(fasta_dict.keys()):
    for j, key2 in enumerate(fasta_dict.keys()):

        # Only calculate scores when i < j to ensure the same sequence does
        # not get scored with itself or that duplicate scores are reported.
        if i < j:

            # Assign sequences to variables.
            seq_a, seq_b = fasta_dict[key1], fasta_dict[key2]

            # Initialize alignment scoring variables.
            identity = 0
            gaps = 0
            score = 0
            alignment_len = len(seq_a)
            disregarded_positions = 0

            # Zip the sequences to compare alignment at each position and loop
            # through nucleotide by nucleotide.
            for a, b in zip(seq_a, seq_b):

                # Check if either nucleotide is the character N. If so, do not
                # change any scores but update the disregarded_positions 
                # variable
                if 'N' in (a, b):
                    disregarded_positions += 1

                # Check if characters are identical.
                elif a == b:

                    # If they represent matching nucleotides, update score and
                    # identity count.
                    if '-' not in a:
                        score += parameters_dict['match_score']
                        identity += 1
                    
                    # If the match is a gap, do not update any scores but 
                    # update the disregarded_positions variable.
                    else:
                        disregarded_positions += 1
                        
                # If characters are not identical:
                else:

                    # Check if the unidentical characters represent a
                    # nucleotide transition and update the score.
                    if (len({a, b} & {'A', 'G'}) == 2
                        or len({a, b} & {'C', 'T'}) == 2):
                        score += parameters_dict['transition']

                    # The unidentical characters do not represent a nucleotide
                    # transition.
                    else:

                        # Check if the unidentical characters represent a gap
                        # and update the score and gap count.
                        if '-' in (a, b):
                            score += parameters_dict['gap_penalty']
                            gaps += 1

                        # The last possible combinatinon of characters
                        # represents a transversion. Update the score.
                        else:
                            score += parameters_dict['transversion']

            # Add the score summary to the string score_summary.
            denominator = alignment_len - disregarded_positions
            identity_pct = round(identity / denominator * 100, 1)
            gaps_pct = round(gaps / denominator * 100, 1)
            
            score_summary += (f'{key1[1:]}-{key2[1:]}: '
                              f'Identity: {identity}/{denominator} '
                              f'({identity_pct}%), '
                              f'Gaps: {gaps}/{denominator} '
                              f'({gaps_pct}%), '
                              f'Score={score}\n')

# Write output file.
with open(output_path, 'w') as f:
    f.write(score_summary)