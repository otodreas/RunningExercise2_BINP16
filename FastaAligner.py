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
    
Non-standard modules: None.

Procedure:
    1: Import libraries
    2: Gather user inputs from the command line
    3: Define user-defined functions
    4: Run program, iterating through sequences and calculating their alignment
scores.

Input: input file, parameters [optional], output file [optional]

Usage: ./FastaAligner.py input_file

Version 1.0
Date: 2025-10-16
Name: Oliver Todreas
'''


##### ==================== ###################################################
##### 1: Library importing ###################################################
##### ==================== ###################################################
import sys
import os

##### ============== #########################################################
##### 2: User inputs #########################################################
##### ============== #########################################################

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
                    'output file, "no"/"n" for parameters file): ')
    
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


### 3.2: Define parameter importer function.
### ----------------------------------------

# Define parameter_importer that takes a file path as a string and returns a
# dictionary used for scoring.
def parameter_importer(path: str = None) -> dict:
        
    '''
    This function creates a dictionary of alignment scoring parameters. The
    user can pass a file path to a text document with custom parameter values,
    in which case those parameters will be used to populate the dictionary
    instead. The format for the parameters file is as follows:
        
        match_score = <value>\n
        transition = <value>\n
        transversion = <value>\n
        gap_penalty = <value>
        
    The function ensures that the following conditions are met if a parameter
    file is passed as an argument:
        
        The parameter file exists
        The parameter file is a text file
        Each line of the parameter file contains exactly 1 "="
        A valid parameter name is used to the left of "=" on each line
        Each parameter has a numerical value
    '''
    
    # Assign scores to variables in a dictionary using default values.
    parameters_dict = {'match_score': 1,
                       'transition': -1, # a<->g or c<->t
                       'transversion': -2, # other mutations
                       'gap_penalty': -1}
        
    # Check if the user imported a path
    if path is not None:
            
        # Import libraries used in the function.
        import os
        
        # Check that the parameters file is of the correct filetype.
        if not path.endswith('.txt'):
            raise ValueError('The parameters file must be a text file.')
            
        # Check that the path passed exists.
        if not os.path.exists(path):
            raise ValueError('The parameters file does not exist.')
            
        
    
        # If a parameter document is passed to the program, import the preferred
        # scores and update the dictionary with their values.
        if parameters_path is not None:
            
            # Open the file to read it.
            with open(path, 'r') as file:
                
                # Initiate a while loop to read each line individually.
                while True:
                    
                    # Assign the line to the variable line.
                    line = file.readline()
                    
                    # Check that the end of the file has not been reached.
                    if line:
                        
                        # Find the index of the equals sign in the line.
                        try:
                            equals_sign_index = line.index('=')
                            
                        # If no equals sign is found, raise an error.
                        except ValueError:
                            print('The parameter file is misconfigured. Every '
                                  'line must contain an equals sign.')
    
                        # Check for redundant equals signs.
                        if line.count('=') > 1:
                            raise ValueError('The parameter file is '
                                             'misconfigured. Every line must '
                                             'contain no more than 1 equals sign')
                            
                        # assign the characters before the equals sign to the
                        # variable param_read, stripping leading whitespaces.
                        param_read = line[: equals_sign_index].strip()
                        
                        # Check if the parameter has a match in the parameter
                        # dictionary. If not, raise a formatting error.
                        if param_read not in parameters_dict.keys():
                            raise ValueError('Parameter not found error: the '
                                             'parameter file is misconfigured. '
                                             'Please use the following format:\n'
                                             'match_score = <value>\n'
                                             'transition = <value>\n'
                                             'transversion = <value>\n'
                                             'gap_penalty = <value>')
                        
                        # If the parameter has a match in the parameter dictionary,
                        # try adding the value to the right of the equals sign to
                        # that dictionary key.
                        else:
                            
                            # Assign the value following the equals sign to the 
                            # variable value as a float.
                            try:
                                value = float(line[equals_sign_index
                                                   + 1:].strip())
                            
                            # Throw an error if the value is not convertable to 
                            # float.
                            except ValueError:
                                print('Parameter value error: all parameter '
                                      'values must be numerical.')
                                
                            # Update the dictionary with the new value.
                            parameters_dict[param_read] = value
                        
                    # Break the while loop when there are no more lines to be
                    # read.
                    else:
                        break
            
    # Return the dictionary
    return parameters_dict


##### ================ #######################################################
##### 4: Program logic #######################################################
##### ================ #######################################################

### 4.1: Import data and parameters if necessary.
### ---------------------------------------------

# Run the fasta_importer function on the file passed to the program by the 
# user.
fasta_dict = fasta_importer(input_path)
parameters_dict = parameter_importer(parameters_path)

        
### 4.2: Calculate scores for each pairwise alignment.
### --------------------------------------------------

# Assign an empty list to the variable score_summary.
score_summary = ''

# Loop through the keys in the dictionary, calculating each possible pairwise
# score.
for i, key1 in enumerate(fasta_dict.keys()):
    for j, key2 in enumerate(fasta_dict.keys()):

        # Only calculate scores when i < j to ensure the same sequence does
        # not get scored with itself or that duplicate scores are reported.
        if i < j:

            # Assign sequences to variables and normalize case. Ensure that
            # they are the same length, otherwise interrupt the program.
            seq_a, seq_b = fasta_dict[key1], fasta_dict[key2]
            if len(seq_a) != len(seq_b):
                raise ValueError('sequences are not of the same length. '
                                 'program interrupted')

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


### 4.3: Write output file.
### -----------------------
with open(output_path, 'w') as f:
    f.write(score_summary)