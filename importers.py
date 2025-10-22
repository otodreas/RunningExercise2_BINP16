#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -----------------------
# FASTA importer function
# -----------------------

# Define fasta_importer.
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

    # Initialize variables.
    first_line = True
    file_startswith_GT = False
    header_read = False
    invalid_characters_found = False
    head = ''
    seq = ''

    # Define the set valid_chars containing the permitted characters.
    valid_chars = ['A', 'C', 'G', 'T', '-', 'N']

    # Open the FASTA file.
    with open(path, 'r') as f:

        # Initiate while loop to read lines one at a time.
        while True:

            # Assign stripped line to the variable line.
            line = f.readline().strip()

            # Check if the line starts with '>', in which case it is a header.
            if line.startswith('>'):

                # Update file_startswith_GT.
                if first_line:
                    file_startswith_GT = True

                # Update dictionary fasta_dict once a sequence has been read.
                else:
                    fasta_dict[head] = seq
            
                # Assign head and header_read.
                head = line
                header_read = True

            # If the line does not start with '>' it is a sequence.
            else:
                
                # Convert line to uppercase.
                line_upper = line.upper()
                
                # Check if the sequence contains only valid characters.
                if sum(line_upper.count(c) for c in valid_chars) != len(line):
                    
                    # Replace invalid characters with 'N' using a string.
                    line_list = list(line_upper)
                    for i, l in enumerate(line_list):
                        if l not in valid_chars:
                            line_list[i] = 'N'
                            invalid_characters_found = True
                    
                    # Convert line_list back to a string.
                    line_upper = ''.join(line_list)
                
                # Assign the variable seq based on if a header was read.
                if header_read:
                    seq = line_upper
                    
                else:
                    seq += line_upper

                # Update header_read following a sequence read.
                header_read = False
                
            # Check if the file starts with '>'.
            if first_line and not file_startswith_GT:
                raise ValueError('FASTA file corrupted. The FASTA file must '
                                 'start with ">".')
                
            # Update first_line to False.
            first_line = False
                
            # Check if there are no more lines to read.
            if not line:
                
                # Add the final head and seq pair to fasta_dict.
                fasta_dict[head] = seq
                
                # Warn user if invalid characters were found.
                if invalid_characters_found:
                    print('Warning: invalid characters found in the input '
                          'file. These have been converted to "N".')
                    
                # Break the while loop.
                break

    # Assign all sequence lengths to the list seq_lengths.
    seq_lengths = []
    for i in fasta_dict.values():
        seq_lengths.append(len(i))
        
    # Check that more than one sequence was read.
    if len(seq_lengths) < 2:
        raise ValueError('The FASTA file does not contain more than 1 '
                         'sequence.')
        
    # Ensure that all sequences are of the same length.
    if len(np.unique(seq_lengths)) > 1:
        raise ValueError('Not all sequences in the FASTA file are of the same '
                         'length.')
        
    # Return the dictionary fasta_dict.
    return fasta_dict


# ---------------------------
# Parameter importer function
# ---------------------------

# Define parameter_importer.
def parameter_importer(path: str = None) -> dict:
        
    '''
    This function creates a dictionary of alignment scoring parameters. The
    user can pass a file path to a text document with custom parameter values,
    in which case those parameters will be used to populate the dictionary
    instead. The format for the parameters file is as follows:
        
        match_score = <value>
        transition = <value>
        transversion = <value>
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
    
        # Import parameters passed to the function.
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
                                  'line must contain "=".')
    
                        # Check for redundant equals signs.
                        if line.count('=') > 1:
                            raise ValueError('The parameter file is '
                                             'misconfigured. Every line must '
                                             'contain no more than 1 "=".')
                            
                        # Assign characters before "=" to param_read.
                        param_read = line[: equals_sign_index].strip()
                        
                        # Check if the parameter has a matching key.
                        if param_read.lower() not in parameters_dict.keys():
                            raise ValueError('The parameter file is '
                                             'misconfigured. Please use the '
                                             'following format:'
                                             '\nmatch_score = <value>'
                                             '\ntransition = <value>'
                                             '\ntransversion = <value>'
                                             '\ngap_penalty = <value>')
                        
                        else:
                            
                            # Try converting the value right of "=" to float.
                            try:
                                value = float(line[equals_sign_index
                                                   + 1:].strip())
                            
                            except ValueError:
                                print('Parameter value error: all parameter '
                                      'values must be numerical.')
                                
                            # Update the dictionary with the new value.
                            parameters_dict[param_read.lower()] = value
                        
                    # Break the while loop.
                    else:
                        break
            
    # Return the dictionary.
    return parameters_dict