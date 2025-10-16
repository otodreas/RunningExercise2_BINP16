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
    raise IndexError('Too few arguments passed. You must pass the path to the input file to the program.')
    
# If four or more, raise an error.
elif len(sys.argv) > 4:
    raise IndexError('Too many arguments passed. Please pass no more than three arguments to the program.')

# Check that the first file passed exists.
elif not os.path.exists(sys.argv[1]):
    raise ValueError('The input file does not exist. Ensure that the input file is passed as the first argument!')
    
# If one argument is passed after the program name, assign it to input_path.
elif len(sys.argv) == 2:
    input_path = sys.argv[1]

# If two, assign the first to input_path and ask the user which type of file
# the second argument is referring to.
elif len(sys.argv) == 3:
    input_path = sys.argv[1]
    argtype = input('The identity of the last argument is ambiguous. Are you providing a path to an output file? ("yes"/"y" for output file, "no"/"n" for parameters file): ')
    if argtype == 'yes' or argtype == 'y':
        output_path = sys.argv[2]
    elif argtype == 'no' or argtype == 'n':
        parameters_path = sys.argv[2]
    else:
        raise ValueError('Invalid answer. Please start over.')
    
# If three, assign them to input_path, parameters_path, and output_path.
elif len(sys.argv) == 4:
    input_path, parameters_path, output_path = sys.argv[1:4]
    

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
                    line_list = list(line_upper)
                    for i, l in enumerate(line_list):
                        if l not in valid_chars:
                            line_list[i] = 'N'
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

### 3.1: Import parameters if necessary.
### ------------------------------------

# Run the fasta_importer function on the file passed to the program by the 
# user.
fasta_dict = fasta_importer(input_path)

# Assign scores to variables in a dictionary using default values.
parameters_dict = {'match_score': 1,
                   'transition': -1, # a<->g or c<->t
                   'transversion': -2, # other mutations
                   'gap_penalty': -1}

# If a parameter document is passed to the program, import the preferred
# scores and update the dictionary with their values.
if parameters_path is not None:
    
    # Open the file to read it.
    with open(parameters_path, 'r') as file:
        
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
                    print('The parameter file is misconfigured. Every line must contain an equals sign.')
                    
                # assign the characters before the equals sign to the variable
                # param_read, stripping leading whitespaces.
                param_read = line[: equals_sign_index].strip()
                
                # Check if the parameter has a match in the parameter
                # dictionary.
                if param_read in parameters_dict.keys():
                    
                    # Assign the value following the equals sign to the 
                    # variable value as a float.
                    try:
                        value = float(line[equals_sign_index + 1: ].strip())
                    
                    # Throw an error if the value is not convertable to float.
                    except ValueError:
                        print('Parameter value error: all parameters must be numbers')
                    # Update the dictionary with the new value.
                    parameters_dict[param_read] = value
                    
                # If the parameter is not found, throw an error and interrupt
                # the program. Print the accepted format for the parameter 
                # file.
                else:
                    raise ValueError('Parameter not found error: the parameter file is misconfigured. '
                                     'Please use the following format:\n'
                                     'match_score = <value>\n'
                                     'transition = <value>\n'
                                     'transversion = <value>\n'
                                     'gap_penalty = <value>\n')
        
            # Break the while loop when there are no more lines to be read.
            else:
                break
        
### 3.2: Calculate scores for each pairwise alignment.
### --------------------------------------------------

# Assign an empty list to the variable score_summary.
score_summary = ''

# Loop through the keys in the dictionary, calculating each possible pairwise
# score.
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
            disregarded_positions = 0

            # Zip the sequences to compare alignment at each position and loop through nucleotide by nucleotide.
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


### 3.3: Write output file.
### -----------------------
with open(output_path, 'w') as f:
    f.write(score_summary)