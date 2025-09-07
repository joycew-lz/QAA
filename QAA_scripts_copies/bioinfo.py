#!/usr/bin/env python

# Author: Joyce Wang <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score.'''
    return ord(letter) - 33
    pass

def qual_score(phred_score: str) -> float:
    '''Calculates the average phred quality score of a phred string.
    For each letter in a string, convert to a phred score.
    Store this phred score to the variable total_score, and
    divide the total_score with the number of letters looped over
    in the string.'''
    total_score = 0
    for letter in phred_score:
        total_score += convert_phred(letter)
    average_score = total_score/len(phred_score)
    return average_score
    pass

DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)
    pass

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)
    pass

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list as a float.'''
    length = len(lst)
    # get the middle index
    mid = length // 2

    if length % 2 == 0: # if even number of elements, average the two middle values
         return float((lst[mid - 1] + lst[mid])/2)
    else: # if odd number of elements, the middle index's element is the median
        return float(lst[mid])
    pass

def oneline_fasta(input_file: str, output_file: str):
    '''
    Takes a fasta file and writes a new file where each sequence is on a single line
    following its  ">" header.
    '''
    with open(input_file, 'r') as i_fh:
        with open(output_file, 'w') as o_fh:
            sequence_lines = ''
            for line in i_fh:
                line = line.strip('\n')
                if ">" in line:
                    if sequence_lines == '':
                        o_fh.write(line + '\n')
                    else:
                        o_fh.write(sequence_lines + '\n')
                        sequence_lines = ''
                        o_fh.write(line + '\n')
                else:
                    sequence_lines += line
            o_fh.write(sequence_lines + '\n')
    pass

if __name__ == "__main__":
    # These tests are run when you execute this file directly (./bioinfo.py in the terminal) (instead of importing it)
    # Remember, assert statements are meant to be true, and if not, it moves on to the next line of code (displays the string after)

    # tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job.")

    # tests for qual_score
    assert qual_score("III") == 40
    assert qual_score("JG") == 39.5
    assert qual_score("#C") == 18
    print("You calcluated the correct average phred score.")

    # tests for validate_base_seq
    assert validate_base_seq("GCATTACC") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AGCUAGCU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("I love bunnies!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("I love bunnies!", True) == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("gctataca"), "Validate base seq does not work on lowercase DNA"
    assert validate_base_seq("ugacu", True), "Validate base seq does not work on lowercase RNA"
    print("Passed DNA and RNA tests.")

    # tests for gc_content
    assert gc_content("GCGCGC") == 1
    assert gc_content("TATAAATT") == 0
    assert gc_content("GCTAGCTAGCTA") == 0.5
    print("GC content successfully calculated.")

    # tests for calc_median
    assert calc_median([4,100,200]) == 100, "calc_median function does not work for odd length list"
    assert calc_median([8,9]) == 8.5, "calc_median function does not work for even length list"
    assert calc_median([7,8,9,1000,40901]) == 9
    assert calc_median([2,3,6,9,10,15,17,18]) == 9.5
    print("Median successfully calculated.")

    # no need to test oneline_fasta

    # make sure these unit tests work by executing this file!

    # git add, commit, and push to see if it passes Leslie's autograders.
