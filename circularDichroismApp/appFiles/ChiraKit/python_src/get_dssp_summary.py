import os
import pandas as pd

def letter_to_secondary_structure(letter):

    letter_to_secondary_structure = {
    'H':'Alpha-helix',
    'G':'Alpha-helix',
    'E':'Beta-sheet',
    'T':'Turns',
    '.':'Other',
    'I':'Other',
    'B':'Other',
    'S':'Other',
    'C':'Other',
    'P':'Other'
    }

    return letter_to_secondary_structure[letter]

def split_alpha_helix(lst):

    """
    Given a list of strings, if there are subsequences of consecutive elements equal to 'Alpha-helix',
    replace the first and last two ocurrences with 'Alpha-d' and the middle ones with 'Alpha-r'
    If the sequence length of 'Alpha-helix's is less than 5 elements, replace all of them with 'Alpha-d'

    Return the modified list 
    """

    count_a = 0

    lst = [ 'Alpha-d' if l == 'Alpha-helix' else l for l in lst ]

    lst.append('.')

    for i in range(len(lst)-1):

        if lst[i] == 'Alpha-d':
            count_a += 1

            if count_a > 4 and lst[i+1] != 'Alpha-d':
                # Replace the first two and last two 'Alpha-d's with 'Alpha-r's
                
                init = i-count_a+3
                end  = i-1

                lst[init:end] = ['Alpha-r'] * (end-init) 
                
                count_a = 0

        else:
            
            count_a = 0

    return lst[:-1]

def split_beta_sheet(lst):
    
    """
    Given a list of strings, if there are subsequences of consecutive elements equal to 'Beta-sheet',
    replace the first and last ocurrence with 'Beta-d' and the middle ones with 'Beta-r'
    If the sequence length of 'Beta-sheet's is less than 3 elements, replace all of them with 'Beta-d'

    Return the modified list 
    """

    count_b = 0

    lst = [ 'Beta-d' if l == 'Beta-sheet' else l for l in lst ]

    lst.append('.')

    for i in range(len(lst)-1):

        if lst[i] == 'Beta-d':
            count_b += 1

            if count_b > 2 and lst[i+1] != 'Beta-d':
                # Replace the first two and last two 'Alpha-d's with 'Alpha-r's
                
                init = i-count_b+2
                end  = i

                lst[init:end] = ['Beta-r'] * (end-init) 
                
                count_b = 0

        else:
            
            count_b = 0

    return lst[:-1]

def read_secondary_structure(dssp_generated_file):

    """
    Open the DSSP generated file and compute the six desired secondary structure components:

        Alpha-regular   Alpha-distorted     Beta-regular    Beta-distorted  Turns   Other

    """

    with open(dssp_generated_file,'r') as f:
        ls = f.read().splitlines()

        dsspMetadaLine    = [i for i,l in enumerate(ls) if '_dssp_struct_summary.secondary_structure' in l][0]
        dsspDataFirstLine = [i for i,l in enumerate(ls[dsspMetadaLine:]) if len(l.split()) > 5][0]
        dsspDataLastLine  = [i for i,l in enumerate(ls[(dsspMetadaLine + dsspDataFirstLine):]) if len(l.split()) < 5][0]

        ls = ls[(dsspMetadaLine + dsspDataFirstLine):(dsspMetadaLine + dsspDataFirstLine + dsspDataLastLine)]

        alpha_regular   = []
        alpha_distorted = []
        beta_regular    = []
        beta_distorted  = []
        turns           = []
        other           = []

        secondary_structure_elements  = [letter_to_secondary_structure(l.split()[4]) for l in ls]

        secondary_structure_elements = split_alpha_helix(secondary_structure_elements)
        secondary_structure_elements = split_beta_sheet(secondary_structure_elements)

        # Create a Series and use value_counts to count unique values
        value_counts_series = pd.Series(secondary_structure_elements).value_counts()

        # Calculate percentages
        percentages = round(value_counts_series / len(secondary_structure_elements) * 100,1)

        # Create a DataFrame from the Series and percentages
        df = pd.DataFrame({'Sec_str_component': value_counts_series.index, 'Count': value_counts_series.values, 'Percentage': percentages.values})

        return df

def run_dssp_workflow(pdb_file):

    """

    Run the DSSP algorithm for secondary structure calculation using as input a PDB file. 
    Requires the script 'mkdssp-4.4.0-linux-x64' to be available from the console (e.g., in the $PATH folders)
    
    The input file format can be either PDB or mmCIF
    The output is a dataframe with six secondary structure components:

        Alpha-regular
        Alpha-distorted
        Beta-regular
        Beta-distorted
        Turns
        Other

    The total alpha-helix content corresponds to the DSSP fraction 'H' and 'G'.
    The total beta-sheet  content corresponds to the DSSP fraction 'E'.
    Turns are derived from the DSSP fraction 'T'. 
    The first and last two residues of an alpha-helix are considered 'distorted'. 
    The first and last     residue  of a  beta-sheet  are considered 'distorted'.

    """

    df = None

    try:

        os.system('[ -e dssp_input_file ] && rm dssp_input_file')
        os.system('mkdssp-4.4.0-linux-x64 ' + pdb_file + ' > dssp_input_file')
        df = read_secondary_structure('dssp_input_file')

        # Define the custom sort order
        custom_order = ['Alpha-r','Alpha-d', 'Beta-r', 'Beta-d', 'Turns', 'Other']

        # Convert the 'Sec_str_component' column to a Categorical type with custom sort order
        df['Sec_str_component'] = pd.Categorical(df['Sec_str_component'], categories=custom_order, ordered=True)

        # Sort the DataFrame based on the 'Name' column
        df = df.sort_values(by='Sec_str_component')

    except:

        pass

    return df