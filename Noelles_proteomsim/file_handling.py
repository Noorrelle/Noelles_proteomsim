def read_fasta(filepath):
    """
    this functions turns Faster Files into dictiornaries with the key as the ID and the value as the Seqeunce of the protein.

    Parameters
    filepath : str
        your path to the file
    Returns
    -------
    dict : 
        Dictionary with protein IDs as keys and aminoacid sequence as values
    """

    fasta_path = filepath
    protein_map = {} 
    current_id = None 
    current_sequence = [] 

    with open(fasta_path, 'r') as fasta_handle: #opens fasta file in read mode and names it "fasta_handle"
        for line in fasta_handle: #reads line by line
            stripped = line.strip() #removes whitespace from beginning and end of lines
            if not stripped: #checks if line is empty
                continue #skips empty lines
            if stripped.startswith('>'): #selects header lines
                if current_id is not None: #looks whether an ID is there
                    protein_map[current_id] = ''.join(current_sequence) #assign each ID it's sequence fragments joined togethere, starts with whitespace
                    current_sequence = [] #empties current sequence?
                current_id = stripped.split('|')[1] #extracts protein ID from header by spliting text fragmenrs by the character "|" and taking the second element
            else:
                current_sequence.append(stripped) #if it is not a header current sequence is elongated

    if current_id is not None: #saves the last protein as there is no next header anymore
        protein_map[current_id] = ''.join(current_sequence)

    return protein_map #displays dcitionary
