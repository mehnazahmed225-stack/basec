"""
FASTA PARSER NOTES:

When using the Fasta File Parser, the parser returns current ID as the keys and the name of the fasta entry as the description which is index zero of the list. 
The index 1 of the list is the actual sequence in the entry. 

"""

def fasta_parser(filepath):
    sequences = {}
    with open(filepath, "r") as f:
        seq = 1
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_description = line[1:]
                current_id = f"SEQ_{seq}"
                seq+=1
                sequences[current_id] = [current_description, ""]
            elif current_id is not None:
                sequences[current_id][1] += line

    return sequences
