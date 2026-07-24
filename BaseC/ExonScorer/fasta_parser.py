"""
FASTA PARSER NOTES:

When using the Fasta File Parser, make sure that there is only one entry. This is because it can process only one gene/mrna at a time. 
The fasta parser only returns the nucleotides of the entry, not the descriptor/name

"""

def fasta_parser(filepath):
    sequences = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:]
                sequences[current_id] = ""
            elif current_id is not None:
                sequences[current_id] += line

    return sequences