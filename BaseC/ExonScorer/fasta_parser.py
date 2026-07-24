"""
FASTA PARSER NOTES:

The Fasta Parser returns a dictionary from the fasta file. I'll need to figure out how to make it easier to identify the different sequences in a Fasta File

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
