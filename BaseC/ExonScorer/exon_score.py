"""
Exon_proclivity algorithm
"""

def exon_score(rna, dna):
    exon_proclivity = {}
    for n, item in enumerate(dna.string):
        exon_proclivity[f"{n}_{item}"] = []
    for item in dna.suffix_index_filtered:
        if item[1][0] == rna.string[0]:
            boolean = has_subsequence(rna.string, item[1])
            if boolean == True:
                positions = subsequence_match_positions(rna.string, item[1])
                if positions is not None:
                    for i, pos in enumerate(positions):
                        abs_pos = item[0] + pos
                        nt = dna.string[abs_pos]
                        weight = len(rna.string)/(positions[-1] - positions[0] + 1)
                        exon_proclivity[f"{abs_pos}_{nt}"].append(weight)
        else:
            exon_proclivity[f"{item[0]}_{item[1][0]}"].append(0)
    for key, values in exon_proclivity.items():
        if values:
            exon_proclivity[key] = sum(values) / len(values)
        else:
            exon_proclivity[key] = 0
    return exon_proclivity

def has_subsequence(query, target):
    q = 0
    for char in target:
        if q < len(query) and char == query[q]:
            q += 1
    return q == len(query)

def subsequence_match_positions(query, target):
    positions = []
    q = 0
    for i, char in enumerate(target):
        if q < len(query) and char == query[q]:
            positions.append(i)      
            q += 1
    if q == len(query):
        return positions
    return None