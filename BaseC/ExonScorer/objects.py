"""
These are the objects to be used in the ExonFinder algorithm
"""

class dna:
    def __init__(self, string):
        self.string = string
        self.suffix_index_filtered = []
        self.suffix_index = [
            (i, string[i:])
            for i in range(len(string))
        ]
        self.rc_string = self.reverse_complement()
        self.rc_suffix_index = [
            (i, self.rc_string[i:])
            for i in range(len(self.rc_string))
        ]
    def filter(self, rna):
        self.suffix_index_filtered = []
        for item in self.suffix_index:
            if len(item[1]) >= len(rna.string):
                self.suffix_index_filtered.append(item)
        return self
    def reverse_complement(self):
        complement = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
        }
        return "".join(
            complement[base]
            for base in reversed(self.string)
        )
    def reverse_suffix_index(self):
        rc = self.reverse_complement()
        return [
            (i, rc[i:])
            for i in range(len(rc))
        ]

class mrna:
    def __init__(self, string):
        self.string = string.upper().replace("U", "T")
