#Parsing fasta files to yield a dictionary of sequences with the IDs as the keys
def parser_fastq(filename):
    sequences = {}
    with open(filename, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            if not header.startswith('@'):
                raise ValueError(f"Invalid FASTQ header: {header}")
            if plus != '+':
                raise ValueError("Invalid FASTQ format")
            if len(seq) != len(qual):
                raise ValueError("Sequence and quality lengths differ")
            sequences[header[1:]] = (seq, qual)
    return sequences

#Parsing fasta files
def parser_fasta(filename):
    sequences = {}
    with open(filename, "r") as f:
        current_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)
    return sequences

#Kmer objects that can be generated from a string with n steps and offset accordingly
class Kmer(object):

    def __init__(self, id, sequence, k, quality=None, step=1, offset=0):
        self.id = id
        self.sequence = sequence
        self.quality = quality
        self.offset = offset
        self.k = k
        self.step = step
        self.index = {}
        for i in range(offset, len(sequence) - k + 1, step):
            kmer = sequence[i:i+k]
            if quality is not None:
                qmer = quality[i:i+k]
                entry = (i, qmer)
            else:
                entry = (i,)
            if kmer not in self.index:
                self.index[kmer] = []
            self.index[kmer].append(entry)
    #This is so that it doesn't just give you a text wrapper when trying to access the object
    def __repr__(self):
        return (
            f"Kmer id='{self.id}', "
            f"k={self.k}, "
            f"length={len(self.sequence)}, "
            f"unique_kmers={len(self.index)})"
        )
    #Gives you the number of kmers obtained
    def __len__(self):
        return len(self.index)
    #So this enables you to check if a specific kmer is in the dictionary
    def __contains__(self, p):
        return p in self.index
    def __str__(self):
        return str(self.sequence)
    #Shows you the frequency of the kmers(hashtable)
    def frequencies(self):
        return {
            kmer: len(pos)
            for kmer, pos in self.index.items()
        }
    #Unfortunately Binary search wasn't used because python's dictionary lookup is more efficient
    def query(self, p):
        return self.index.get(p, [])
    #This returns the count of positions with a certain kmer
    def count(self, p):
        return len(self.index.get(p, []))
    #This returns a kmer with its positions
    def positions(self, p):
        return self.index.get(p, [])
    #This returns the Kmers used. len(kmers) returns the number of unique kmers. 
    @property
    def n_kmers(self):
        return sum(len(v) for v in self.index.values())

#Making a dat file for storage with a kmer
def kmer_to_dat(list_kmers, filename):
    if isinstance(list_kmers, Kmer):
        list_kmers = [list_kmers]
    with open(f"{filename}.dat", "w") as f:
        for obj in list_kmers:
            f.write("BEGIN_KMER\n")
            f.write(f"id={obj.id}\n")
            f.write(f"Sequence={obj.sequence}\n")
            f.write(f"k={obj.k}\n")
            f.write(f"Quality={obj.quality}\n")
            f.write(f"step={obj.step}\n")
            f.write(f"offset={obj.offset}\n")
            for kmer, positions in sorted(obj.index.items()):
                f.write(f"{kmer}: {positions}\n")
            f.write("END_KMER\n\n")

#To parse a kmer.dat file and make a list of kmer objects
def dat_to_kmer(dat_file):
    kmers = []
    with open(f"{dat_file}.dat", "r") as f:
        lines = [line.strip() for line in f]
    i = 0
    while i < len(lines):
        if lines[i] == "BEGIN_KMER":
            id_ = lines[i + 1].split("=", 1)[1]
            sequence = lines[i + 2].split("=", 1)[1]
            k = int(lines[i + 3].split("=", 1)[1])
            quality = lines[i + 4].split("=", 1)[1]
            step = int(lines[i + 5].split("=", 1)[1])
            offset = int(lines[i + 6].split("=", 1)[1])
            obj = Kmer(
                id=id_,
                sequence=sequence,
                k=k,
                quality=None if quality == "None" else quality,
                step=step,
                offset=offset
            )
            kmers.append(obj)
            while i < len(lines) and lines[i] != "END_KMER":
                i += 1
        i += 1
    # Return a single Kmer if only one was found
    if len(kmers) == 1:
        return kmers[0]
    return kmers

def naive_algorithm(pattern, text):
    occurrences = []
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern:
            occurrences.append(i)
    return occurrences

def naive_hamming(p,t,maxDistance):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        nmm = 0
        for j in range(len(p)):
            if t[i+j] != p[j]:
                nmm +=1
                if nmm > maxDistance:
                    break
        if nmm <= maxDistance:
            occurrences.append(i)
    return occurrences

#Suffix Tree

#Suffix index

#FM index

#Boyer-Moore Matching/Alignment

#ALignment(perfect match)

#Subsequence search

#Data cleaning tools

#Alignment(imperfect)

#Edit distances

#Hamming distances