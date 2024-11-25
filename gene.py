from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import json
import os

class Gene():
    def __init__(self, nts, generation, parent, mutation):
        self.nts = nts
        self.generation = generation
        self.parent = parent
        self.mutation = mutation
        
    def get_protein(self):
        return Seq(self.nts).translate(stop_symbol="")
    
    def point_insertion(self):
        base = random.choice(('A', 'T', 'C', 'G'))
        pos = random.randint(0, len(self.nts))
        self.nts.insert(pos, base)

    def point_deletion(self):
        if len(self.nts) > 1:
            pos = random.randint(0, len(self.nts)-1)
            self.nts = self.nts[:pos] + self.nts[pos+1:]

    def point_mutation(self):
        pos = random.randint(0, len(self.nts)-1)
        base = random.choice(list({'A', 'T', 'C', 'G'} - {self.nts[pos]}))
        self.nts[pos] = base
    
    def partial_deletion(self):
        if len(self.nts) > 1:
            start = random.randint(0, len(self.nts)-1)
            end = random.randint(start+1, len(self.nts))
            self.nts = self.nts[start:end]
    
    def random_insertion(self):
        pos = random.randint(0, len(self.nts))
        seq_len = random.randint(1, len(self.nts))
        new_seq = ''.join([random.choice(('A', 'T', 'C', 'G')) for _ in range(seq_len)])
        self.nts = self.nts[:pos] + new_seq + self.nts[pos:]

    def partial_duplication(self):
        pos = random.randint(0, len(self.nts)-1)
        seq_len = random.randint(1, len(self.nts))
        new_seq = self.nts[pos:pos+seq_len]
        self.nts = self.nts[:pos] + new_seq + self.nts[pos:]

    def circular_permutation(self):
        pos = random.randint(1, len(self.nts)-1)
        self.nts = self.nts[pos:]  + self.nts[:pos]

    def full_duplication(self):
        self.nts = self.nts = self.nts

    

class Evolution():
    def __init__(self):
        all_genes = []

    def generate_gene(self, length):
        nts = ''.join([random.choice(('A', 'T', 'C', 'G')) for _ in range(length)])
        return Gene(nts, 0, None, None)
    
    def 