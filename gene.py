from Bio.Seq import Seq
from copy import deepcopy
import random

class Gene():
    def __init__(self, length):
        self.nts = [random.choice(('A', 'T', 'C', 'G')) for _ in range(length)]
        self.generation = 0
        self.mutation = None

    def __str__(self):
        return f'nts: {"".join(self.nts)}\n generation: {self.generation}\n mutation: {self.mutation}'

    def get_protein(self):
        return Seq(''.join(self.nts)).translate(stop_symbol="")

    def clone(self):
        new_gene = deepcopy(self)
        new_gene.generation += 1
        return new_gene

    def point_insertion(self):
        base = random.choice(('A', 'T', 'C', 'G'))
        pos = random.randint(0, len(self.nts))
        self.nts.insert(pos, base)
        self.mutation = 'point_insertion'

    def point_deletion(self):
        pos = random.randint(0, len(self.nts)-1)
        del self.nts[pos]
        self.mutation = 'point_deletion'

    def point_mutation(self):
        pos = random.randint(0, len(self.nts)-1)
        base = random.choice(list({'A', 'T', 'C', 'G'} - {self.nts[pos]}))
        self.nts[pos] = base
        self.mutation = 'point_mutation'
    
    def partial_deletion(self):
        start = random.randint(0, len(self.nts)-1)
        end = random.randint(start+3, len(self.nts))
        self.nts = self.nts[start:end]
        self.mutation = 'partial_deletion'
    
    def random_insertion(self):
        pos = random.randint(0, len(self.nts))
        seq_len = random.randint(1, len(self.nts))
        new_seq = [random.choice(('A', 'T', 'C', 'G')) for _ in range(seq_len)]
        self.nts[pos:pos] = new_seq
        self.mutation = 'random_insertion'

    def partial_duplication(self):
        pos = random.randint(0, len(self.nts)-1)
        seq_len = random.randint(1, len(self.nts))
        new_seq = self.nts[pos:pos+seq_len]
        self.nts[pos:pos] = new_seq
        self.mutation = 'partial_duplication'

    def circular_permutation(self):
        pos = random.randint(1, len(self.nts)-1)
        self.nts = self.nts[pos:] + self.nts[:pos]
        self.mutation = 'circular_permutation'

    def full_duplication(self):
        self.nts = self.nts + self.nts
        self.mutation = 'full_duplication'

    def random_mutation(self):
        mutations = [
            self.point_insertion,
            self.point_mutation,
            self.random_insertion,
            self.partial_duplication,
            self.circular_permutation,
            self.full_duplication,
        ]
        if len(self.nts) > 3:
            mutations += [self.point_deletion, self.partial_deletion]
        mutation = random.choice(mutations)
        mutation()

# test all mutations
if __name__ == '__main__':
    genes = [Gene(10) for _ in range(8)]
    all_mutations = [
        'point_insertion',
        'point_mutation',
        'random_insertion',
        'partial_duplication',
        'circular_permutation',
        'full_duplication',
        'point_deletion',
        'partial_deletion',
    ]
    new_genes = []
    for gene, mutation in zip(genes, all_mutations):
        mutation = getattr(Gene, mutation)
        new_gene = gene.clone()
        mutation(new_gene)
        print(new_gene)
        print(len(new_gene.nts))