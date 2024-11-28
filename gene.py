from Bio.Seq import Seq
from copy import deepcopy
import random

class Gene():
    def __init__(self, length):
        self.bases = [random.choice(('A', 'T', 'C', 'G')) for _ in range(length)]
        self.mutations = []

    def __str__(self):
        return f'bases: {"".join(self.bases)}\n mutations: {self.mutations}'

    def get_protein(self):
        return Seq(''.join(self.bases)).translate(stop_symbol="")

    def clone(self):
        new_gene = deepcopy(self)
        return new_gene

    def point_insertion(self):
        base = random.choice(('A', 'T', 'C', 'G'))
        pos = random.randint(0, len(self.bases))
        self.bases.insert(pos, base)
        self.mutations.append('point_insertion')

    def point_deletion(self):
        pos = random.randint(0, len(self.bases)-1)
        del self.bases[pos]
        self.mutations.append('point_deletion')

    def point_mutation(self):
        pos = random.randint(0, len(self.bases)-1)
        base = random.choice(list({'A', 'T', 'C', 'G'} - {self.bases[pos]}))
        self.bases[pos] = base
        self.mutations.append('point_mutation')
    
    def partial_deletion(self):
        start = random.randint(0, len(self.bases)-3)
        end = random.randint(start+3, len(self.bases))
        self.bases = self.bases[start:end]
        self.mutations.append('partial_deletion')
    
    def random_insertion(self):
        pos = random.randint(0, len(self.bases))
        seq_len = random.randint(1, len(self.bases))
        new_seq = [random.choice(('A', 'T', 'C', 'G')) for _ in range(seq_len)]
        self.bases[pos:pos] = new_seq
        self.mutations.append('random_insertion')

    def partial_duplication(self):
        pos = random.randint(0, len(self.bases)-1)
        seq_len = random.randint(1, len(self.bases))
        new_seq = self.bases[pos:pos+seq_len]
        self.bases[pos:pos] = new_seq
        self.mutations.append('partial_duplication')

    def circular_permutation(self):
        pos = random.randint(1, len(self.bases)-1)
        self.bases = self.bases[pos:] + self.bases[:pos]
        self.mutations.append('circular_permutation')

    def full_duplication(self):
        self.bases = self.bases + self.bases
        self.mutations.append('full_duplication')

    def random_mutation(self):
        new_gene = self.clone()
        mutations = [
            new_gene.point_insertion,
            new_gene.point_mutation,
            new_gene.random_insertion,
            new_gene.partial_duplication,
            new_gene.circular_permutation,
            new_gene.full_duplication,
        ]
        if len(self.bases) > 3:
            mutations += [new_gene.point_deletion, new_gene.partial_deletion]
        mutation = random.choice(mutations)
        mutation()
        return new_gene

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
        print(len(new_gene.bases))