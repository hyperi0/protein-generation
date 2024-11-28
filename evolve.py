import random
import os
import shutil
import glob
from pathlib import Path
from tqdm import trange

from gene import Gene

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from colabfold.batch import get_queries, run
from colabfold.download import default_data_dir
from colabfold.utils import setup_logging

def clear_dir(filepath):
    for root, dirs, files in os.walk(filepath):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))

def save_proteins(genes, folder):
    os.mkdir(folder)
    for i, gene in enumerate(genes):
        protein = gene.get_protein()
        record = SeqRecord(protein, id=str(i), description=str(gene))
        SeqIO.write(record, f'{folder}/gene{i}.fasta', 'fasta')
    
def fold(input_dir, result_dir, kwargs):
    setup_logging(Path(result_dir).joinpath("log.txt"))
    kwargs['queries'], kwargs['is_complex'] = get_queries(input_dir)
    kwargs['result_dir'] = result_dir
    results = run(**kwargs)
    return results

def calculate_fitnesses(results):
    fitnesses = []
    for metric in results['metric']:
        mean_plddt, ptm = metric[0]['mean_plddt'], metric[0]['ptm']
        fitness = (mean_plddt / 100) * 0.5 + ptm * 0.5
        fitnesses.append(fitness)
    return fitnesses

def log_evolution(generations, logfile):
    with open(logfile, 'w') as f:
        for i, generation in enumerate(generations):
            f.write(f'-------------Generation {i}---------------\n')
            for gene, fitness in generation.items():
                f.write(str(gene) + ": " + str(fitness))
            f.write('\n')

    
# TODO: setup command line argument parsing
n_genes = 10
n_bases = 30
n_generations = 10
proteins_filepath = 'proteins'
folds_filepath = 'folds'
logfile = 'log.txt'

colabfold_params = {
    'msa_mode': "single_sequence", #@param ["MMseqs2 (UniRef+Environmental)", "MMseqs2 (UniRef only)","single_sequence","custom"]
    'num_models': 1, #@param [1,2,3,4,5] {type:"raw"}
    'num_recycles': 3, #@param [1,3,6,12,24,48] {type:"raw"}
    'stop_at_score': 100, #@param {type:"string"}
    'num_relax': 0, #@param [0, 1, 5] {type:"raw"}
    'use_amber': False,
    'relax_max_iterations': 200, #@param [0,200,2000] {type:"raw"}
    'use_templates': False, #@param {type:"boolean"}
    'keep_existing_results': True, #@param {type:"boolean"}
    'zip_results': False, #@param {type:"boolean"}
    'model_type': "auto",
    'model_order': [1, 2, 3, 4, 5],
    'data_dir': default_data_dir,
    'keep_existing_results': True,
    'rank_by': "auto",
    'pair_mode': "unpaired+paired",
    'user_agent': "colabfold/google-colab-batch",
}


if __name__ == '__main__':
    clear_dir(proteins_filepath)
    clear_dir(folds_filepath)

    generations = []
    generation = {}
    for i in trange(n_generations):
        if i == 0:
            mutations = [Gene(n_bases) for _ in range(n_genes * 2)]
        else:
            mutations = [gene.random_mutation() for gene in generation.keys()]
        save_proteins(mutations, f'{proteins_filepath}/gen{i}')
        results = fold(f'{proteins_filepath}/gen{i}', f'{folds_filepath}/gen{i}', colabfold_params)
        mutation_fitnesses = calculate_fitnesses(results)
        genes = list(generation.keys()) + mutations
        fitnesses = list(generation.values()) + mutation_fitnesses
        selection = np.random.choice(genes, size=n_genes, replace=False, p=np.divide(fitnesses, sum(fitnesses)))
        generation = {gene: fitness for gene, fitness in zip(genes, fitnesses) if gene in selection}
        generations.append(generation)

    log_evolution(generations, logfile)