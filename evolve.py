from gene import Gene
import os
import shutil
import json
from pathlib import Path
from tqdm import trange
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from colabfold.batch import get_queries, run
from colabfold.download import default_data_dir
from colabfold.utils import setup_logging
from matplotlib import pyplot as plt

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
        filename = f'{folder}/gene{i}.fasta'
        gene.fasta = filename
        SeqIO.write(record, f'{folder}/gene{i}.fasta', 'fasta')
    
def fold(input_dir, result_dir, logfile, kwargs):
    # setup_logging(Path(result_dir).joinpath(logfile))
    kwargs['queries'], kwargs['is_complex'] = get_queries(input_dir)
    kwargs['result_dir'] = result_dir
    results = run(**kwargs)
    return results

def calculate_fitnesses(genes, results):
    for gene, metric in zip(genes, results['metric']):
        mean_plddt, ptm = metric[0]['mean_plddt'], metric[0]['ptm']
        fitness = (mean_plddt / 100) * 0.5 + ptm * 0.5 + len(gene.get_protein()) * 0.01
        gene.fitness = fitness

def log_evolution(generations, logfile):
    with open(logfile, 'w') as f:
        data = [[gene.__dict__ for gene in gen] for gen in generations]
        json.dump(data, f, indent=4)

def read_evolution(logfile):
    with open(logfile, 'r') as f:
        data = json.load(f)
        generations = [[Gene(**d) for d in gen] for gen in data]
    return generations

def plot_fitnesses(generations):
    fig = plt.figure()
    x = np.arange(len(generations))
    y = [np.mean([gene.fitness for gene in gen]) for gen in generations]
    yerr = [np.std([gene.fitness for gene in gen]) for gen in generations]
    plt.errorbar(x, y, yerr=yerr)
    plt.show()

if __name__ == '__main__':

    # i cant run this from the command line bc idk how localcolabfold conda setup works oopsie
    #parser = argparse.ArgumentParser()
    #parser.add_argument('--n_genes', type=int, default=10)
    #parser.add_argument('--n_bases', type=int, default=30)
    #parser.add_argument('--n_generations', type=int, default=10)
    #parser.add_argument('--proteins_path', type=str, default='proteins')
    #parser.add_argument('--folds_path', type=str, default='folds')
    #parser.add_argument('--logfile', type=str, default='log.txt')
    #args = parser.parse_args()

    n_genes = 3
    n_bases = 50
    n_generations = 3
    proteins_path = 'proteins'
    folds_path = 'folds'
    logfile = 'log.txt'
    evolution_file = 'evolution.json'

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

    clear_dir(proteins_path)
    clear_dir(folds_path)

    for i in trange(n_generations+1):
        # first generation
        if i == 0:
            generation = [Gene.generate(n_bases) for _ in range(n_genes)]
            # filter empty proteins (begins with stop codon)
            generation = [gene for gene in generation if len(gene.get_protein()) > 0]
            save_proteins(generation, f'{proteins_path}/gen{0}')
            results = fold(f'{proteins_path}/gen{0}', f'{folds_path}/gen{0}', logfile, colabfold_params)
            calculate_fitnesses(generation, results)
            generations = [generation]
        # subsequent generations
        else:
            mutations = [gene.random_mutation() for gene in generation]
            # remove empty proteins
            mutations = [gene for gene in mutations if len(gene.get_protein()) > 0]

            # fold mutated proteins and calculate fitnesses as a batch
            save_proteins(mutations, f'{proteins_path}/gen{i}')
            results = fold(f'{proteins_path}/gen{i}', f'{folds_path}/gen{i}', logfile, colabfold_params)
            calculate_fitnesses(mutations, results)
            generation += mutations

            # stochastic selection
            fitnesses = [gene.fitness for gene in generation]
            generation = np.random.choice(generation, size=n_genes, replace=False, p=np.divide(fitnesses, sum(fitnesses))).tolist()
            generations.append(generation)

    log_evolution(generations, evolution_file)