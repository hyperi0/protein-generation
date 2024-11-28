import random
import os
import shutil
import glob

from gene import Gene

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from colabfold.batch import get_queries, run
from colabfold.download import default_data_dir
from colabfold.utils import setup_logging
from pathlib import Path

def clear_dir(filepath):
    for root, dirs, files in os.walk(filepath):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))

def save_protein(protein, ident, desc, filename):
    record = SeqRecord(protein, id=ident, description=desc)
    SeqIO.write(record, filename, 'fasta')

def first_generation(n_genes, n_bases, filepath='evolution'):
    try:
        os.mkdir('evolution/gen0/')
    except FileExistsError: # clear directory contents
        files =- glob.glob
    genes = [Gene(n_bases) for _ in range(n_genes)]
    for i, gene in enumerate(genes):
        ident = str(i)
        save_protein(gene.get_protein(), ident, str(gene), f'{filepath}/gen0/{ident}.fasta')
    return genes

def generate_mutations(genes, filepath, n_generation):
    # make folder or clear folder from previous run
    gen_path = f'{filepath}/gen{n_generation}'
    try:
        os.mkdir(gen_path)
    except FileExistsError: 
        files = glob.glob(gen_path)
        for f in files:
            os.remove(f)
    # mutate genes and save new proteins
    new_genes = [gene.clone() for gene in genes]
    for i, gene in enumerate(new_genes):
        gene.random_mutation()
        ident = str(i)
        save_protein(gene.get_protein(), ident, str(gene), f'{gen_path}/{ident}.fasta')
    return new_genes
    
def fold(input_dir, result_dir, kwargs):
    setup_logging(Path(result_dir).joinpath("log.txt"))
    kwargs['queries'], kwargs['is_complex'] = get_queries(input_dir)
    kwargs['input_dir'], kwargs['result_dir'] = input_dir, result_dir
    results = run(**kwargs)
    return results

def calculate_fitnesses(results):
    fitnesses = 0
    for metric in results['metric']:
        mean_plddt, ptm = metric[0]['mean_plddt'], metric[0]['ptm']
        fitness = (mean_plddt / 100) * 0.5 + ptm * 0.5
        fitnesses.append(fitness)
    return fitnesses

def select(genes, fitnesses, n, method='weak'):
    if method == 'weak':
        return random.choices(genes, weights=fitnesses, k=n)
    else:
        print('not implemented hehe :3')
        return None
    
# TODO: setup command line argument parsing
n_genes = 10 # start with 10 genes, add 10 mutations then select 10
n_bases = 30 # number of bases in initial population
proteins_filepath = 'proteins'
folds_filepath = 'folds'

colabfold_params = {
    'msa_mode': "single_sequence", #@param ["MMseqs2 (UniRef+Environmental)", "MMseqs2 (UniRef only)","single_sequence","custom"]
    'num_models': 1, #@param [1,2,3,4,5] {type:"raw"}
    'num_recycles': 3, #@param [1,3,6,12,24,48] {type:"raw"}
    'stop_at_score': 100, #@param {type:"string"}
    'use_custom_msa': False,
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

    # initialize random first generation of genes
    # save and fold proteins 
    # calculate fitnesses
    # for each generation:
    #   mutate genes
    #   save and fold proteins from mutated genes
    #   calculate fitnesses for mutated proteins
    #   select genes from mutations and old genes



    for i in range(n_generations):
        save_proteins(genes, proteins_filepath)
        results = fold(proteins_filepath)
        fitnesses = calculate_fitnesses(results)
        mutations = generate_mutations(genes)
        


    # actual code
    generations = []
    genes = first_generation(n_genes, n_bases, proteins_filepath)

    results = fold(f'{genes_filepath}/gen0', 'folds', colabfold_params)
    fitnesses = calculate_fitnesses(results)
    next_gen = select(genes, fitnesses, n_genes)

    mutations = generate_mutations(next_gen)
    results = fold(f'{genes_filepath}/gen1', 'folds', colabfold_params)
