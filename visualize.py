from pymol import cmd
from PIL import Image
import os
import shutil

def clear_dir(filepath):
    for root, dirs, files in os.walk(filepath):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))

def visualize_proteins(images_path, filename, n_generations, n_proteins, width, height):

    clear_dir(images_path)

    protein_width = int(width / n_proteins)
    protein_height = int(height / n_generations)

    for gen in range(n_generations):
        os.mkdir(f'images/gen{gen}')
        pdb_files = [f for f in os.listdir(f'folds/gen{gen}') if f.endswith('pdb')]
        for f in pdb_files:
            cmd.load(f'folds/gen{gen}/' + f)
            cmd.show('cartoon', 'all')
            cmd.png(f'images/gen{gen}/' + f, width=protein_width, height=protein_height)
            cmd.delete('all')

    grid_image = Image.new('RGBA', (width, height), (255, 255, 255, 255))

    for gen in range(n_generations):
        image_files = [f'images/gen{gen}/'+f for f in os.listdir(f'images/gen{gen}')]
        for i, image_file in enumerate(image_files):
            image = Image.open(image_file)
            x = i * protein_width
            y = gen * protein_height
            grid_image.paste(image, (x, y))

    grid_image.save(f'{images_path}/{filename}.png')

if __name__ == '__main__':

    images_path = 'images'
    filename = 'evolution'
    n_generations = len(os.listdir('proteins'))
    n_proteins = max([len(os.listdir(f'proteins/{dir}')) for dir in os.listdir('proteins')])
    width = 1000
    height = 1000

    visualize_proteins(images_path, filename, n_generations, n_proteins, width, height)