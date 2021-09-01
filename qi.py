__author__ = 'florian'
__Filename__ = 'qi'
__Creationdate__ = '22/07/2021'

from cross_validation import CrossValidation as cv
path_file = ''

file_atlas_reference = 'd_prior_360_cell_none'

fichier_position = "spatial.xlsx"
fichier_dataset = 'df_transcriptomique_3000_genes'

time_dataset = 170
# %%  Parameters of novosparc


num_neighbors_s = num_neighbors_t = 20
alpha_linear = 0.5
epsilon = 1e-2

# %%


import novosparc
import time
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import altair as alt
from scipy.spatial.distance import cdist, squareform, pdist
from scipy.stats import ks_2samp
from scipy.stats import pearsonr
import io
import traduction_gene_marqueur as tgm


def position_cell_excel(t, file_name, sheet='Feuil1'):
    df_space = pd.read_excel(io=file_name, sheet_name=sheet)
    df_space = df_space.dropna()
    position = df_space[df_space['time'] == t]
    return position


def gene_cell(t, file_name_gene, file_name_position, sheet_position='Feuil1'):
    position = position_cell_excel(t, file_name_position, sheet_position)
    df = pd.read_csv(file_name_gene)
    df = df.rename(columns={'name': 'cell_name'})
    name_cell = position['cell_name']
    gene_mat = pd.merge(df, name_cell)
    return gene_mat


def embedding(dataset, color, title=None, size_x=None, size_y=None,
              tit_size=15):
    """
    Plots fields (color) of Scanpy AnnData object on spatial coordinates
    dataset -- Scanpy AnnData with 'spatial' matrix in obsm containing the spatial coordinates of the tissue
    color -- a list of fields - gene names or columns from obs to use for color
    """
    title = color if title is None else title
    ncolor = len(color)
    per_row = 3
    per_row = ncolor if ncolor < per_row else per_row
    nrows = int(np.ceil(ncolor / per_row))
    size_x = 5 * per_row if size_x is None else size_x
    size_y = 3 * nrows if size_y is None else size_y
    size_z = 2 * nrows if size_y is None else size_y
    fig = plt.figure(figsize=plt.figaspect(nrows))
    xy = dataset.obsm['spatial']
    x = xy[:, 0]
    y = xy[:, 1] if xy.shape[1] > 1 else np.ones_like(x)
    z = xy[:, 2]

    for i, g in enumerate(color):
        if g in dataset.var_names:
            values = dataset[:, g].X
        elif g in dataset.obs.columns:
            values = dataset.obs[g]
        else:
            continue
        ax = fig.add_subplot(1, 2, i + 1, projection='3d')
        ax.scatter(x, y, z, c=np.array(values), linewidth=0.5)
        ax.set_title(title[i], size=tit_size)

    plt.show()
    plt.tight_layout()


# %% load data


dataset = gene_cell(time_dataset, fichier_dataset, fichier_position)

cell_name_dataset = list(dataset['cell_name'])
np.save('cell_name_dataset', cell_name_dataset)
del dataset['Unnamed: 0']
del dataset['cell_name']
# %%
gene_names = dataset.columns.values

num_cells, num_genes = dataset.shape

dataset.to_csv('modif_transcriptomie.txt', sep=' ', index=False)
data_path = os.path.join('modif_transcriptomie.txt')
dataset = sc.read(data_path)
# %%


atlas = pd.read_csv(path_file + file_atlas_reference)
value_name = list(atlas['name'])

# %%
num_locations = len(atlas)
locations_apriori = atlas[:num_locations][['X', 'Y', 'Z']].values
# %%
tgm.traduction_gene_marqueur(file_atlas_reference, time_dataset, fichier_dataset, fichier_position)
# %%

atlas_path = os.path.join("gene_marker_news.txt")
atlas = sc.read(atlas_path)
atlas_genes = atlas.var.index.tolist()
atlas.obsm['spatial'] = locations_apriori


novosparc.pl.embedding(atlas, ['hmg-11'])
tissue = novosparc.cm.Tissue(dataset=dataset, locations=locations_apriori)

markers = atlas_genes[:99]
atlas_matrix = atlas.to_df()[markers].values
markers_idx = pd.DataFrame({'markers_idx': np.arange(num_genes)}, index=gene_names)
# %%
markers_to_use = np.concatenate(markers_idx.loc[markers].values)

# %%
atlas_matrix = atlas_matrix.astype(object)
# %%


# tissue.markers_to_use = tissue.dge[:, list(markers_to_use)] / np.amax(tissue.dge[:, list(markers_to_use)])

cell_expression = tissue.dge[:, markers_to_use] / np.amax(tissue.dge[:, markers_to_use])
atlas_expression = atlas_matrix / np.amax(atlas_matrix)

# %%
cell_expression = cell_expression.astype(float)
atlas_expression = atlas_expression.astype(float)

# %%
tissue.costs['markers'] = cdist(cell_expression, atlas_expression)
tissue.markers_to_use = markers_to_use
tissue.num_markers = len(tissue.markers_to_use)
# %%
tissue.setup_smooth_costs(num_neighbors_s=num_neighbors_s, num_neighbors_t=num_neighbors_t, verbose=True)

# %%
# tissue.setup_reconstruction(atlas_matrix=atlas_matrix,markers_to_use=list(markers_to_use), num_neighbors_s=num_neighbors_s,num_neighbors_t=num_neighbors_t)
# %%
# tissue.setup_reconstruction(num_neighbors_s=num_neighbors_s,num_neighbors_t=num_neighbors_t)

tissue.reconstruct(alpha_linear=alpha_linear, epsilon=epsilon)

sdge = tissue.sdge
dataset_reconst = sc.AnnData(pd.DataFrame(sdge.T, columns=gene_names))
dataset_reconst.obsm['spatial'] = locations_apriori

novosparc.pl.embedding(dataset_reconst,['hmg-11'])



