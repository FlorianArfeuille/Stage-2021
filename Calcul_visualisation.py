__author__ = 'florian'
__Filename__ = 'code_propre'
__Creationdate__ = '02/07/2021'

# %% Parameters data:
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
import traduction_gene_marqueur
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
sdge = sdge.T
sdge = pd.DataFrame(sdge, columns=gene_names)
# %%


valueX = []
valueY = []
valueZ = []

for i in range(len(locations_apriori)):
    valueX.append(locations_apriori[i][0])
    valueY.append(locations_apriori[i][1])
    valueZ.append(locations_apriori[i][2])
sdge.insert(0, 'X', valueX, allow_duplicates=False)
sdge.insert(1, 'Y', valueY, allow_duplicates=False)
sdge.insert(2, 'Z', valueZ, allow_duplicates=False)
sdge.insert(3, 'name', value_name, allow_duplicates=False)
sdge.to_csv('matrice gene position gène')
# %%


gw = tissue.gw

# %%

np.savetxt("matrice cellule position.txt", gw, fmt='%.18f')

# %%
embedding(dataset_reconst, ['cul-2', 'spc-1'])




# %%
# nom_cell1 = ['Z2', 'Z3', 'MSppapp', 'Caapa', 'MSaaapp', 'MSapaap', 'MSpaapp', 'MSppaap', 'MSpappa', 'MSpappp', 'MSaappa', 'MSaappp', 'MSappaa', 'MSappap', 'MSpppaa', 'MSpppap', 'MSppppa', 'MSppppp', 'MSapppa', 'MSapppp', 'Dapa', 'Dapp', 'Dppa', 'Dppp', 'Dpaa', 'Dpap', 'ABplpappaa', 'ABplpappap', 'Daaa', 'Daap', 'ABarapappa', 'ABarapappp', 'ABalpapppa', 'ABalpapppp', 'ABalapappa', 'ABalapappp', 'ABalappppa', 'ABalappppp', 'ABalpappaa', 'ABalpappap', 'ABalppappa', 'ABalppappp', 'ABalpppppa', 'ABalpppppp', 'ABarapapaa', 'ABarapapap', 'ABarpaappa', 'ABarpaappp', 'ABplaapppa', 'ABplaapppp', 'ABplapapaa', 'ABplapapap', 'ABplapappa', 'ABplapappp', 'ABplappaaa', 'ABplappaap', 'ABplappapa', 'ABplappapp', 'ABplpaappa', 'ABplpaappp', 'ABpraaappa', 'ABpraaappp', 'ABprpaappa', 'ABprpaappp', 'ABprpappaa', 'ABprpappap', 'Caaaaa', 'Caaaap', 'ABprappaaa', 'ABprappaap', 'ABprappapa', 'ABprappapp', 'Cpaapa', 'Cpaapp', 'ABalapaaaa', 'ABalapaaap', 'ABalppapaa', 'ABalppapap', 'ABarappapa', 'ABarappapp', 'ABarappppa', 'ABarappppp', 'ABarpaapaa', 'ABarpaapap', 'ABarpapppa', 'ABarpapppp', 'ABplapaaaa', 'ABplapaaap', 'ABplapaapa', 'ABplapaapp', 'ABplpaaaaa', 'ABplpaaaap', 'ABplpaaapa', 'ABplpaaapp', 'ABplpapaaa', 'ABplpapaap', 'ABplpapapa', 'ABplpapapp', 'ABplppppaa', 'ABplppppap', 'ABpraapppa', 'ABpraapppp', 'ABprapapaa', 'ABprapapap', 'ABprapappa', 'ABprapappp', 'ABprpapapa', 'ABprpapapp', 'Caaapa', 'Caaapp', 'Cpaaaa', 'Cpaaap', 'ABalapaapa', 'ABalapaapp', 'ABalapapaa', 'ABalapapap', 'ABalapppaa', 'ABalapppap', 'ABalpaapaa', 'ABalpaapap', 'ABalppppaa', 'ABalppppap', 'ABaraaapaa', 'ABaraaapap', 'ABarapppaa', 'ABarapppap', 'ABplpapppa', 'ABplpapppp', 'ABplpppppa', 'ABplpppppp', 'ABprapaaaa', 'ABprapaaap', 'ABprpaaaaa', 'ABprpaaaap', 'ABprpaaapa', 'ABprpaaapp', 'ABprpaapaa', 'ABprpaapap', 'ABprpapaaa', 'ABprpapaap', 'ABprpapppa', 'ABprpapppp', 'ABprppppaa', 'ABprppppap', 'ABprpppppa', 'ABprpppppp', 'ABalaapppa', 'ABalaapppp', 'ABalappapa', 'ABalappapp', 'ABalpapaaa', 'ABalpapaap', 'ABalpapapa', 'ABalpapapp', 'ABalpppapa', 'ABalpppapp', 'ABarpppppa', 'ABarpppppp', 'ABplaaappa', 'ABplaaappp', 'ABplpaapaa', 'ABplpaapap', 'ABplpppaaa', 'ABplpppaap', 'ABplpppapa', 'ABplpppapp', 'ABpraaaapa', 'ABpraaaapp', 'ABpraaapaa', 'ABpraaapap', 'ABprapaapa', 'ABprapaapp', 'ABprpppaaa', 'ABprpppaap', 'ABprpppapa', 'ABprpppapp', 'ABalaaapal', 'ABalaaapar', 'ABalaaappl', 'ABalaaappr', 'ABalpaappa', 'ABalpaappp', 'ABalppaapa', 'ABalppaapp', 'ABaraapaaa', 'ABaraapaap', 'ABaraapapa', 'ABaraapapp', 'ABaraappaa', 'ABaraappap', 'ABaraapppa', 'ABaraapppp', 'ABarapaapa', 'ABarapaapp', 'ABarappaaa', 'ABarappaap', 'ABarpapapa', 'ABarpapapp', 'ABarppaapa', 'ABarppaapp', 'ABarppapaa', 'ABarppapap', 'ABarppappa', 'ABarppappp', 'ABarpppapa', 'ABarpppapp', 'ABarppppaa', 'ABarppppap', 'ABplaaaaaa', 'ABplaaaaap', 'ABplaaaapa', 'ABplaaaapp', 'ABplaapaaa', 'ABplaapaap', 'ABplaapapa', 'ABplaapapp', 'ABplaappaa', 'ABplaappap', 'ABplapppaa', 'ABplapppap', 'ABplappppa', 'ABplappppp', 'ABplppaaaa', 'ABplppaaap', 'ABplppapaa', 'ABplppapap', 'ABplppappa', 'ABplppappp', 'ABpraapaaa', 'ABpraapaap', 'ABpraapapa', 'ABpraapapp', 'ABpraappaa', 'ABpraappap', 'ABprapppaa', 'ABprapppap', 'ABprappppa', 'ABprappppp', 'ABprppaaaa', 'ABprppaaap', 'ABprppapaa', 'ABprppapap', 'ABprppappa', 'ABprppappp', 'Cpapaa', 'Cpapap', 'MSaaaaaa', 'MSaaaaap', 'MSpaaaaa', 'MSpaaaap', 'ABarpapaaa', 'ABarpapaap', 'ABaraaappa', 'ABaraaappp', 'ABalpaaaaa', 'ABalpaaaap', 'ABalpaaapa', 'ABalpaaapp', 'ABalppaaaa', 'ABalppaaap', 'ABalpppaaa', 'ABalpppaap', 'ABarapaaaa', 'ABarapaaap', 'ABarppaaaa', 'ABarppaaap', 'ABarpppaaa', 'ABarpppaap', 'ABplppaapa', 'ABplppaapp', 'ABpraaaaaa', 'ABpraaaaap', 'ABprppaapa', 'ABprppaapp', 'Caappd', 'Caappv', 'Cpappd', 'Cpappv', 'Eplaa', 'Eplap', 'Epraa', 'Eprap', 'MSaaaapa', 'MSaaaapp', 'MSaapaaa', 'MSaapaap', 'MSpaaapa', 'MSpaaapp', 'Ealaa', 'Ealap', 'Ealpa', 'Ealpp', 'Earaa', 'Earap', 'Earpa', 'Earpp', 'ABalaaaala', 'ABalaaaalp', 'ABalaaaarl', 'ABalaaaarr', 'ABalaapaaa', 'ABalaapaap', 'ABalaapapa', 'ABalaapapp', 'ABalaappaa', 'ABalaappap', 'ABalappaaa', 'ABalappaap', 'ABaraaaaaa', 'ABaraaaaap', 'ABarpaaaaa', 'ABarpaaaap', 'ABarpaaapa', 'ABarpaaapp', 'ABarpappaa', 'ABarpappap', 'ABplaaapaa', 'ABplaaapap', 'Cappaa', 'Cappap', 'Capppa', 'Capppp', 'Cpppaa', 'Cpppap', 'Cppppa', 'Cppppp', 'MSaaapaa', 'MSaaapap', 'MSaapapa', 'MSaapapp', 'MSapaaaa', 'MSapaaap', 'MSapapaa', 'MSapapap', 'MSpapapa', 'MSpapapp', 'MSppaaaa', 'MSppaaap', 'MSppapaa', 'MSppapap', 'MSpapaaa', 'MSpapaap', 'ABaraaaapa', 'ABaraaaapp', 'Capaaa', 'Capaap', 'Capapa', 'Capapp', 'Cppaaa', 'Cppaap', 'Eplpa', 'Eplpp', 'Eprpa', 'Eprpp', 'MSpaapaa', 'MSpaapap', 'Cppapa', 'Cppapp', 'MSapappa', 'MSapappp']
# data set
# nom_cell = ['Caapa', 'ABpraapppp', 'ABalapppap', 'ABplaapaap', 'Ealpp', 'ABprppappp', 'Cpppaa', 'ABprappaap', 'MSaaapp', 'MSpaapp', 'MSpapapa', 'ABalaapppa', 'ABarppaapp', 'ABprapaapa', 'ABprpapaap', 'ABplaapapp', 'ABalpapppa', 'ABalaapppp', 'ABarpaaapa', 'Epra', 'ABarpapaaa', 'ABalaaappr', 'MSpapaaa', 'ABprpaaaap', 'ABplaaaaaa', 'ABprapppap', 'ABalaaaarl', 'ABplppappp', 'ABarappppa', 'ABalaaappl', 'ABarpaaaaa', 'ABplpppaap', 'ABprapppaa', 'ABplappaaa', 'MSaaaaaa', 'Cappap', 'Cppaa', 'ABarpppaaa', 'ABalappaa', 'Earpa', 'Cpppa', 'ABarapppap', 'ABplppppap', 'ABplpapppp', 'Ealaa', 'ABprpaaaaa', 'ABplapaaap', 'ABalppapaa', 'ABplpppaaa', 'MSapapa', 'ABarppppaa', 'Capppp', 'Cpppap', 'ABplppapaa', 'ABalpappaa', 'MSaappp', 'ABarpaaapp', 'ABalaaaalp', 'ABarapapap', 'ABprpaappp', 'Caapp', 'Earaa', 'Epla', 'MSpaapaa', 'MSpapaap', 'ABalaaapar', 'ABplaaapaa', 'MSaapaaa', 'ABplaapapa', 'ABarppaaap', 'ABplpaapaa', 'ABarpapaap', 'ABplpaappa', 'MSppapaa', 'ABarapaapa', 'ABpraaappa', 'MSaapaap', 'MSaaapap', 'ABplpapppa', 'ABarappaaa', 'ABplapaapp', 'ABplpaappp', 'ABprappapp', 'ABplapapap', 'Cpappd', 'ABalpapaap', 'ABarapaaap', 'MSapapap', 'ABaraapapa', 'ABalaapaaa', 'ABarpppppp', 'ABalppaaaa', 'ABprppaaaa', 'ABplpaaapa', 'MSaapapa', 'Cpaapp', 'ABaraaaap', 'MSppaaaa', 'ABprpppapa', 'Eplpp', 'Caappd', 'ABpraapppa', 'ABpraapapp', 'ABaraapaaa', 'MSppaap', 'ABprppppaa', 'ABplpapaap', 'MSappaa', 'ABarappppp', 'ABarppaapa', 'Cpppp', 'ABaraaapap', 'ABarapaaaa', 'ABalppapap', 'ABpraaaaap', 'Z2', 'ABplpapaaa', 'ABprppapap', 'Earap', 'ABalapappa', 'ABalaappaa', 'ABplppaaap', 'MSappap', 'ABplaaappp', 'ABplapaaaa', 'MSaaaapa', 'ABplpaaaap', 'ABpraaaapp', 'ABarpapppa', 'ABarapapaa', 'ABalpppapp', 'ABalappaaa', 'ABprppaaap', 'Ealpa', 'ABplappppa', 'ABpraaappp', 'Dpap', 'Eara', 'Capaa', 'ABpraapapa', 'ABalaapapa', 'ABarappapp', 'ABprppappa', 'Cpappv', 'Cppaap', 'ABalpaaapp', 'ABplaaaaap', 'ABpraappap', 'ABaraaapaa', 'MSpappp', 'ABalaappap', 'Cappa', 'ABpraaapap', 'MSpaaaaa', 'ABplppapap', 'Cpaapa', 'Cppapp', 'ABprpapapa', 'Cppppa', 'ABprapappa', 'ABarppaaaa', 'ABalpaappp', 'MSpaapap', 'ABplaappaa', 'ABalaaapal', 'Capapa', 'ABarpaaaap', 'MSpppaa', 'ABprapappp', 'MSppapa', 'ABaraappap', 'ABalapaapp', 'ABalpppppp', 'ABplppaapp', 'ABprpppapp', 'ABarpaapaa', 'Earpp', 'MSapaaap', 'Caaapa', 'ABalaapaap', 'ABplaappap', 'ABarpappap', 'Caaaap', 'ABalapapaa', 'ABprpaappa', 'ABarpaaap', 'ABalppppap', 'ABalpapppp', 'ABalapaaaa', 'Daap', 'ABprpapppa', 'ABalpppppa', 'ABplaapppp', 'ABalppaapp', 'Cappp', 'ABalaaaarr', 'ABplappapp', 'ABalppappa', 'ABalapppaa', 'ABplapppap', 'ABalpaaapa', 'ABpraaaaaa', 'ABarpppppa', 'ABprppppap', 'Cppap', 'ABarpapapa', 'ABplpppapp', 'Eprap', 'ABalapaaap', 'ABarpappaa', 'ABplppaapa', 'ABplpappaa', 'ABplapappa', 'Cppppp', 'ABplppaaaa', 'MSpaaaap', 'MSppppp', 'ABalpaaaaa', 'ABalpppaap', 'ABalpappap', 'ABalpaapap', 'Daaa', 'Eprp', 'ABalpaaaap', 'ABarppappa', 'ABplpppapa', 'ABpraappaa', 'Dapp', 'Eprpa', 'ABpraaaapa', 'MSppppa', 'ABprappaaa', 'Caaaaa', 'ABalappppa', 'ABplaaaapa', 'ABprpappaa', 'ABplaaaapp', 'ABplaapaaa', 'ABprpppaaa', 'ABprpaapap', 'ABplappppp', 'ABalappapa', 'Capap', 'Cppapa', 'ABplappapa', 'ABplaapppa', 'Epraa', 'ABalppaapa', 'ABalaaaala', 'ABarpaappp', 'MSapaap', 'ABalpapaaa', 'ABalapappp', 'ABarappaap', 'ABprapaaaa', 'ABaraaaaap', 'ABalappapp', 'ABprpapapp', 'ABaraapapp', 'ABarpapapp', 'ABplppppaa', 'MSpaaapa', 'Cpapaa', 'MSppapap', 'Cpaaap', 'ABplpaaapp', 'ABprpapaaa', 'ABplappaap', 'ABprapaapp', 'Eprpp', 'Cppaaa', 'ABarppapa', 'MSapaaaa', 'ABplpppppp', 'ABprapaaap', 'ABpraapaap', 'Capppa', 'MSaaapaa', 'ABaraaaapa', 'ABalpppapa', 'MSapppa', 'ABplpapapp', 'Capaap', 'ABplpaaaaa', 'MSpapapp', 'ABalpapapp', 'ABalappaap', 'ABaraaappp', 'ABarpaapap', 'ABplppappa', 'ABprpppppp', 'ABprapapap', 'ABplapppaa', 'ABprappapa', 'ABaraaappa', 'ABalpaappa', 'ABprpapppp', 'ABprappppp', 'Dppp', 'ABplapapaa', 'ABarppappp', 'ABaraappaa', 'ABalppaaap', 'MSapapaa', 'ABprpaapaa', 'MSapaaa', 'Dapa', 'ABpraaapaa', 'ABprpppaap', 'ABalpppaaa', 'Eplaa', 'MSppaaap', 'ABplpaapap', 'ABplpappap', 'MSapapp', 'ABarppppap', 'Z3', 'ABarpppapa', 'ABarpapppp', 'ABprpppppa', 'ABpraapaaa', 'ABplaaapap', 'Dpaa', 'MSapppp', 'ABplapappp', 'Cpaaaa', 'ABarppapp', 'ABplpppppa', 'MSpaaapp', 'ABprppaapa', 'ABprpaaapp', 'ABprapapaa', 'ABarpppaap', 'Eplp', 'ABaraapppa', 'ABalpapapa', 'Eplap', 'Cappaa', 'Earp', 'MSpppap', 'ABarapappa', 'ABalppappp', 'ABalppppaa', 'Capaaa', 'ABarapppaa', 'ABplaaappa', 'ABprpaaapa', 'ABarppapaa', 'Caappv', 'ABaraaaaaa', 'MSaappa', 'ABalappppp', 'MSppaaa', 'MSppapp', 'ABplpapapa', 'ABarappapa', 'Caaapp', 'ABalapapap', 'ABarapaapp', 'MSaaaapp']
# genes marqueurs

d_visu = pd.read_csv(path_file+'matrice gene position gène')
M = np.loadtxt(path_file+'matrice cellule position.txt')

# %%
nom_cell = list(d_visu['name'])
nom_cell1 = list(np.load(path_file + 'cell_name_dataset.npy'))

# %%
M = M * np.shape(M)[0]

import plotly.express as px
import plotly.graph_objects as go

liste_gene_e = list(d_visu.columns)

# %%

fig = px.scatter_3d(d_visu, x='X', y='Y', z='Z', color=liste_gene_e[20],
                    # animation_frame="time", animation_group="cell",
                    hover_name=nom_cell, range_x=[100, 600], range_y=[150, 350], range_z=[40, 100])

nb_cell = np.shape(M)[0]

from sklearn.decomposition import PCA

l1, l2, l3 = 400, 400, 80


def pca_value(mat):
    X = mat
    pca = PCA(n_components=3)
    pca.fit(X)
    vp = pca.explained_variance_
    vs = pca.singular_values_
    axis = pca.components_.T
    # axis /= axis.std()
    return vp, axis


def pca(mat, d):
    v, w = pca_value(mat)
    # w = np.transpose(wt)
    p = np.asarray(d[['X', 'Y', 'Z']])
    px = np.mean(p[:, 0])
    py = np.mean(p[:, 1])
    pz = np.mean(p[:, 2])
    # pmean = np.array([px, py, pz])
    # d1 = v[0]*w[0,:]
    # d2 = v[1]*w[1,:]
    # d3 = v[2]*w[2,:]
    # d1 = (v[0]/2)*w[:,0]
    # d2 = (v[1]/2)*w[:,1]
    # d3 = (v[2]/2)*w[:,2]

    d1 = w[:, 0]
    d2 = w[:, 1]
    d3 = w[:, 2]

    d1x = [l1 * d1[0] + px, px, -d1[0] * l1 + px]
    d1y = [l1 * d1[1] + py, py, -d1[1] * l1 + py]
    d1z = [l1 * d1[2] + pz, pz, -d1[2] * l1 + pz]

    d2x = [d2[0] * l2 + px, px, px - d2[0] * l2]
    d2y = [d2[1] * l2 + py, py, py - d2[1] * l2]
    d2z = [d2[2] * l2 + pz, pz, pz - d2[2] * l2]

    d3x = [d3[0] * l3 + px, px, px - d3[0] * l3]
    d3y = [d3[1] * l3 + py, py, py - d3[1] * l3]
    d3z = [d3[2] * l3 + pz, pz, pz - d3[2] * l3]

    return [d1x, d1y, d1z], [d2x, d2y, d2z], [d3x, d3y, d3z]


def trace_pca_e(mat, d):
    # fig = go.Figure()
    axe = pca(mat, d)
    axe1 = axe[0]
    axe2 = axe[1]
    axe3 = axe[2]

    trace1 = go.Scatter3d(
        x=axe1[0],
        y=axe1[1],
        z=axe1[2],
        mode='lines',
        name='AP')

    trace2 = go.Scatter3d(
        x=axe2[0],
        y=axe2[1],
        z=axe2[2],
        mode='lines',
        name='LR')

    trace3 = go.Scatter3d(
        x=axe3[0],
        y=axe3[1],
        z=axe3[2],
        mode='lines',
        name='DV')
    return trace1, trace2, trace3


tr1, tr2, tr3 = trace_pca_e(np.asarray(d_visu[['X', 'Y', 'Z']]), d_visu)
fig = go.Figure()

fig.add_trace(tr1)
fig.add_trace(tr2)
fig.add_trace(tr3)

fig_ref = go.Figure()
fig_ref.add_trace(tr1)
fig_ref.add_trace(tr2)
fig_ref.add_trace(tr3)


def update_fig_ref(celule):
    fig = go.Figure()
    fig.add_trace(tr1)

    fig.add_trace(
        go.Scatter3d(
            x=[-l1, 0, l1],
            y=[0, 0, 0],
            z=[0, 0, 0],
            mode="text",
            text=["Posterior", "", "Anterior"],
            textposition="bottom center",
            name='Anterior / Posterior'
        ))

    fig.add_trace(tr2)

    fig.add_trace(
        go.Scatter3d(
            x=[0, 0, 0],
            y=[-l2, 0, l2],
            z=[0, 0, 0],
            mode="text",
            text=["Left", "", "Right"],
            textposition="bottom center",
            name='Left / Right'
        ))

    fig.add_trace(tr3)
    fig.add_trace(
        go.Scatter3d(
            x=[0, 0, 0],
            y=[0, 0, 0],
            z=[-l3, 0, l3],
            mode="text",
            text=["Dorsal", "", "Ventral"],
            textposition="bottom center",
            name='Dorsal / Ventral'
        ))

    # fig_fct = fig
    indice_cellule = nom_cell.index(nom_cell1[celule])
    couleur = []
    for c in range(nb_cell):
        if c != indice_cellule:
            couleur.append(0)
        if c == indice_cellule:
            couleur.append(1)
    fig.add_trace(go.Scatter3d(x=d_visu['X'], y=d_visu['Y'], z=d_visu['Z'],
                               mode='markers',
                               marker={'color': couleur,
                                       'opacity': 0.5,

                                       'colorscale': "Viridis"
                                       }

                               ,  # animation_frame="time", animation_group="cell",
                               # size = 'size',
                               name='Cells',
                               hovertext=nom_cell,

                               # hoverinfo = 'name'
                               # ,  range_x=[-330,290], range_y=[-200,210], range_z=[-100,90]
                               ))
    return fig


# %%


## A incorporer avec la fonction de print


# %%


# %%
from cross_validation import CrossValidation as cv

# %%
cross_val = cv(atlas=atlas,
               location_a_priori=locations_apriori,
               alpha=alpha_linear,
               dataset=dataset,
               resultat_novosparc=M)

# %%
gene_to_delete = 'hlh-16'

sdge_cross_val = cross_val.calcul_novosparc(
    dataset=cross_val.dataset,
    locations_apriori=cross_val.locations_apriori,
    num_neighbors=20, alpha_i=0.5, epsilon_i=1e-2,
    atlas_genes=cross_val.atlas_gene,
    atlas=cross_val.atlas,
    num_genes=cross_val.num_genes,
    gene_names=cross_val.gene_names,
    num_locations=cross_val.num_locations,
)

# %%


qi, sdge_comp = cross_val.calcul_novosparc_tilt(
    dataset=cross_val.dataset,
    nom_gene_to_delete=cross_val.atlas_gene[0],
    locations_apriori=cross_val.locations_apriori,
    num_neighbors=20, alpha_i=0.5, epsilon_i=1e-2,
    atlas_genes=cross_val.atlas_gene,
    atlas=cross_val.atlas,
    num_genes=cross_val.num_genes,
    gene_names=cross_val.gene_names,
    num_locations=cross_val.num_locations,
    Xatlas=cross_val.Xatlas)

# %% validation score


cross_val.calcul_cross_validation_score(dataset=cross_val.dataset,
                                        liste_gene_loop=cross_val.atlas_gene,
                                        locations_apriori=cross_val.locations_apriori,
                                        num_neighbors=20, alpha_i=0.8, epsilon_i=1e-2,
                                        atlas_genes=cross_val.atlas_gene,
                                        atlas=cross_val.atlas,
                                        num_genes=cross_val.num_genes,
                                        gene_names=cross_val.gene_names,
                                        num_locations=cross_val.num_locations,
                                        Xatlas=cross_val.Xatlas)

# %% example of validation score according to alpha

Q_according_alpha = []
qi_according_alpha = []
alpha_range = [0, 0.1,0.2,0.3,0.4, 0.5,0.6,0.7,0.8,0.9, 1]
for alpha in alpha_range :
    cross_val_alpha  = cv(atlas=atlas,
               location_a_priori=locations_apriori,
               alpha=alpha,
               dataset=dataset,
               resultat_novosparc=M)
    liste_qi_alpha, Q_alpha = cross_val.calcul_cross_validation_score(dataset=cross_val.dataset,
                                        liste_gene_loop=cross_val.atlas_gene,
                                        locations_apriori=cross_val.locations_apriori,
                                        num_neighbors=20, alpha_i=alpha, epsilon_i=1e-2,
                                        atlas_genes=cross_val.atlas_gene,
                                        atlas=cross_val.atlas,
                                        num_genes=cross_val.num_genes,
                                        gene_names=cross_val.gene_names,
                                        num_locations=cross_val.num_locations,
                                        Xatlas=cross_val.Xatlas)
    qi_according_alpha.append(liste_qi_alpha)
    Q_according_alpha.append(Q_alpha)
#%%
fig_Q_according_alpha = px.line(x = alpha_range, y = Q_according_alpha)
fig_Q_according_alpha.update_layout(xaxis={
       'title' : 'alpha '},
    yaxis={'title' : 'Q'},
    title = 'Evolution of Q according alpha')


# %%

print(cross_val.Q)
fighisto = cross_val.histogramme_qi(cross_val.liste_qi)

# %%

import dash
import dash_html_components as html
import dash_core_components as dcc
from plotly.subplots import make_subplots
from dash.dependencies import Input, Output

config = {'modeBarButtonsToAdd': ['drawline',
                                  'drawopenpath',
                                  'drawclosedpath',
                                  'drawcircle',
                                  'drawrect',
                                  'eraseshape'
                                  ]}
app1 = dash.Dash(__name__)
app1.layout = html.Div([
    html.Div([
        html.Div(['''
3D prediction: you can choose between cell which make a color plot of probability embedding or
gene which embed the selected gene distribution
''']),
        dcc.Dropdown(
            id='slct_color',

            options=[
                {'label': 'gene', 'value': 'gene'},
                {'label': 'cell', 'value': 'cell'},

            ],
            value='gene',
            style={'width': "50%"}
        ),

        dcc.Dropdown(
            id='slct_gene_prior',
            options=[{'label': liste_gene_e[i], 'value':
                liste_gene_e[i]} for i in range(len(liste_gene_e))],
            value=liste_gene_e[10]
        ),

        dcc.Dropdown(
            id='slct_celule',
            options=[{'label': nom_cell1[i], 'value': i} for i in range(nb_cell)],
            value=10
        ),

        dcc.Graph(id='prior', style={'width': '160vh',
                                     'height': '90vh'},
                  figure=fig,
                  config=config,
                  ),
        html.Div([''' 3D reference : highlight cell selected true position
''']),
        dcc.Graph(id='ref', style={'width': '160vh',
                                   'height': '90vh'},
                  figure=fig,
                  config=config,
                  ),
        dcc.Graph(id='histo', style={'width': '160vh',
                                     'height': '90vh'},
                  figure=fighisto,
                  config=config,
                  ),
        html.Div(['''
Select gene to delete (in atlas reference) to compute new novo sparc prediction :
''']),
        dcc.Dropdown(
            id='slct_gene_delete',
            options=[{'label': cross_val.atlas_gene[i], 'value':
                cross_val.atlas_gene[i]} for i in range(len(cross_val.atlas_gene))],
            value=cross_val.atlas_gene[0]
        ),

        html.Div(['''
Select gene to visualize for comparaison
''']),
        dcc.Dropdown(
            id='slct_gene_prior_comparaison',
            options=[{'label': liste_gene_e[5:][i], 'value':
                liste_gene_e[5:][i]} for i in range(len(liste_gene_e[5:]))],
            value='cul-2'
        ),

        dcc.Graph(id='comparaison', style={'width': '160vh',
                                     'height': '90vh'},
                  figure=fighisto,
                  config=config,
                  ),

        dcc.Graph(id='Q_alpha', style={'width': '160vh',
                                     'height': '90vh'},
                  figure=fig_Q_according_alpha,
                  config=config,
                  ),
        html.Div(['''
Select alpha
''']),
        dcc.Dropdown(
            id='slct_alpha',
            options=[{'label': alpha_range[i], 'value':
                i} for i in range(len(alpha_range))],
            value=3
        ),

        dcc.Graph(id='qi_alpha', style={'width': '160vh',
                                     'height': '90vh'},
                  figure=fighisto,
                  config=config,
                  )
    ])
])


@app1.callback(
    [
        dash.dependencies.Output('prior', 'figure'),
        dash.dependencies.Output('ref', 'figure')

    ],
    [dash.dependencies.Input('slct_gene_prior', 'value'),
     dash.dependencies.Input('slct_color', 'value'),
     dash.dependencies.Input('slct_celule', 'value')])
def update_prior(gene, choix, celule):
    # print(file)
    # d_prior_360 = fonction_create_df_prior(180, 10),
    if choix == 'gene':
        couleur = list(d_visu[gene])
    if choix == 'cell':
        couleur = list(M[celule, :])

    # liste_gene_e = list(d_visu.columns)
    # liste_gene_e = liste_gene_e[4:]

    # Création d'une nouvelle figure
    fig = go.Figure()

    # ajout des axes principaux
    fig.add_trace(tr1)

    fig.add_trace(
        go.Scatter3d(
            x=[-l1, 0, l1],
            y=[0, 0, 0],
            z=[0, 0, 0],
            mode="text",
            text=["Posterior", "", "Anterior"],
            textposition="bottom center",
            name='Anterior / Posterior'
        ))

    fig.add_trace(tr2)

    fig.add_trace(
        go.Scatter3d(
            x=[0, 0, 0],
            y=[-l2, 0, l2],
            z=[0, 0, 0],
            mode="text",
            text=["Left", "", "Right"],
            textposition="bottom center",
            name='Left / Right'
        ))

    fig.add_trace(tr3)
    fig.add_trace(
        go.Scatter3d(
            x=[0, 0, 0],
            y=[0, 0, 0],
            z=[-l3, 0, l3],
            mode="text",
            text=["Dorsal", "", "Ventral"],
            textposition="bottom center",
            name='Dorsal / Ventral'
        ))
    # fin des axes principaux
    # fig_fct = fig
    fig.add_trace(go.Scatter3d(x=d_visu['X'], y=d_visu['Y'], z=d_visu['Z'],
                               mode='markers',
                               marker={'color': couleur,

                                       'colorscale': "Viridis",
                                       'colorbar': {'thickness': 10}
                                       }

                               ,  # animation_frame="time", animation_group="cell",
                               # size = 'size',
                               name='Cells',
                               hovertext=nom_cell,

                               # hoverinfo = 'name'
                               # ,  range_x=[-330,290], range_y=[-200,210], range_z=[-100,90]
                               ))
    fig.update_layout(coloraxis=dict(colorscale='Viridis'), showlegend=False)
    fig.update_coloraxes(showscale=True)
    fig.update_layout(coloraxis_showscale=True)

    fig_ref = update_fig_ref(celule)

    return [fig, fig_ref]


@app1.callback(
    [
        dash.dependencies.Output('comparaison', 'figure')

    ],
    [dash.dependencies.Input('slct_gene_delete', 'value'),
        dash.dependencies.Input('slct_gene_prior_comparaison', 'value')]
     )
def update_prior_comparaison(gene_to_delete, gene):

    q, sdge_comp = cross_val.calcul_novosparc_tilt(
    dataset=cross_val.dataset,
    nom_gene_to_delete=gene_to_delete,
    locations_apriori=cross_val.locations_apriori,
    num_neighbors=20, alpha_i=0.5, epsilon_i=1e-2,
    atlas_genes=cross_val.atlas_gene,
    atlas=cross_val.atlas,
    num_genes=cross_val.num_genes,
    gene_names=cross_val.gene_names,
    num_locations=cross_val.num_locations,
    Xatlas=cross_val.Xatlas)

    couleur = list(d_visu[gene]/np.max(d_visu[gene]))
    gene_index_cross = cross_val.gene_names.index(gene)
    couleur_comparaison = list(sdge_comp[:,gene_index_cross]/np.max(sdge_comp[:,gene_index_cross]))
    couleur_comparaison_value = list(sdge_comp[:,gene_index_cross])

    # Création d'une nouvelle figure
    fig = make_subplots(rows=1, cols=2,
    subplot_titles = ("Novo sparc prediction", "Prediction computed without  "+str(gene_to_delete) +" in reference atlas." ),
                        specs=[[{"type": "scene"}, {"type": "scene"}]]

                        )

    fig.add_trace(go.Scatter3d(x=d_visu['X'], y=d_visu['Y'], z=d_visu['Z'],
                               mode='markers',
                               marker={'color': couleur,

                                       'colorscale': "Viridis",
                                       'colorbar': {'thickness': 10},
                                       'coloraxis': "coloraxis"
                                       }

                               ,  # animation_frame="time", animation_group="cell",
                               # size = 'size',
                               name=str(gene)+' expression ',
                               hovertext=nom_cell,

                               # hoverinfo = 'name'
                               # ,  range_x=[-330,290], range_y=[-200,210], range_z=[-100,90]
                               ), row=1, col = 1 )

    fig.add_trace(go.Scatter3d(x=d_visu['X'], y=d_visu['Y'], z=d_visu['Z'],
                               mode='markers',
                               marker={'color': couleur_comparaison,

                                       'colorscale': "Viridis",
                                       'colorbar': {'thickness': 10},
                                       'coloraxis' : "coloraxis"
                                       }

                               ,  # animation_frame="time", animation_group="cell",
                               # size = 'size',
                               name=str(gene)+' expression ',
                               hovertext= nom_cell,
                               #hovertemplate= np.array([nom_cell, couleur_comparaison]),

                               # hoverinfo = 'name'
                               # ,  range_x=[-330,290], range_y=[-200,210], range_z=[-100,90]
                               ), row=1, col=2)
    fig.update_layout(coloraxis=dict(colorscale='Viridis'), showlegend=False)

    fig.update_traces(hovertemplate="<b>%{text}</b><br><br>" +
                                    '<i>Cell name</i> : %{text} <br>' +
                                    'Normalized gene expression: %{marker.color}', text=nom_cell, row=1,
                      col=1)
    fig.update_traces(hovertemplate="<b>%{text}</b><br><br>" +
                                    '<i>Cell name</i> : %{text} <br>'+
                                    'Normalized gene expression: %{marker.color}<br>'
                                    #'True expression of gene expression : %{hoverinfo}'
                                    , text=nom_cell,  row=1,
                      col=2)

    fig.update_coloraxes(showscale=True)
    fig.update_layout(coloraxis_showscale=True)

    #fig_ref = update_fig_ref(celule)

    return [fig]


@app1.callback(
    [
        dash.dependencies.Output('qi_alpha', 'figure')

    ],
    [
     dash.dependencies.Input('slct_alpha', 'value')
        ]
     )
def update_prior_comparaison( select_alpha):
    #index_gene_to_delete = cross_val.gene_names.index(gene_to_delete)
    liste_qi_selected = qi_according_alpha[select_alpha]
    fig = px.histogram(liste_qi_selected)
    fig.update_layout(xaxis={
        'title' :'Range value for qi '},
    yaxis={'title':'count'},
    title = 'Histogram of qi according to alpha ='+str(alpha_range[select_alpha]))
    return [fig]

# %%


if __name__ == '__main__':
    app1.run_server()
