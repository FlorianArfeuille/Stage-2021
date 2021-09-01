
import pandas as pd
import numpy as np
import scanpy as sc
import novosparc
from scipy.spatial import distance
from scipy.spatial.distance import cdist
from sklearn.preprocessing import normalize
import plotly.express as px



def update_atlas_variable(atlas_matrix, tissue, markers_to_use ) :
    atlas_matrix = atlas_matrix.astype(object)

    cell_expression = tissue.dge[:, markers_to_use] / np.amax(tissue.dge[:, markers_to_use])

    cell_expression = cell_expression.astype(float)

    atlas_expression = atlas_matrix / np.amax(atlas_matrix)
    atlas_expression = atlas_expression.astype(float)

    tissue.costs['markers'] = cdist(cell_expression, atlas_expression)
    tissue.markers_to_use = markers_to_use
    tissue.num_markers = len(tissue.markers_to_use)
    return tissue



#tissue.markers_to_use = tissue.dge[:, list(markers_to_use)] / np.amax(tissue.dge[:, list(markers_to_use)])



#%%


#%%



class CrossValidation ():

    def __init__(self, atlas, location_a_priori, alpha, dataset, resultat_novosparc ):
        self.atlas = atlas
        self.Xatlas = atlas.X  #
        self.atlas_gene = atlas.var.index.tolist()  # nom des genes de l'atlas
        self.len_atlas_gene = len(atlas.var.index.tolist())
        self.locations_apriori = location_a_priori
        self.dataset = dataset # données rna seq (nb_cell*nb_genes)
        self.num_genes = dataset.shape[1]
        self.gene_names = dataset.var.index.to_list()
        self.num_locations = len(location_a_priori)
        self.alpha = alpha
        self.range_alpha = [0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9]
        self.result = []
        self.resultat_novosparc = resultat_novosparc

        #array for save all the time we compyte novosparc
        self.tot = []

        #array to save all the sdge from novosparc when we remove one gene in computation
        #3dimension :
        self.novo_array_all = np.zeros(shape=(self.len_atlas_gene, np.shape(resultat_novosparc)[0], len(location_a_priori)))

        ### normalisation des données
        self.Xatlas = normalize(self.Xatlas)
        self.Xdataset = dataset.X
        self.Xdataset = normalize(self.Xdataset)

    #def normaliser(self, dataset):

    def formule_quadratique(self, g , g_tilt, num_locations):
        ### g : np.array() dim nb_locations*1 valeur du genes sur les positions (issue de l'atlas reference)
        ### g_tilt :  np.array() dim nb_locations*1 valeur du genes prédite (après novosparc)
        ### num_locations : nombre de localisation spatiale
        return distance.euclidean(g, g_tilt)**2/num_locations

    def calcul_novosparc(self, dataset, locations_apriori, num_neighbors,
                              alpha_i, epsilon_i,
                              atlas_genes, atlas, num_genes, gene_names, num_locations):


        #gene_i = atlas_genes.index(nom_gene_to_delete)
        tissue = novosparc.cm.Tissue(dataset=dataset, locations=locations_apriori)

        num_neighbors_s = num_neighbors_t = num_neighbors

        markers = atlas_genes
        #gene_enleve = atlas_genes[gene_i]
        atlas_matrix = atlas.to_df()[markers].values
        markers_idx = pd.DataFrame({'markers_idx': np.arange(num_genes)}, index=gene_names)
        markers_to_use = np.concatenate(markers_idx.loc[markers].values)

        tissue = update_atlas_variable(atlas_matrix, tissue, markers_to_use)

        tissue.setup_smooth_costs( num_neighbors_s=num_neighbors_s, num_neighbors_t=num_neighbors_t,
        verbose=False)
        #tissue.setup_smooth_costs(num_neighbors_s=num_neighbors_s, num_neighbors_t=num_neighbors_t)
        alpha_linear = alpha_i
        epsilon = epsilon_i

        tissue.reconstruct(alpha_linear=alpha_linear, epsilon=epsilon)

        sdge = tissue.sdge
        #dataset_reconst = sc.AnnData(pd.DataFrame(sdge.T, columns=self.gene_names))

        #dataset_reconst.obsm['spatial'] = locations_apriori

        sdge = sdge.T

        sdge = sdge * num_locations

        #self.novo_array_all = sdge
        return sdge

    def calcul_novosparc_tilt(self, dataset, nom_gene_to_delete,locations_apriori, num_neighbors,
                              alpha_i, epsilon_i,
                              atlas_genes, atlas, num_genes, gene_names, num_locations, Xatlas  ):
        ### calculer novosparc en prenant 1 gene marqueur de moins dans l'atlas de reference

        gene_i = atlas_genes.index(nom_gene_to_delete)
        tissue = novosparc.cm.Tissue(dataset=dataset, locations=locations_apriori)

        num_neighbors_s = num_neighbors_t = num_neighbors

        new_markers = atlas_genes[:gene_i] + atlas_genes[gene_i + 1:]


        gene_enleve = atlas_genes[gene_i]
        atlas_matrix = atlas.to_df()[new_markers].values
        markers_idx = pd.DataFrame({'markers_idx': np.arange(num_genes)}, index=gene_names)
        markers_to_use = np.concatenate(markers_idx.loc[new_markers].values)

        sdge = self.calcul_novosparc(dataset, locations_apriori,
                                     num_neighbors, alpha_i, epsilon_i,
                                     new_markers, atlas, num_genes,
                                     gene_names, num_locations)
        '''
        tissue = update_atlas_variable(atlas_matrix, tissue, markers_to_use)

        tissue.setup_smooth_costs( num_neighbors_s=num_neighbors_s, num_neighbors_t=num_neighbors_t,
        verbose=False)
        #tissue.setup_smooth_costs(num_neighbors_s=num_neighbors_s, num_neighbors_t=num_neighbors_t)
        alpha_linear = alpha_i
        epsilon = epsilon_i

        tissue.reconstruct(alpha_linear=alpha_linear, epsilon=epsilon)

        sdge = tissue.sdge
        #dataset_reconst = sc.AnnData(pd.DataFrame(sdge.T, columns=self.gene_names))

        #dataset_reconst.obsm['spatial'] = locations_apriori
        
        sdge = sdge.T

        sdge = sdge * num_locations
        
        '''
        gw = tissue.gw
        self.novo_array_all[gene_i,:,:] = gw
        
        sdge_not_normalize = sdge


        sdge = normalize(sdge)
        sdge = pd.DataFrame(sdge, columns=gene_names)
        sum = 0
        for j in range(num_locations):
            #if sdge[gene_enleve][j] != 0:
            sum += (sdge[gene_enleve][j] - Xatlas[j][gene_i])**2   #/ (sdge[gene_enleve][j])) ** 2
                # print(sum)
                # print((i+j)/(num_locations+97)*100,'%')

        tot = sum / num_locations
        #self.tot.append(tot)
        return tot, sdge_not_normalize


    def get_sdge_from_gw(self, gw, atlas_genes, nom_gene_to_delete ):
        gene_i = atlas_genes.index(nom_gene_to_delete)
        self.gw_all[gene_i,:,:] = gw


    def calcul_cross_validation_score(self, dataset, liste_gene_loop,locations_apriori, num_neighbors,
                              alpha_i, epsilon_i,
                              atlas_genes, atlas, num_genes, gene_names, num_locations, Xatlas ):

        Q = []
        for gene_to_delete in liste_gene_loop :
            result, sd = self.calcul_novosparc_tilt(dataset, gene_to_delete, locations_apriori, num_neighbors,
                                                alpha_i, epsilon_i,
                                                atlas_genes, atlas, num_genes, gene_names, num_locations, Xatlas)
            Q.append(result)

        self.liste_qi = Q
        self.Q = sum(Q)/len(liste_gene_loop)
        print(Q)
        return Q, sum(Q)/len(liste_gene_loop)


    def histogramme_qi(self, liste_qi):

        fig = px.histogram(liste_qi)
        return fig

'''    def visualisation_comparaison(self, novo_array_all, gene_i_to_remove, atlas):

        try novo_array_all[gene_i_to_remove] is not None :
            sdge = novo_array_all[gene_i_to_remove]
        except :
            sdge =  calcul_novosparc()

'''


#%%









































''' 
    def calcul_qi(self, nom_gene):
        




Xatlas=atlas.X


alpha=[0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9]
result=[]
for a in alpha:
    tot = 0
    for i in range(len(atlas_genes)):
        tissue = novosparc.cm.Tissue(dataset=dataset, locations=locations_apriori)

        num_neighbors_s = num_neighbors_t = 10

        markers = atlas_genes[:i]+atlas_genes[i+1:]
        gene_enleve=atlas_genes[i]
        atlas_matrix = atlas.to_df()[markers].values
        markers_idx = pd.DataFrame({'markers_idx': np.arange(num_genes)}, index=gene_names)
        markers_to_use = np.concatenate(markers_idx.loc[markers].values)
        tissue.setup_reconstruction(atlas_matrix=atlas_matrix,markers_to_use=markers_to_use,num_neighbors_s=num_neighbors_s,num_neighbors_t=num_neighbors_t)
        tissue.setup_reconstruction(num_neighbors_s=num_neighbors_s,num_neighbors_t=num_neighbors_t)
        alpha_linear = a
        epsilon = 1e-2

        tissue.reconstruct(alpha_linear=alpha_linear, epsilon=epsilon)


        sdge = tissue.sdge
        dataset_reconst = sc.AnnData(pd.DataFrame(sdge.T, columns=gene_names))

        dataset_reconst.obsm['spatial'] = locations_apriori

        sdge=sdge.T

        sdge = sdge *360
        sdge=pd.DataFrame(sdge, columns=gene_names)
        sum=0
        for j in range(num_locations):
            if sdge[gene_enleve][j]+Xatlas[j][i]!=0:
                sum+=((sdge[gene_enleve][j]-Xatlas[j][i])/(sdge[gene_enleve][j]))**2
            #print(sum)
            #print((i+j)/(num_locations+97)*100,'%')

        tot+=sum/num_locations
    tot=tot/97
    result.append(tot)

for i in range(len(result)):
    print("alpha =",alpha[i],"quadratique=",result[i])
    
'''