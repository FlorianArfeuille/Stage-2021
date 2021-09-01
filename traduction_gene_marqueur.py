__author__ = 'florian'
__Filename__ = 'traduction gene marqueur'
__Creationdate__ = '06/07/2021'
import numpy as np
import pandas as pd

def position_cell_excel(t, file_name, sheet='Feuil1'):
    df_space = pd.read_excel(io=file_name, sheet_name=sheet)
    df_space = df_space.dropna()
    position = df_space[df_space['time'] == t]
    return position

def gene_cell(t, file_name_gene, file_name_position, sheet_position='Feuil1'):
    position=position_cell_excel(t, file_name_position, sheet_position)
    df = pd.read_csv(file_name_gene)
    df=df.rename(columns = {'name': 'cell_name'})
    name_cell=position['cell_name']
    gene_mat=pd.merge(df,name_cell)
    return gene_mat

def traduction_gene_marqueur(fichier, time_dataset, fichier_dataset, fichier_position):
    dataset=gene_cell(time_dataset,fichier_dataset, fichier_position)
    nom_gene_dataset= dataset.columns[2:]

    file_name=fichier
    df_space = pd.read_csv(file_name)
    del df_space['size']
    nom=df_space.columns
    nom_gene=nom[4:]

    dictio={'fkh-4_11D3b24': 'fkh-4', 'nob-1_L2': 'nob-1', 'pha-4_3E3C5_1xx': 'pha-4', 'mnm-2_7_L1': 'mnm-2', 'hnd-1_6_L1': 'hnd-1', 'nob-1_b26_L1': 'nob-1', 'elt-2_7E1_5_L1': 'elt-2', 'egl-5_7E3_3_L1': 'egl-5', 'skr-8_1_L1': 'skr-8', 'glp-1_4_L1': 'glp-1', 'pgp-3_1': 'pgp-3', 'pha-4_3E3C5_1yy_L1': 'pha-4', 'die-1_3': 'die-1', 'eor-1_9G1_13_L2': 'eor-1', 'fkh-4_11D3b24_L1': 'fkh-4', 'elt-6_8': 'elt-6', 'pal-1_7A2_9_L1': 'pal-1', 'elt-6_RW10178_L2': 'elt-6', 'alr-1_10A2_3_L2': 'alr-1', 'dpl-1_8_L2': 'dpl-1', 'die-1_3e4_10': 'die-1', 'isw-1_3_L1': 'isw-1', 'cnd-1_3C3_6': 'cnd-1', 'ref-2a_14_L2': 'ref-2', 'sea-1_5_L1': 'sea-1', 'egl-5_7E3_12_L1': 'egl-5', 'med-2_2': 'med-2', 'irx-1_b1_L1': 'irx-1', 'tbx-8_10': 'tbx-8', 'tbx-11_10': 'tbx-11', 'cnd-1_3C3_11': 'cnd-1', 'pal-1_7A2_2_L2': 'pal-1', 'pha-4_3E3C5_1': 'pha-4', 'ttx-3_8G1_14_L1': 'ttx-3', 'ref-1_2_L1': 'ref-1', 'elt-6_RW10178_L1': 'elt-6', 'mab-5_3F4_10_L1': 'mab-5', 'tbx-8_6': 'tbx-8', 'tlp-1_RW10609_L1': 'tlp-1', 'mab-5_3F4_10_L2': 'mab-5', 'hnd-1_6': 'hnd-1', 'end-3_2': 'end-3', 'pha-4_3E3C5_1_L2': 'pha-4', 'lir-2_OP175_L1': 'lir-2', 'die-1_sd1566': 'die-1', 'moe-3_8_L1': 'moe-3', 'med-2_6B3_1_L2': 'med-2', 'tbx-11_RW10249_L1': 'tbx-11', 'tbx-8_7G2_1_L2': 'tbx-8', 'ama-1_3A3_5': 'ama-1', 'hnd-1_F396': 'hnd-1', 'tps-2_12_L1': 'tps-2', 'elt-6_5': 'elt-6', 'tps-2_2_L1': 'tps-2', 'end-1_1': 'end-1', 'F47H4_1_14_L1': 'F47H4.2', 'dmd-4b_L1': 'dmd-4', 'sdc-2_3_L1': 'sdc-2', 'glp-1_5_L2': 'glp-1', 'cwn-1_5': 'cwn-1', 'fkh-4_11D3b24_L2': 'fkh-4', 'vab-7_xl1': 'vab-7', 'tbx-11_8': 'tbx-11', 'eor-1_9G1_13_L1': 'eor-1', 'alr-1_10A2_3_L1': 'alr-1', 'mml-1_5': 'mml-1', 'dpy-7_1': 'dpy-7', 'elt-2_sk': 'elt-2', 'tbx-11_RW10249_L2': 'tbx-11', 'pal-1_3': 'pal-1', 'pha-4_3e3_9': 'pha-4', 'cnd-1_3C3_6yy_L1': 'cnd-1', 'cnd-1_3C3_6yy_L2': 'cnd-1', 'dve-1_15_L1': 'dve-1', 'hnd-1_6_L2': 'hnd-1', 'cep-1_b_1_L1': 'cep-1', 'ref-1_2_L2': 'ref-1', 'F47H4_1_14_L2': 'F47H4.2', 'sma-9k_1_L2': 'sma-9', 'pgp-2_5_L2': 'pgp-2', 'pha-4_3E3C5_1yy_L2': 'pha-4', 'mep-1_11C1_8_L1': 'mep-1', 'vab-7_6': 'vab-7', 'tbx-9_8': 'tbx-9', 'lin-39_3d3_6': 'lin-39', 'dpy-31_3_L2': 'dpy-31', 'tbx-38_3': 'tbx-38', 'tbx-35_2x6': 'tbx-35', 'hlh-16_12_L1': 'hlh-16', 'mir_57': 'miro-1', 'nhr-57_9_L1': 'nhr-57', 'egl-27b_1_L1': 'egl-27', 'lin-13_6B12_1_L1': 'lin-13', 'tbx-37_b12_L1': 'tbx-37', 'lin-39_9': 'lin-39', 'sdz-28_4': 'sdz-28', 'ceh-32_4_L1': 'ceh-32', 'ceh-21_6_L2': 'ceh-21', 'lin-32_1_L1': 'lin-32', 'nhr-49b_5_L1': 'nhr-49', 'ceh-16_b7_L1': 'ceh-16', 'ceh-21_5_L1': 'ceh-21', 'ceh-43_RW10154_L1': 'ceh-43', 'sdz-28_6': 'sdz-28', 'sdz-38_1_L1': 'sdz-38', 'lin-13_6B12_1_L2': 'lin-13', 'tbx-38_9': 'tbx-38', 'lin-39_10': 'lin-39', 'ceh-27_1_L1': 'ceh-27', 'rad-26_6_L2': 'rad-26', 'egl-27_b': 'egl-27', 'nhr-79_1_L1': 'nhr-79', 'hlh-26_9': 'hlh-26', 'ceh-36b4_L1': 'ceh-36', 'ceh-14_8F2_14_L1': 'ceh-14', 'nhr-67_3_L2': 'nhr-67', 'ceh-36_b13_L1': 'ceh-36', 'hmg-11_7_L1': 'hmg-11', 'mel-28_2_L1': 'mel-28', 'cnd1_4-2': 'cnd-1', 'eft3': 'eftu-2', 'nhr-57_8_L1': 'nhr-57', 'ceh-36_RW10598_L1': 'ceh-36', 'cnd1_3': 'cnd-1', 'tbx-35_6': 'tbx-35', 'ceh-41_2_L1': 'ceh-41', 'nhr-57_5_L1': 'nhr-57', 'egl-27c_L1': 'egl-27', 'C05D10_1b_3_L1': 'C05D10.4', 'ceh-43_RW10345_L1': 'ceh-43', 'mml1_b': 'mml-1', 'nhr-68_13_L2': 'nhr-68', 'ceh-43_RW10345_L2': 'ceh-43', 'T23H4_2_2_L1': 'T23H2.3', 'nhr-68_14_L1': 'nhr-68', 'nhr-67_3_L1': 'nhr-67', 'ceh-36B4_L1': 'ceh-36', 'ceh-41_14_L1': 'ceh-41', 'ceh-32_4_L2': 'ceh-32', 'nhr-79_3_L2': 'nhr-79', 'mir61_e16': 'miro-1', 'lin-26_1_L1': 'lin-26', 'ceh-432x3': 'ceh-43', 'T23H4_2_10_L1': 'T23H2.3', 'tbx-37_b12_L2': 'tbx-37', '_nhr-23_gsIs269_L1': 'gene_short_name', 'F38C2_7_12_L2': 'F38C2.7', 'F21A10_2_9_L1': 'F21A9.1', 'nhr-25_3H4_8': 'nhr-25', 'pax-3_3_L1': 'pax-3', 'B0310_2_8_L1': 'B0310.3', 'nhr-2_9F1_10': 'nhr-2', 'F28C6_1_10_L2': 'F28C6.10', 'T28H10_3_4_L2': 'T28H10.3', 'W10D9_4_4_L1': 'W10D9.1', 'hlh-1_6B4_3_L1': 'hlh-1', 'T28H10_3_4_L1': 'T28H10.3', 'ceh-26_10b1_3_L1': 'ceh-26', 'F21A10_2a_10_L1': 'F21A9.2', 'F58D2_1_2_L1': 'F58D2.2', 'D1081_8_L2': 'D1081.7', 'hlh_3_5': 'hlh-3', 'nhr-2_1': 'nhr-2', 'F39B2_1_7_L1': 'F39B2.3', 'R144_3_3_RW10881_L1': 'R144.3', 'B0310_2_13_L1': 'B0310.2', 'dsl-1_5': 'dsl-1', 'T22C8_3_1_L2': 'T22C8.3', 'T23G4_1_1_L1': 'T23G4.2', 'unc-130_9E1_4_L1': 'unc-130', 'hlh1_4': 'hlh-1', 'F21D5_9_1_L1': 'F21D5.9', 'nhr-2_9F1_10_L1': 'nhr-2', 'nhr-25_h4_8': 'nhr-25', 'F09G2_9_1_L1': 'F09G2.8', 'K02G10_1_4_L1': 'K02G10.15', 'B0310_2_8_L2': 'B0310.2', 'tag-185_b_6_L1': 'tag-185', 'Y106G6H_4_9_L2': 'Y106G6H.1', 'elt-1_3': 'elt-1', 'F23F12_9_3_L1': 'F23F12.8', 'nhr-171_17_L1': 'nhr-171', 'nhr-171_16_L2': 'nhr-171', 'F28C6_1_2_L1': 'F28C6.10', 'C08B11_3_6': 'C08B11.8', 'ceh-6_9H1_B_L1': 'ceh-6', 'hlh-1_4': 'hlh-1', 'F09G2_9_6_L1': 'F09G2.8', 'C01B7_1_6': 'C01B7.3', 'F16B12_6_2_L1': 'W10D9.1', 'F21D5_9_5_L1': 'F21D5.9', 'T23G4_1_1_L2': 'T23G4.2', 'F17C11_5_12_L1': 'F17C11.4', 'C08B11_3b': 'C08B11.9', 'pes-1_7_L1': 'pes-1', 'pes-1_7_L2': 'pes-1', 'T23G5_6_15_L1': 'T23G5.2', 'C25D7_10_1_L2': 'C25D7.10', 'D1081_8_L1': 'D1081.9', 'unc-130_9E1_8_L2': 'unc-130', 'R144_3_12_L1': 'R144.3', 'R144_3_12_L2': 'R144.3', 'T22C8_3_11_L1': 'T22C8.3', 'ZK185_1_3_L1': 'ZK185.1', 'egl5_sk': 'egl-5', 'hlh26_1_L2': 'hlh-26', 'W10D9_4_4_L2': 'W10D9.3', 'B0310_2_2_L1': 'B0310.2', 'T23G5_6_15_L2': 'T23G5.3', 'pax-3_1_L1': 'pax-3', 'F21D5_9_9_L2': 'F21D5.9', 'B0336_3_13_L1': 'B0336.3', 'F16B12_6_1_L1': 'F16B12.6', 'pha4_b2': 'pha-4', 'hsp3_a1a': 'hsp-3', 'pha4I2L_6': 'pha-4', 'pha4_I2L_3': 'pha-4', 'his72_D1A': 'his-72', 'elt7_c2a': 'elt-7', 'pha4I2L_11': 'pha-4', 'dyf7pJIM20': 'dyf-7', 'end3_2': 'end-3', 'lin1_10': 'lin-1', 'end3_2xx': 'end-3', 'pha4A1L_20_10026': 'pha-4', 'pha4_I1LBBB': 'pha-4', 'ob-1GFP': 'gob-1', 'lh1AF16': 'hlh-1', 'rw10029_norfpxx': 'srw-100', 'rw10029_norfp': 'srw-100', 'ob1GFP': 'gob-1', 'rw10026_hlh1_8': 'hlh-1', 'mi57': 'mig-17', 'rw10029_norfpyy': 'srw-100', 'noRFP': 'nog-1', 'Ges1-4': 'ges-1', '_tbx-8_7G2_1': 'tbx-8', '_lin-11_7H3_1_L1': 'lin-11', 'GL-5_AF16_L1': 'egl-5', 'RW10029_cnd1_9': 'cnd-1', 'c50f7_5_7': 'C50F7.5', 'f7red': 'dyf-7'}


    cle_sup=[]
    cle1=[]
    for e in nom_gene:
        for b in nom_gene:
            if e in dictio.keys() and b in dictio.keys():
                if b!=e and dictio[e]==dictio[b] and b not in cle1:
                    cle_sup.append(b)
                    cle1.append(e)
    del dictio['B0310_2_13_L1']
    del dictio['R144_3_3_RW10881_L1']
    del dictio['cnd-1_3C3_6yy_L1']
    del dictio['cnd-1_3C3_6yy_L2']
    del dictio['elt-6_8']
    del dictio['elt-6_RW10178_L2']
    del dictio['hlh-1_6B4_3_L1']
    del dictio['lh1AF16']
    del dictio['hlh1_4']
    del dictio['nhr-57_5_L1']
    del dictio['die-1_sd1566']
    del dictio['pha-4_3E3C5_1yy_L2']
    del dictio['pha4I2L_11']
    del dictio['pha4I2L_6']
    del dictio['pha-4_3e3_9']
    del dictio['pha-4_3E3C5_1']
    del dictio['pha4_I1LBBB']
    del dictio['hnd-1_F396']
    del dictio['hnd-1_6_L2']
    del dictio['tbx-8_6']
    del dictio['tbx-11_8']
    del dictio['tbx-11_10']
    del dictio['ceh-43_RW10345_L1']
    del dictio['ceh-43_RW10154_L1']
    del dictio['pal-1_7A2_2_L2']
    del dictio['cnd-1_3C3_11']
    del dictio['cnd1_4-2']
    del dictio['F21D5_9_1_L1']
    del dictio['R144_3_12_L1']
    del dictio['ceh-36_b13_L1']
    del dictio['ceh-36_RW10598_L1']
    del dictio['rw10029_norfpxx']
    for i in cle_sup:
        if i!='rw10029_norfpxx' and i!='ceh-36_RW10598_L1' and i!='ceh-36_b13_L1' and i!='R144_3_12_L1' and i!='F21D5_9_1_L1' and i!='cnd1_4-2' and i!='cnd-1_3C3_11' and i!='pal-1_7A2_2_L2' and i!='ceh-43_RW10154_L1' and i!='ceh-43_RW10345_L1' and i!='tbx-11_10' and i!='tbx-11_8' and i!='tbx-8_6' and i!='hnd-1_6_L2' and i!='hnd-1_F396' and i!='pha4_I1LBBB' and i!= 'pha-4_3E3C5_1' and i!='pha-4_3e3_9' and i!='pha4I2L_6' and i!='pha4I2L_11' and i!='pha-4_3E3C5_1yy_L2' and i!='die-1_sd1566' and i!='nhr-57_5_L1' and i!='hlh1_4' and i !='B0310_2_13_L1' and i!='R144_3_3_RW10881_L1' and i!='cnd-1_3C3_6yy_L1' and i!='cnd-1_3C3_6yy_L2' and i!='cnd1_3' and i!='elt-6_8' and i!='elt-6_RW10178_L2' and i!='hlh-1_6B4_3_L1' and i!='lh1AF16':
            del dictio[i]

    for e in nom_gene:
        if e in dictio.keys():
            df_space=df_space.rename(columns = {e : dictio[e]})
            if dictio[e] not in nom_gene_dataset:
                del df_space[dictio[e]] # repaire les genes qui ne sont pas dans le fichier transcriptomie_3000
        else:
            del df_space[e]#suppresion des genes qui ne sont pas présent dans le dictionnaire

    nom=df_space.columns # actualise avec les nouveaux nom de genes présent
    nom_gene=nom[5:]
    df_space[nom_gene].to_csv("gene_marker_news.txt", sep=' ',index=False)
