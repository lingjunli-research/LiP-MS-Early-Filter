# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:04:18 2022

@author: lawashburn
"""

import pandas as pd
import numpy as np
import csv
from scipy import stats
import math

working_directory = r"C:\Users\lawashburn\Documents\LiP-MS\EarlyFilter\Working_Directory"
#protein_path = input('Input path to Proteingroup.csv file: ')
#peptides_path = r"C:\Users\lawashburn\Documents\LiP-MS\EarlyFilter\15_peptidesgroup.csv"
#protein_path = r"C:\Users\lawashburn\Documents\LiP-MS\EarlyFilter\15_Proteingroup.csv"
#peptides_path = r"C:\Users\lawashburn\Documents\LiP-MS\15_peptidesgroup.csv"
#working_directory = r"C:\Users\lawashburn\Documents\LiP-MS\Out_update"
updated_common = r"C:\Users\lawashburn\Documents\LiP-MS\LiP-MS_EarlyFilter_Distribution\Example_Data\common_updated.csv"

df_common = pd.read_csv(updated_common, na_values= '#DIV/0!' )

df_common.dropna(subset = ["Protein name"], inplace=True)

df_common = df_common.rename(columns={'LFQ intensity CTRL_351':'Pep Ctrl 1','LFQ intensity CTRL_375':'Pep Ctrl 2','LFQ intensity CTRL_275':'Pep Ctrl 3',
                                  'LFQ intensity CTRL_019':'Pep Ctrl 4','LFQ intensity CTRL_111':'Pep Ctrl 5','LFQ intensity AD_059':'Pep AD 1',
                                  'LFQ intensity AD_214':'Pep AD 2','LFQ intensity AD_075':'Pep AD 3','LFQ intensity AD_082':'Pep AD 4',
                                  'LFQ intensity AD_078':'Pep AD 5','LFQ intensity MCI_154':'Pep MCI 1','LFQ intensity MCI_068':'Pep MCI 2',
                                  'LFQ intensity MCI_219':'Pep MCI 3','LFQ intensity MCI_242':'Pep MCI 4','LFQ intensity MCI_003':'Pep MCI 5'})

print(df_common)
df_common = df_common.drop(df_common[df_common['Pep Ctrl 1'] == 'CTRL'].index)
df_common = df_common.drop(df_common[df_common['Pep AD 1'] == 'AD'].index)
df_common = df_common.drop(df_common[df_common['Pep MCI 1'] == 'MCI'].index)

ad_ctrl_test = df_common["AD/Ctrl fold change"].astype(float)
ad_ctrl_p = df_common['p-value_ctrl_ad'].astype(float)
ad_ctrl_test = np.where(ad_ctrl_test > 1.5, df_common["AD/Ctrl fold change"],1)
ad_ctrl_test = np.where(ad_ctrl_test < (2/3), df_common["AD/Ctrl fold change"],1)
ad_ctrl_test = np.where(ad_ctrl_p < 0.05, df_common["AD/Ctrl fold change"],1)

ADvCtrl = pd.DataFrame()
ADvCtrl['AD 1 vs Ctrl normalized'] = (df_common['Pep AD 1'].astype(float))/(ad_ctrl_test)
ADvCtrl['AD 2 vs Ctrl normalized'] = (df_common['Pep AD 2'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 3 vs Ctrl normalized'] = (df_common['Pep AD 3'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 4 vs Ctrl normalized'] = (df_common['Pep AD 4'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 5 vs Ctrl normalized'] = (df_common['Pep AD 5'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD vs Ctrl normalized average'] = ADvCtrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Pep Ctrl 4'] = df_common['Pep Ctrl 4'].astype(float)
Ctrl['Pep Ctrl 5'] = df_common['Pep Ctrl 5'].astype(float)
Ctrl['Ctrl mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
ADvCtrl['Pep Ctrl 1'] = Ctrl['Pep Ctrl 1']
ADvCtrl['Pep Ctrl 2'] = Ctrl['Pep Ctrl 2']
ADvCtrl['Pep Ctrl 3'] = Ctrl['Pep Ctrl 3']
ADvCtrl['Pep Ctrl 4'] = Ctrl['Pep Ctrl 4']
ADvCtrl['Pep Ctrl 5'] = Ctrl['Pep Ctrl 5']
ADvCtrl['Ctrl mean'] = Ctrl['Ctrl mean']

ad_mci_test = df_common["AD/MCI fold change"].astype(float)
ad_mci_p = df_common['p-value_ad_mci'].astype(float)
ad_mci_test = np.where(ad_mci_test > 1.5, df_common["AD/MCI fold change"],1)
ad_mci_test = np.where(ad_mci_test < (2/3), df_common["AD/MCI fold change"],1)
ad_mci_test = np.where(ad_mci_p < 0.05, df_common["AD/MCI fold change"],1)

ADvMCI = pd.DataFrame()
ADvMCI['AD 1 vs MCI normalized'] = (df_common['Pep AD 1'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 2 vs MCI normalized'] = (df_common['Pep AD 2'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 3 vs MCI normalized'] = (df_common['Pep AD 3'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 4 vs MCI normalized'] = (df_common['Pep AD 4'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 5 vs MCI normalized'] = (df_common['Pep AD 5'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD vs MCI normalized average'] = ADvMCI.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
MCI = pd.DataFrame()
MCI['Pep MCI 1'] = df_common['Pep MCI 1'].astype(float)
MCI['Pep MCI 2'] = df_common['Pep MCI 2'].astype(float)
MCI['Pep MCI 3'] = df_common['Pep MCI 3'].astype(float)
MCI['Pep MCI 4'] = df_common['Pep MCI 4'].astype(float)
MCI['Pep MCI 5'] = df_common['Pep MCI 5'].astype(float)
MCI['MCI mean'] = MCI.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
ADvMCI['Pep MCI 1'] = MCI['Pep MCI 1']
ADvMCI['Pep MCI 2'] = MCI['Pep MCI 2']
ADvMCI['Pep MCI 3'] = MCI['Pep MCI 3']
ADvMCI['Pep MCI 4'] = MCI['Pep MCI 4']
ADvMCI['Pep MCI 5'] = MCI['Pep MCI 5']
ADvMCI['MCI mean'] = MCI['MCI mean']

mci_ctrl_test = df_common["MCI/Ctrl fold change"].astype(float)
mci_ctrl_test = df_common['p-value_ctrl_mci'].astype(float)
mci_ctrl_test = np.where(ad_mci_test > 1.5, df_common["MCI/Ctrl fold change"],1)
mci_ctrl_test = np.where(ad_mci_test < (2/3), df_common["MCI/Ctrl fold change"],1)
mci_ctrl_test = np.where(ad_mci_p < 0.05, df_common["MCI/Ctrl fold change"],1)

MCIvCtrl = pd.DataFrame()
MCIvCtrl['MCI 1 vs Ctrl normalized'] = (df_common['Pep MCI 1'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 2 vs Ctrl normalized'] = (df_common['Pep MCI 2'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 3 vs Ctrl normalized'] = (df_common['Pep MCI 3'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 4 vs Ctrl normalized'] = (df_common['Pep MCI 4'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 5 vs Ctrl normalized'] = (df_common['Pep MCI 5'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI vs Ctrl normalized average'] = MCIvCtrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Pep Ctrl 4'] = df_common['Pep Ctrl 4'].astype(float)
Ctrl['Pep Ctrl 5'] = df_common['Pep Ctrl 5'].astype(float)
Ctrl['Ctrl mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
MCIvCtrl['Pep Ctrl 1'] = Ctrl['Pep Ctrl 1']
MCIvCtrl['Pep Ctrl 2'] = Ctrl['Pep Ctrl 2']
MCIvCtrl['Pep Ctrl 3'] = Ctrl['Pep Ctrl 3']
MCIvCtrl['Pep Ctrl 4'] = Ctrl['Pep Ctrl 4']
MCIvCtrl['Pep Ctrl 5'] = Ctrl['Pep Ctrl 5']
MCIvCtrl['Ctrl mean'] = Ctrl['Ctrl mean']

Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Pep Ctrl 4'] = df_common['Pep Ctrl 4'].astype(float)
Ctrl['Pep Ctrl 5'] = df_common['Pep Ctrl 5'].astype(float)
Ctrl['Ctrl average'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

ADvCtrl['Protein IDs'] = df_common['Protein ID']
ADvCtrl['Protein Names'] = df_common['Protein name']
ADvCtrl['Gene Names'] = df_common['Gene names']
ADvCtrl['Fasta headers'] = df_common['Fasta header']
ADvCtrl['Sequence'] = df_common['Sequence']

ADvMCI['Protein IDs'] = df_common['Protein ID']
ADvMCI['Protein Names'] = df_common['Protein name']
ADvMCI['Gene Names'] = df_common['Gene names']
ADvMCI['Fasta headers'] = df_common['Fasta header']
ADvMCI['Sequence'] = df_common['Sequence']

MCIvCtrl['Protein IDs'] = df_common['Protein ID']
MCIvCtrl['Protein Names'] = df_common['Protein name']
MCIvCtrl['Gene Names'] = df_common['Gene names']
MCIvCtrl['Fasta headers'] = df_common['Fasta header']
MCIvCtrl['Sequence'] = df_common['Sequence']

ADvCtrl['AD vs Ctrl fold change'] = ADvCtrl['AD vs Ctrl normalized average'] / Ctrl['Ctrl average']
ADvMCI['AD vs MCI fold change'] = ADvMCI['AD vs MCI normalized average'] / Ctrl['Ctrl average']
MCIvCtrl['MCI vs Ctrl fold change'] = MCIvCtrl['MCI vs Ctrl normalized average'] / Ctrl['Ctrl average']

#ADvCtrl.dropna(subset = ['AD vs Ctrl fold change'], inplace=True)
ADvMCI.dropna(subset = ['AD vs MCI fold change'], inplace=True)
MCIvCtrl.dropna(subset = ['MCI vs Ctrl fold change'], inplace=True)

ADvCtrl_out = working_directory + '\\ADvCtrl_nofilter.csv' 
ADvMCI_out = working_directory + '\\ADvMCI_nofilter.csv' 
MCIvCtrl_out = working_directory + '\\MCIvCtrl_nofilter.csv' 

with open(ADvCtrl_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvCtrl.to_csv(filed,index=False)
    
with open(ADvMCI_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    ADvMCI.to_csv(filed,index=False)

with open(MCIvCtrl_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    MCIvCtrl.to_csv(filed,index=False)

    
print('Unfiltered data has been exported to working directory')