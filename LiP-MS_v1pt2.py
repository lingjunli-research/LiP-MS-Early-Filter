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

#p-value cutoff
p_cutoff = 0.05
#upregulated fold change cutoff
up_fc_cutoff = 1.5
#downregulated fold change cutoff
down_fc_cutoff = (2/3)


#Path to output directory
working_directory = r"C:\Users\lawashburn\Documents\LiP-MS\EarlyFilter\Working_Directory"
#Import .csv from v1pt1 with p-value calculated
updated_common = r"C:\Users\lawashburn\Documents\LiP-MS\LiP-MS_EarlyFilter_Distribution\Example_Data\common_updated.csv"
#Import results, sub #DIV/0! errors with nan
df_common = pd.read_csv(updated_common, na_values= '#DIV/0!' )
#Drop rows where protein name is empty
df_common.dropna(subset = ["Protein name"], inplace=True)
#reformatting columns
df_common = df_common.rename(columns={'LFQ intensity CTRL_351':'Pep Ctrl 1','LFQ intensity CTRL_375':'Pep Ctrl 2','LFQ intensity CTRL_275':'Pep Ctrl 3',
                                  'LFQ intensity CTRL_019':'Pep Ctrl 4','LFQ intensity CTRL_111':'Pep Ctrl 5','LFQ intensity AD_059':'Pep AD 1',
                                  'LFQ intensity AD_214':'Pep AD 2','LFQ intensity AD_075':'Pep AD 3','LFQ intensity AD_082':'Pep AD 4',
                                  'LFQ intensity AD_078':'Pep AD 5','LFQ intensity MCI_154':'Pep MCI 1','LFQ intensity MCI_068':'Pep MCI 2',
                                  'LFQ intensity MCI_219':'Pep MCI 3','LFQ intensity MCI_242':'Pep MCI 4','LFQ intensity MCI_003':'Pep MCI 5'})

#Drop rows with empty values
df_common = df_common.drop(df_common[df_common['Pep Ctrl 1'] == 'CTRL'].index)
df_common = df_common.drop(df_common[df_common['Pep AD 1'] == 'AD'].index)
df_common = df_common.drop(df_common[df_common['Pep MCI 1'] == 'MCI'].index)

ad_ctrl_test = df_common["AD/Ctrl fold change"].astype(float)
ad_ctrl_p = df_common['p-value_ctrl_ad'].astype(float)

#Identify rows where fold change is greater than 1.5 or less than 2/3, and where p value is less than 0.05
ad_ctrl_test = np.where(ad_ctrl_test > up_fc_cutoff, df_common["AD/Ctrl fold change"],1)
ad_ctrl_test = np.where(ad_ctrl_test < down_fc_cutoff, df_common["AD/Ctrl fold change"],1)
ad_ctrl_test = np.where(ad_ctrl_p < p_cutoff, df_common["AD/Ctrl fold change"],1)

#Normalize all AD:Ctrl replicates
ADvCtrl = pd.DataFrame()
ADvCtrl['AD 1 vs Ctrl normalized'] = (df_common['Pep AD 1'].astype(float))/(ad_ctrl_test)
ADvCtrl['AD 2 vs Ctrl normalized'] = (df_common['Pep AD 2'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 3 vs Ctrl normalized'] = (df_common['Pep AD 3'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 4 vs Ctrl normalized'] = (df_common['Pep AD 4'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
ADvCtrl['AD 5 vs Ctrl normalized'] = (df_common['Pep AD 5'].astype(float)) / (df_common['AD/Ctrl fold change'].astype(float))
#Average of AD:Ctrl normalized replicates
ADvCtrl['AD vs Ctrl normalized average'] = ADvCtrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

#DataFrame of Ctrl replicates
Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Pep Ctrl 4'] = df_common['Pep Ctrl 4'].astype(float)
Ctrl['Pep Ctrl 5'] = df_common['Pep Ctrl 5'].astype(float)
#Mean of Ctrl replicates
Ctrl['Ctrl mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

#Organization of ctrl intensities
ADvCtrl['Pep Ctrl 1'] = Ctrl['Pep Ctrl 1']
ADvCtrl['Pep Ctrl 2'] = Ctrl['Pep Ctrl 2']
ADvCtrl['Pep Ctrl 3'] = Ctrl['Pep Ctrl 3']
ADvCtrl['Pep Ctrl 4'] = Ctrl['Pep Ctrl 4']
ADvCtrl['Pep Ctrl 5'] = Ctrl['Pep Ctrl 5']
ADvCtrl['Ctrl mean'] = Ctrl['Ctrl mean']

ad_mci_test = df_common["AD/MCI fold change"].astype(float)
ad_mci_p = df_common['p-value_ad_mci'].astype(float)

#Filter AD:MCI based on fold change and p-value cutoffs
ad_mci_test = np.where(ad_mci_test > up_fc_cutoff, df_common["AD/MCI fold change"],1)
ad_mci_test = np.where(ad_mci_test < down_fc_cutoff, df_common["AD/MCI fold change"],1)
ad_mci_test = np.where(ad_mci_p < p_cutoff, df_common["AD/MCI fold change"],1)

#Normalize AD:MCI intensity
ADvMCI = pd.DataFrame()
ADvMCI['AD 1 vs MCI normalized'] = (df_common['Pep AD 1'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 2 vs MCI normalized'] = (df_common['Pep AD 2'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 3 vs MCI normalized'] = (df_common['Pep AD 3'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 4 vs MCI normalized'] = (df_common['Pep AD 4'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
ADvMCI['AD 5 vs MCI normalized'] = (df_common['Pep AD 5'].astype(float)) / (df_common['AD/MCI fold change'].astype(float))
#Mean of all AD:MCI normalized intensities
ADvMCI['AD vs MCI normalized average'] = ADvMCI.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

#Organization of MCI intensities
MCI = pd.DataFrame()
MCI['Pep MCI 1'] = df_common['Pep MCI 1'].astype(float)
MCI['Pep MCI 2'] = df_common['Pep MCI 2'].astype(float)
MCI['Pep MCI 3'] = df_common['Pep MCI 3'].astype(float)
MCI['Pep MCI 4'] = df_common['Pep MCI 4'].astype(float)
MCI['Pep MCI 5'] = df_common['Pep MCI 5'].astype(float)
#Mean of MCI intensities
MCI['MCI mean'] = MCI.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
ADvMCI['Pep MCI 1'] = MCI['Pep MCI 1']
ADvMCI['Pep MCI 2'] = MCI['Pep MCI 2']
ADvMCI['Pep MCI 3'] = MCI['Pep MCI 3']
ADvMCI['Pep MCI 4'] = MCI['Pep MCI 4']
ADvMCI['Pep MCI 5'] = MCI['Pep MCI 5']
ADvMCI['MCI mean'] = MCI['MCI mean']

mci_ctrl_test = df_common["MCI/Ctrl fold change"].astype(float)
mci_ctrl_test = df_common['p-value_ctrl_mci'].astype(float)
#Filter MCI:Ctrl entries by fc and p-values
mci_ctrl_test = np.where(ad_mci_test > up_fc_cutoff, df_common["MCI/Ctrl fold change"],1)
mci_ctrl_test = np.where(ad_mci_test < down_fc_cutoff, df_common["MCI/Ctrl fold change"],1)
mci_ctrl_test = np.where(ad_mci_p < p_cutoff, df_common["MCI/Ctrl fold change"],1)

#Normalized intensities of MCI:ctrl intensities
MCIvCtrl = pd.DataFrame()
MCIvCtrl['MCI 1 vs Ctrl normalized'] = (df_common['Pep MCI 1'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 2 vs Ctrl normalized'] = (df_common['Pep MCI 2'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 3 vs Ctrl normalized'] = (df_common['Pep MCI 3'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 4 vs Ctrl normalized'] = (df_common['Pep MCI 4'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
MCIvCtrl['MCI 5 vs Ctrl normalized'] = (df_common['Pep MCI 5'].astype(float)) / (df_common['MCI/Ctrl fold change'].astype(float))
#Mean of MCI:ctrl intensities
MCIvCtrl['MCI vs Ctrl normalized average'] = MCIvCtrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

#Ctrl intensitiy organization
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

#Peptide fold change calculation
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