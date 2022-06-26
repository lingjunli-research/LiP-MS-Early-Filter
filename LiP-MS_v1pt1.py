# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:48:03 2022

@author: lawashburn
"""

import pandas as pd
import numpy as np
import csv
from scipy import stats
import math

#user input prompts at command line
working_directory = input('Input path to working directory, where all results will be stored: ')
protein_path = input('Input path to Proteingroup.csv file: ')
peptides_path = input('Input path to Peptidegroup.csv file: ')

protein = pd.read_csv(protein_path)  #read Proteingroup.csv
peptide = pd.read_csv(peptides_path) #read Peptidegroup.csv

#Renaming columns for easier handling
protein = protein.rename(columns={'LFQ intensity MCI_003':'MCI1','LFQ intensity MCI_154':'MCI2','LFQ intensity MCI_242':'MCI3',
                                  'LFQ intensity MCI_219':'MCI4','LFQ intensity MCI_068':'MCI5','LFQ intensity CTRL_019':'Ctrl1',
                                  'LFQ intensity CTRL_275':'Ctrl2','LFQ intensity CTRL_375':'Ctrl3','LFQ intensity CTRL_351':'Ctrl4',
                                  'LFQ intensity CTRL_111':'Ctrl5','LFQ intensity AD_214':'AD1','LFQ intensity AD_082':'AD2',
                                  'LFQ intensity AD_075':'AD3','LFQ intensity AD_078':'AD4','LFQ intensity AD_059':'AD5'})

peptide = peptide.rename(columns={'Protein names':'Protein name'})

#DataFrame of AD protein values
AD = pd.DataFrame()
AD['AD1'] = protein['AD1']
AD['AD2'] = protein['AD2']
AD['AD3'] = protein['AD3']
AD['AD4'] = protein['AD4']
AD['AD5'] = protein['AD5']
AD = AD.iloc[1:,:] #remove first row of dataframe
AD['mean'] = AD.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)  #calculate mean intensity of all AD protein replicates

#DataFrame of MCI protein values
MCI = pd.DataFrame()
MCI['MCI1'] = protein['MCI1']
MCI['MCI2'] = protein['MCI2']
MCI['MCI3'] = protein['MCI3']
MCI['MCI4'] = protein['MCI4']
MCI['MCI5'] = protein['MCI5']
MCI = MCI.iloc[1:,:]  #remove first row of dataframe
MCI['mean'] = MCI.astype(float).replace(0, np.nan).mean(axis=1, skipna=True) #calculate mean intensity of all MCI protein replicates

#DataFrame of Ctrl protein values
Ctrl = pd.DataFrame()
Ctrl['Ctrl1'] = protein['Ctrl1']
Ctrl['Ctrl2'] = protein['Ctrl2']
Ctrl['Ctrl3'] = protein['Ctrl3']
Ctrl['Ctrl4'] = protein['Ctrl4']
Ctrl = Ctrl.iloc[1:,:] #remove first row of dataframe
Ctrl['mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True) #calculate mean intensity of all Ctrl protein replicates

#Generate DataFrame to store all calculations
protein_calc = pd.DataFrame()
protein_calc['Protein ID'] = protein['Protein IDs']
protein_calc['Protein name'] = protein['Protein names']
protein_calc['Gene name'] = protein['Gene names']
protein_calc['Fasta header'] = protein['Fasta headers']
protein_calc['AD avg'] = AD['mean']
protein_calc['MCI avg'] = MCI['mean']
protein_calc['Ctrl avg'] = Ctrl['mean']
protein_calc['AD/Ctrl fold change'] = protein_calc['AD avg'] / protein_calc['Ctrl avg'] #ADvCtrl fold change
protein_calc['AD/MCI fold change'] = protein_calc['AD avg'] / protein_calc['MCI avg'] #ADvMCI fold change
protein_calc['MCI/Ctrl fold change'] = protein_calc['MCI avg'] / protein_calc['Ctrl avg'] #MCIvCtrl fold  change

df_common = pd.merge(protein_calc, peptide, on=['Protein name'], how='inner') #matching each peptide to the protein

#Export results
MCIvCtrl_out2 = working_directory + '\\common.csv'
with open(MCIvCtrl_out2,'w',newline='') as filed:
    writerd = csv.writer(filed)
    df_common.to_csv(filed,index=False)

