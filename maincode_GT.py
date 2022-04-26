# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 15:23:45 2022

@author: Ehsan Saghapour
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

## Pre-Processing data

# Read TCGA dataset
TCGA_Data = pd.read_csv('TCGA1.csv')
idx_sample=TCGA_Data.sampleID
TCGA_Data=TCGA_Data.drop('sampleID', axis=1)
idx_Genes=TCGA_Data.columns
TCGA_Data=pd.DataFrame(TCGA_Data.T.values,columns=idx_sample)
TCGA_Data.index=idx_Genes

# Read Phenotype dataset
TCGA_clinical=pd.read_csv('TCGA_clinical_metadata.csv')

# Read Layout 
organic_layout= pd.read_csv('Layout.txt',header=None,sep='	')
organic_layout.index= organic_layout.iloc[:,0]
organic_layout=organic_layout.drop(0, axis=1)
organic_layout.columns=['X','Y','PPI']

# Merge TCGA Data and Layout
Merge_dataset=pd.merge(organic_layout, TCGA_Data, left_index=True, right_index=True)

standard_embedding= Merge_dataset[['X','Y']]
TCGA_Data=Merge_dataset[idx_sample]

# Create GeneTerrain for each patient
res=200
x_ = np.linspace(0, 1, res)
y_ = np.linspace(0, 1, res)
X, Y = np.meshgrid(x_, y_)
pi=3.14

for name_sample in TCGA_Data.columns:
      gaussian=np.zeros((res,res))

      
      for i in range(np.shape(standard_embedding)[0]):
          sigma= .05
          x=standard_embedding['X'][i]
          y=standard_embedding['Y'][i]
          amp1=TCGA_Data[name_sample][i]
          
          gaussian1 =np.exp(-(((X-x)/sigma)**2+((Y-y)/sigma)**2)) 
          gaussian =amp1* gaussian1 + gaussian

      plt.figure(figsize=(7, 6))
      norm=plt.Normalize(-10,10)
      plt.pcolormesh(gaussian, cmap = 'jet',norm=norm)
      plt.colorbar()
      idxxs=np.where(TCGA_clinical['sampleID']==name_sample)[0]
      category=TCGA_clinical['GeneExp_Subtype'][idxxs].tolist()[0]
 
      path=os.getcwd()
      
      if category =='Classical':
            plt.savefig(path+ '\\desktop\\GeneTerain_GBM\\Classical\\'+ name_sample +'.jpeg',dpi=100)           
      elif category=='Mesenchymal':
            plt.savefig(path+'\\desktop\\GeneTerain_GBM\\Mesenchymal\\'+ name_sample +'.jpeg',dpi=100)
      elif category=='Neural':
            plt.savefig(path+'\\desktop\\GeneTerain_GBM\\Neural\\'+ name_sample +'.jpeg', dpi=100)
      else:
            plt.savefig(path+'\\desktop\\GeneTerain_GBM\\Proneural\\'+ name_sample +'.jpeg',dpi=100)
      plt.close()     
    
   
    
   