#!/usr/bin/env python
# coding: utf-8

# # Dataset for subgroup discovery

# In[5]:


from Serotype_Data import * #Dataset of Serotypes
from Serotype_Functions import *
from PairXPair_dataset import *

#Biopython ver.1.76
from Bio.Seq import Seq #Represent biological sequences with alphabets
from Bio.pairwise2 import format_alignment #Functions to get global and local alignments between two sequences
from Bio import pairwise2

from IPython.display import Image

import numpy as np
import pandas as pd
import tabulate


# In[6]:


import json
import os


# In[7]:


path = os.path.relpath('./serotypes.json')
f = open(path, 'r')
serotypes = json.load(f)


# In[8]:


total_pairgenes_df_remove_fre = pd.DataFrame({
    'Pair of genes': li_store_genes, 
    'Pair of sugars': li_store_sugars,
    #'Num. of occurrence of pair of gene (Written in front of each pair of genes)': li_pair_occur,
    'Num. of occurrence where pair of sugars processing (included observation)': li_pair_processing,
    'Num. of occurrence where pair of sugars processing (without observation)': li_pair_processing_without_one,
    'Serotype(s) where processing pair of genes and sugars': li_ser_pair_processing
})


# In[9]:


eli_li_store_genes = []#li_store_genes 
eli_li_store_sugars = [] #li_store_sugars
eli_li_pair_processing = [] #li_pair_processing
eli_li_pair_processing_without_one = [] #li_pair_processing_without_one
eli_li_ser_pair_processing = []

single_index = []
for i in range(len(li_ser_pair_processing)):
    # Store index of len(ser_pair_processing)==1
    if len(li_ser_pair_processing[i]) == 1:
        single_index.append(i) 
    
li_indexs = list(range(len(li_store_genes)))
for index in single_index:
    if index in li_indexs:
        li_indexs.remove(index)

for indexs in li_indexs:
    eli_li_store_genes.append(li_store_genes[indexs])
    eli_li_store_sugars.append(li_store_sugars[indexs])
    eli_li_pair_processing.append(li_pair_processing[indexs])
    eli_li_pair_processing_without_one.append(li_pair_processing_without_one[indexs])
    eli_li_ser_pair_processing.append(li_ser_pair_processing[indexs])


# In[10]:


eli_total_pairgenes_df_remove_fre = pd.DataFrame({
    'Pair of genes': eli_li_store_genes, 
    'Pair of sugars': eli_li_store_sugars,
    'Num. of occurrence where pair of sugars processing (included observation)': eli_li_pair_processing,
    'Num. of occurrence where pair of sugars processing (without observation)': eli_li_pair_processing_without_one,
    'Serotype(s) where processing pair of genes and sugars': eli_li_ser_pair_processing
})

eli_total_pairgenes_df_remove_fre


# In[11]:


eli_total_pairgenes_df_remove_fre.sort_values(by=['Num. of occurrence where pair of sugars processing (without observation)'], 
                                                ascending = False)


# Let's ignore pairs of genes and corresponding pairs of sugars which appear less than four times.
# 
# Cuz which shown one to three times almost shown in each serogroups, so we eliminate them. 

# In[12]:


condition = eli_total_pairgenes_df_remove_fre['Num. of occurrence where pair of sugars processing (included observation)'] > 3

eli_total_pairgenes_df = eli_total_pairgenes_df_remove_fre.sort_values(by=['Num. of occurrence where pair of sugars processing (without observation)'], 
                                                ascending = False)[condition]


# In[13]:


eli_total_pairgenes_df 


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




