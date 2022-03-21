#!/usr/bin/env python
# coding: utf-8

# # Dataset for subgroup discovery

# In[22]:


from Serotype_Data import * #Dataset of Serotypes
from Serotype_Functions import *

#Biopython ver.1.76
from Bio.Seq import Seq #Represent biological sequences with alphabets
from Bio.pairwise2 import format_alignment #Functions to get global and local alignments between two sequences
from Bio import pairwise2

from IPython.display import Image

import numpy as np
import pandas as pd
import tabulate


# In[23]:


import json
import os


# In[24]:


path = os.path.relpath('./serotypes.json')
f = open(path, 'r')
serotypes = json.load(f)


# In[25]:


total_sugars = []

for key in serotypes:
    for sugar in serotypes[key]['sugars']:
        if sugar in total_sugars:
            continue # to avoid duplicate
        else:
            total_sugars.append(sugar)


# In[26]:


fre_total_sugars = dict()

# Build basic ditionary for frequency of all sugars in dataset
for sugar in total_sugars:
    fre_total_sugars[sugar] = 0

for key in serotypes:
    for sugar in serotypes[key]['sugars']:
        fre_total_sugars[sugar] += 1


# In[27]:


li_total_sugars_keys = list(fre_total_sugars.keys())
li_total_sugars_vals = list(fre_total_sugars.values())


# In[28]:


total_sugars_fre_df = pd.DataFrame({
    "Sugar(s)":li_total_sugars_keys,
    "Num. of occurrence of sugar": li_total_sugars_vals
})
total_sugars_fre_df


# We ignore the sugars which is shown once or twice in all dataset.

# In[29]:


eli_li_total_sugars_keys = []
eli_li_total_sugars_vals = []
temp_indexs = []

for i in range(len(li_total_sugars_vals)):
    if li_total_sugars_vals[i]==1 or li_total_sugars_vals[i]==2:
        continue
    else: # more than 3
        temp_indexs.append(i)

for index in temp_indexs:
    eli_li_total_sugars_keys.append(li_total_sugars_keys[index])
    eli_li_total_sugars_vals.append(li_total_sugars_vals[index])


# In[30]:


eli_total_sugars_fre_df = pd.DataFrame({
    "Sugar(s)": eli_li_total_sugars_keys,
    "Num. of occurrence of sugar": eli_li_total_sugars_vals
})
eli_total_sugars_fre_df.sort_values(by=['Num. of occurrence of sugar'], ascending = False, inplace = True) #Re-store inplace=True


# In[31]:


eli_total_sugars_fre_df


# In[32]:


dict_from_eli_total_sugars_fre_df = eli_total_sugars_fre_df.to_dict('list')

frequent_sugar = [] 

for sugar in dict_from_eli_total_sugars_fre_df['Sugar(s)']:
    frequent_sugar.append(sugar)


# In[33]:


li_sugars = dict()

for i in frequent_sugar:
    li_sugars[i] = [0]*9
    li_sugars[i][frequent_sugar.index(i)] = 1


# In[34]:


li_sugars_df = pd.DataFrame.from_dict(li_sugars)
li_sugars_df


# In[35]:


from PairXPair_dataset import *

eli_total_pairgenes_df = eli_total_pairgenes_df.reset_index()
eli_total_pairgenes_df = eli_total_pairgenes_df.drop(['index'], axis = 1)
eli_total_pairgenes_df


# In[36]:


li_pairgenes = dict()

for pairgenes in eli_total_pairgenes_df['Pair of genes']: 
    li_pairgenes[pairgenes] = [0]*9


# In[37]:


temp_sugar_loc = dict()
sugar_loc = 0

for pairgenes in eli_total_pairgenes_df['Pair of genes']: 
    temp_sugar_loc[pairgenes] = []
    
for i in range(len(eli_total_pairgenes_df['Pair of genes'])):
    pairgenes = []
    pairgenes = eli_total_pairgenes_df['Pair of genes'][i]

    for j in eli_total_pairgenes_df['Pair of sugars'][i]: #j=name of sugar
        sugar_loc = frequent_sugar.index(j)
        temp_sugar_loc[pairgenes].append(sugar_loc)


# In[38]:


for pairgenes in li_pairgenes:
    for j in temp_sugar_loc[pairgenes]:
        if li_pairgenes[pairgenes][j] == 1:
            continue
        else:
            li_pairgenes[pairgenes][j] = 1 


# In[39]:


li_pairgenes


# In[40]:


li_rmlANrmlC = list(li_pairgenes.get(('rmlA', 'rmlC')))
li_rmlCNrmlB = list(li_pairgenes.get(('rmlC', 'rmlB')))
li_rmlBNrmlD = list(li_pairgenes.get(('rmlB', 'rmlD')))
li_wchANwchF = list(li_pairgenes.get(('wchA', 'wchF')))
li_wchJNwchK = list(li_pairgenes.get(('wchJ', 'wchK')))
li_wchANwchJ = list(li_pairgenes.get(('wchA', 'wchJ')))
li_rmlDNglf_b = list(li_pairgenes.get(('rmlD', 'glf-')))
li_wchANwchO = list(li_pairgenes.get(('wchA', 'wchO')))
li_wchONwchP = list(li_pairgenes.get(('wchO', 'wchP')))
li_wchANwciB = list(li_pairgenes.get(('wchA', 'wciB')))
li_wchPNwchQ = list(li_pairgenes.get(('wchP', 'wchQ')))
li_wchMNwchN = list(li_pairgenes.get(('wchM', 'wchN')))
li_wcySNwcrN = list(li_pairgenes.get(('wcyS', 'wcrN')))
li_wchKNwchL = list(li_pairgenes.get(('wchK', 'wchL')))
li_wchLNwchM = list(li_pairgenes.get(('wchL', 'wchM')))


# In[41]:


sub_df = li_sugars_df
sub_df["(rmlA, rmlC)"] = li_rmlANrmlC
sub_df["(rmlC, rmlB)"] = li_rmlCNrmlB
sub_df["(rmlB, rmlD)"] = li_rmlBNrmlD
sub_df["(wchA, wchF)"] = li_wchANwchF
sub_df["(wchJ, wchK)"] = li_wchJNwchK
sub_df["(wchA, wchJ)"] = li_wchANwchJ 
sub_df["(rmlD, glf-)"] = li_rmlDNglf_b
sub_df["(wchA, wchO)"] = li_wchANwchO
sub_df["(wchO, wchP)"] = li_wchONwchP
sub_df["(wchA, wciB)"] = li_wchANwciB
sub_df["(wchP, wchQ)"] = li_wchPNwchQ
sub_df["(wchM, wchN)"] = li_wchMNwchN
sub_df["(wcyS, wcrN)"] = li_wcySNwcrN
sub_df["(wchK, wchL)"] = li_wchKNwchL
sub_df["(wchL, wchM)"] = li_wchLNwchM


# In[42]:


sub_df

