#!/usr/bin/env python
# coding: utf-8

# # Dataset for Subgroup_43x9

# In[1]:


from Final_Functions import *
from Final_PairXPair import *

import numpy as np
import pandas as pd
import tabulate


# In[2]:


import json
import os


# In[3]:


path = os.path.relpath('./serotypes.json')
f = open(path, 'r')
serotypes = json.load(f)


# In[4]:


total_sugars = []

for key in serotypes:
    for sugar in serotypes[key]['sugars']:
        if sugar in total_sugars:
            continue # to avoid duplicate
        else:
            total_sugars.append(sugar)


# In[5]:


fre_total_sugars = dict()

# Build basic ditionary for frequency of all sugars in dataset
for sugar in total_sugars:
    fre_total_sugars[sugar] = 0

for key in serotypes:
    for sugar in serotypes[key]['sugars']:
        fre_total_sugars[sugar] += 1


# In[6]:


li_total_sugars_keys = list(fre_total_sugars.keys())
li_total_sugars_vals = list(fre_total_sugars.values())


# In[7]:


total_sugars_fre_df = pd.DataFrame({
    "Sugar(s)":li_total_sugars_keys,
    "Num. of occurrence of sugar": li_total_sugars_vals
})


# We ignore the sugars which is shown once or twice in all dataset.

# In[8]:


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


# In[9]:


eli_total_sugars_fre_df = pd.DataFrame({
    "Sugar(s)": eli_li_total_sugars_keys,
    "Num. of occurrence of sugar": eli_li_total_sugars_vals
})
eli_total_sugars_fre_df.sort_values(by=['Num. of occurrence of sugar'], ascending = False, inplace = True) #Re-store inplace=True


# In[10]:


eli_total_sugars_fre_df = eli_total_sugars_fre_df.reset_index()
eli_total_sugars_fre_df = eli_total_sugars_fre_df.drop(['index'], axis = 1)

eli_total_sugars_fre_df


# In[11]:


f_sugar = []

for frequentSugar in eli_total_sugars_fre_df['Sugar(s)']:
    #print(frequentSugar)
    f_sugar.append(frequentSugar)


# In[12]:


f_sugar


# Let's ignore pairs of genes and corresponding pairs of sugars which appear less than four times.
# 
# Cuz which shown one to three times almost shown in each serogroups, so we eliminate them. 

# In[13]:


condition = eli_total_pairgenes_df_remove_fre['Num. of occurrence where pair of sugars processing (included observation)'] > 3

eli_total_pairgenes_df = eli_total_pairgenes_df_remove_fre.sort_values(by=['Num. of occurrence where pair of sugars processing (without observation)'], 
                                                ascending = False)[condition]


# In[14]:


eli_total_pairgenes_df = eli_total_pairgenes_df.reset_index()
eli_total_pairgenes_df = eli_total_pairgenes_df.drop(['index'], axis = 1)

eli_total_pairgenes_df


# ------------------------------------------------------------------------------

# In[15]:


total_genes = []
sets = set()

for key in serotypes:
    # Simplify the sequence of genes
    simply_genes = simplify_genes(serotypes[key]['genes']) #list
    for gene in simply_genes:
        sets.add(gene)


# In[16]:


total_pairgenes = []

for key in serotypes: 
    # ignore serotypes which have no structure of sugar
    if len(serotypes[key]['sugars']) == 0:
        continue

    # with serotypes which have own sugar structure
    # simplify the gene sequence
    simply_genes = simplify_genes(serotypes[key]['genes'])
    # copy the [key]['sugars'] as local 
    lin_sugars = serotypes[key]['sugars']
    
    p_genes = [] # Store genes as pair for test
    for i in range(len(simply_genes)-1):
        string = simply_genes[i] + '-' + simply_genes[i+1]
        p_genes.append(string)
        
    for pairgene in p_genes:
        if pairgene in total_pairgenes:
            continue
        else:
            total_pairgenes.append(pairgene)
            


# In[17]:


li_pairgenes = dict()
#sugar_loc = []

for pairgenes in eli_total_pairgenes_df['Pair of genes']: 
    li_pairgenes[pairgenes] = [0]*9


# In[18]:


pairgenes = []
for elem in li_pairgenes.keys():
    pairgenes.append(elem[0]+'-'+elem[1]) 


# In[19]:


data_dict = {}
outcome = {'outcome':[]}
for elem in pairgenes:
    data_dict[elem] = [] #[key]each genes [value]


# In[20]:


outcome_dict = {}
value = 0
for elem in f_sugar:
    outcome_dict[elem] = value
    value += 1


# In[21]:


outcome_dict


# In[22]:


for key in serotypes:
    simply_genes = simplify_genes(serotypes[key]['genes'])
    p_genes = [] # Store genes as pair for test
    for i in range(len(simply_genes)-1):
        p_genes.append(simply_genes[i] +'-' +simply_genes[i+1])
    sugars = serotypes[key]['sugars']
    for sugar in sugars:
        # Only think about frequent sugars
        if sugar in outcome_dict.keys():            
            for elem in list(data_dict.keys()): #elem: each gene
                if elem in p_genes:
                    li = data_dict[elem]
                    li.append(1)
                    data_dict[elem] = li
                else:
                    li = data_dict[elem]
                    li.append(0)
                    data_dict[elem] = li
            li = outcome['outcome']
            li.append(sugar)
            outcome['outcome'] = li


# In[23]:


df = pd.DataFrame().from_dict(data_dict)


# In[24]:


df


# In[25]:


df['outcome'] = outcome['outcome']


# In[26]:


df


# In[ ]:




