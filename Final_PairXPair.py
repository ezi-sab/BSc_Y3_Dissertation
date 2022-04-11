#!/usr/bin/env python
# coding: utf-8

# # Pair of genes x Pair of sugars

# In[1]:


from Final_Functions import *

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


# ##### Build the rules in pairs
# Previously, we build the rules one to one. 
# 
# Now, we will build the rules within the pair of genes and corresponding sugars. 

# In[4]:


simply_genes = [] #Store gene sequences of serotype which is simplified
lin_sugars = [] #Store sugar structure of serotype
# genes_pair = [] #Pair of genes in sequence
# sugars_pair = []  #Pair of sugars in sequence

wholePairs_data = dict() #Store pairgene-pairs of sugars in each serotype

for key in serotypes:
    # ignore serotypes which have no structure of sugar
    if len(serotypes[key]['sugars']) == 0:
        continue
    
    # with serotypes which have own sugar structure
    # simplify the sequence of genes
    simply_genes = simplify_genes(serotypes[key]['genes'])
    # copy the [key]['sugars'] as local 
    lin_sugars = serotypes[key]['sugars']
    
    sugars_pair = []
    for i in range(len(lin_sugars)-1): 
        li_sugars_pair = []
        li_sugars_pair.append(lin_sugars[i])
        li_sugars_pair.append(lin_sugars[i+1])
        
        sugars_pair.append(li_sugars_pair)
        
#    print(key, sugars_pair) #Check pairs of sugars
#    print(key, lin_sugars) #Check the sugars in each serotype
    
    genes_sugars_pairs = dict()
    for x in range(len(simply_genes)-1):
        genes_sugars_pairs[(simply_genes[x], simply_genes[x+1])] = sugars_pair
        
    wholePairs_data[key] = genes_sugars_pairs
      
    #print(key, genes_sugars_pairs)
      
            


# In[5]:


wholePairs_data


# In[6]:


wholePairs = dict() #Store pairs of gene - pairs of sugar

for key in wholePairs_data:
    for pairgene in wholePairs_data[key]:
        for i in range(len(wholePairs_data[key][pairgene])):
            li_draft = []
            if not pairgene in wholePairs:
                wholePairs[pairgene] = []
                wholePairs[pairgene].append(wholePairs_data[key][pairgene][i])
            else:
                if wholePairs_data[key][pairgene][i] in wholePairs[pairgene]:
                    continue                
                else:
                    wholePairs[pairgene].append(wholePairs_data[key][pairgene][i])
                
wholePairs #[key]"tuple" of pair of genes [value]List of "list" of pair of genes


# In[7]:


len(wholePairs.keys())


# In[8]:


pairgene_cols = []
pairgene_cols = list(wholePairs.keys())

pairsugar_rows = []

for pairgene in wholePairs:
    for pairsugar in wholePairs[pairgene]:
        if tuple(pairsugar) in pairsugar_rows:
            continue
        else:
            pairsugar_rows.append(tuple(pairsugar))       


# In[9]:


num_rows = len(pairsugar_rows)
num_cols = len(pairgene_cols)
data = np.zeros(shape=(num_rows, num_cols), dtype=np.int32)

show_compare(pairgene_cols, pairsugar_rows, data)


# In[10]:


for pairgene in wholePairs:
    if not pairgene in pairgene_cols: # pairgene != pairgene_cols[i]
        continue
    else: # pairgene == pairgene_cols[i]
        store_colnum = 0
        store_colnum = pairgene_cols.index(pairgene)
        #print(pairgene, store_colnum)
        
        for pairsugar in wholePairs[pairgene]:
            if not tuple(pairsugar) in pairsugar_rows:
                continue
            else:
                store_rownum = 0
                store_rownum = pairsugar_rows.index(tuple(pairsugar))
                #print(pairsugar, store_rownum)
                data[store_rownum,store_colnum] = 1
                        
show_compare(pairgene_cols, pairsugar_rows, data)


# ##### Check frequency
# Check the frequency of wholePairs

# In[11]:


pairs_fre = dict()

def rules_frequency(pairgene):
    if pairgene in pairs_fre:
        pairs_fre[pairgene] += 1
    else:
        pairs_fre[pairgene] = 1


# In[12]:


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
        p_genes.append((simply_genes[i], simply_genes[i+1]))
    
    p_sugars = [] # Store sugars as pair for test
    for j in range(len(lin_sugars)-1):
         p_sugars.append((lin_sugars[j], lin_sugars[j+1]))
    
    for k in range(len(p_genes)):
        for pairgene in wholePairs:
            
            if not (pairgene == p_genes[k]):
                continue
            else: # pairgene == p_genes[k]
                rules_frequency(pairgene)
                               
pairs_fre


# Clear the data of dataFrame to check frequency

# In[13]:


pairgene_cols = []
pairgene_cols = list(wholePairs.keys())

pairsugar_rows = []

for pairgene in wholePairs:
    for pairsugar in wholePairs[pairgene]:
        if tuple(pairsugar) in pairsugar_rows:
            continue
        else:
            pairsugar_rows.append(tuple(pairsugar))    


# In[14]:


num_rows = len(pairsugar_rows)
num_cols = len(pairgene_cols)
data = np.zeros(shape=(num_rows, num_cols), dtype=np.int32)


# In[15]:


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
        p_genes.append((simply_genes[i], simply_genes[i+1]))
    
    # print(p_genes) #List of tuple of pair of genes
    
    p_sugars = [] # Store sugars as pair for test
    for j in range(len(lin_sugars)-1):
         p_sugars.append((lin_sugars[j], lin_sugars[j+1]))
            
    # print(p_sugars) #List of tuples of pair of sugars
    
    for k in range(len(p_genes)):
        for pairgene in wholePairs:
            if not pairgene == p_genes[k]:
                continue
            else: #pairgene == p_genes[k]
                store_colnum = 0
                store_colnum = pairgene_cols.index(pairgene)
        
                for r in range(len(p_sugars)): 
                    #print(key, k ,r, p_sugars[r])
                    if list(p_sugars[r]) in wholePairs.get(pairgene):
                        store_rownum = 0
                        store_rownum = pairsugar_rows.index(p_sugars[r])
                        
                        data[store_rownum,store_colnum] += 1

                        
show_compare(pairgene_cols, pairsugar_rows, data)


# In[16]:


for n in range(num_rows):
    for m in range(num_cols):
        if data[n][m] == 0:
            continue
        else:
            data[n][m] += -1
            
show_compare(pairgene_cols, pairsugar_rows, data)


# In[17]:


# Set maximum num of rows & cols
#pd.set_option('display.max_rows', None)
#pd.set_option('display.max_columns', None)


# In[18]:


li_pairs_genes = list(wholePairs.keys()) #Pairs of genes


# In[19]:


#Check how many times the pairs of genes appear in dataset
fre_pairgenes = dict()

for pairgenes in li_pairs_genes:
    fre_pairgenes[pairgenes] = 0

for key in serotypes:
    # ignore serotypes which have no structure of sugar
    if len(serotypes[key]['sugars']) == 0:
        continue
        
    # with serotypes which have own sugar structure
    # simplify the gene sequence
    simply_genes = simplify_genes(serotypes[key]['genes'])
#     # copy the [key]['sugars'] as local 
#     lin_sugars = serotypes[key]['sugars']
    
    p_genes = [] # Store genes as pair for test
    for i in range(len(simply_genes)-1):
        p_genes.append((simply_genes[i], simply_genes[i+1]))
    
    # print(p_genes) #List of tuple of pair of genes
    
#     p_sugars = [] # Store sugars as pair for test
#     for j in range(len(lin_sugars)-1):
#          p_sugars.append((lin_sugars[j], lin_sugars[j+1]))
            
    # print(p_sugars) #List of tuples of pair of sugars
    
    for j in range(len(p_genes)):
        for pairgenes in li_pairs_genes:
            if not pairgenes == p_genes[j]:
                continue
            else: #pairgene == p_genes[k]
                fre_pairgenes[p_genes[j]] += 1
                


# In[20]:


li_total_pair_genes_vals = list(fre_pairgenes.values())


# In[21]:


store_genes = []
store_sugars = []
li_store_genes = []
li_store_sugars = []

for pairgenes in li_pairs_genes:     
    store_genes = pairgenes
    
    for pairsugars in wholePairs[pairgenes]:                
        store_sugars = pairsugars
        #print(store_sugars)
        
        li_store_genes.append(store_genes)
        li_store_sugars.append(store_sugars)


# In[22]:


pair_occur = [0]*1194

for genes in fre_pairgenes:
    i = li_store_genes.index(genes)
    
    pair_occur[i] = fre_pairgenes[genes]


# In[23]:


pair_processing = [0]*1194
pair_not_processing = [0]*1194

for key in serotypes: 
    # ignore serotypes which have no structure of sugar
    if len(serotypes[key]['sugars']) == 0:
        continue

    # with serotypes which have own sugar structure
    simply_genes = simplify_genes(serotypes[key]['genes']) # simplify the gene sequence
    lin_sugars = serotypes[key]['sugars'] # copy the [key]['sugars'] as local 
    
    p_genes = [] # Store genes as pair for test
    for i in range(len(simply_genes)-1):
        p_genes.append((simply_genes[i], simply_genes[i+1]))
    
    p_sugars = [] # Store sugars as pair for test
    for j in range(len(lin_sugars)-1):
         p_sugars.append((lin_sugars[j], lin_sugars[j+1]))
            
    for k in range(len(p_genes)):
        for d in range(len(li_store_genes)):
            if not p_genes[k] == li_store_genes[d]:
                continue
            
            else: # p_genes[k] == li_store_genes[d]
                for r in range(len(p_sugars)):
                    if not list(p_sugars[r]) == li_store_sugars[d]:
                        continue
                        #print("!=", p_genes[k], p_sugars[r], li_store_sugars[li_store_genes.index(x)])
                    else: # p_sugars[r] == li_store_sugars[x]
                        pair_processing[d] += 1
                        #print("==", p_genes[k], p_sugars[r], li_store_sugars[li_store_genes.index(x)])                       


# In[24]:


li_pair_occur = list(pair_occur)
li_pair_processing = list(pair_processing)    


# In[25]:


pair_processing_without_one = []
pair_processing_without_one = pair_processing

li_pair_processing_without_one = []


for i in pair_processing_without_one:
    if i == 1:
        i = 0
        li_pair_processing_without_one.append(i)
    else:
        i = i - 1
        li_pair_processing_without_one.append(i)


# In[26]:


ser_pair_processing = ['']*1194
ser_pair_not_processing = []
    
for key in serotypes: 
    # ignore serotypes which have no structure of sugar
    if len(serotypes[key]['sugars']) == 0:
        continue

    # with serotypes which have own sugar structure
    simply_genes = simplify_genes(serotypes[key]['genes']) # simplify the gene sequence
    lin_sugars = serotypes[key]['sugars'] # copy the [key]['sugars'] as local 
    
    p_genes = [] # Store genes as pair for test
    for i in range(len(simply_genes)-1):
        p_genes.append((simply_genes[i], simply_genes[i+1]))
    
    p_sugars = [] # Store sugars as pair for test
    for j in range(len(lin_sugars)-1):
         p_sugars.append((lin_sugars[j], lin_sugars[j+1]))
            
    for k in range(len(p_genes)):
        #Check li_store_genes include p_genes[k] 
        if p_genes[k] in li_store_genes:
            #Store all index where p_genes[k] in li_store_genes
            pairgenes_index = [p for p, x in enumerate(li_store_genes) if x==p_genes[k]]
            #print(p_genes[k], pairgenes_index)
            for r in range(len(p_sugars)):
                for p in pairgenes_index:
                    #print(r,p)
                    if list(p_sugars[r])==li_store_sugars[p]: 
#                         print(p_sugars[r],"==", li_store_sugars[p])
#                     else:
#                        print(p_sugars[r],"!=", li_store_sugars[p])                        
                        if ser_pair_processing[p]==['']:
                            ser_pair_processing[p]=key
                            print(ser_pair_processing[p])
                        else:
                            toy_list=[]
                            toy_list.append(key)                            
                            ser_pair_processing[p] = list(ser_pair_processing[p]) + toy_list


# In[37]:


modify_ser_pair_processing = ['']*1194
#restore_ser_pair_processing = 

for i in range(len(ser_pair_processing)):
    sers = ser_pair_processing[i]
    #print(sers)
    for ser in sers:
        #print((ser))
        if not ser in modify_ser_pair_processing[i]:
            toy_list = []
            toy_list.append(ser)
            modify_ser_pair_processing[i] = list(modify_ser_pair_processing[i]) + toy_list

#modify_ser_pair_processing
li_ser_pair_processing = modify_ser_pair_processing


# In[38]:


total_pairgenes_df_with_ser = pd.DataFrame({
    'Pair of genes': li_store_genes, 
    'Pair of sugars': li_store_sugars,
    'Num. of occurrence of pair of gene (Written in front of each pair of genes)': li_pair_occur,
    'Num. of occurrence where pair of sugars processing (included observation)': li_pair_processing,
    'Num. of occurrence where pair of sugars processing (without observation)': li_pair_processing_without_one,
    #"Num. of occurrence where the pairs of sugars are absent":
    'Serotype(s) where processing pair of genes and sugars': li_ser_pair_processing
})

total_pairgenes_df_with_ser


# In[40]:


total_pairgenes_df_with_ser = total_pairgenes_df_with_ser.sort_values(by=['Num. of occurrence where pair of sugars processing (without observation)'], 
                                                                      ascending = False)
total_pairgenes_df_with_ser = total_pairgenes_df_with_ser.reset_index()
total_pairgenes_df_with_ser = total_pairgenes_df_with_ser.drop(['index'], axis = 1)
total_pairgenes_df_with_ser


# Eliminate the elements which only appear in observation

# In[45]:


total_pairgenes_df_remove_fre = pd.DataFrame({
   'Pair of genes': li_store_genes, 
   'Pair of sugars': li_store_sugars,
   'Num. of occurrence where pair of sugars processing (included observation)': li_pair_processing,
   'Num. of occurrence where pair of sugars processing (without observation)': li_pair_processing_without_one,
   'Serotype(s) where processing pair of genes and sugars': li_ser_pair_processing
})


# In[46]:


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


# In[49]:


eli_total_pairgenes_df_remove_fre = pd.DataFrame({
    'Pair of genes': eli_li_store_genes, 
    'Pair of sugars': eli_li_store_sugars,
    'Num. of occurrence where pair of sugars processing (included observation)': eli_li_pair_processing,
    'Num. of occurrence where pair of sugars processing (without observation)': eli_li_pair_processing_without_one,
    'Serotype(s) where processing pair of genes and sugars': eli_li_ser_pair_processing
})

eli_total_pairgenes_df_remove_fre


# In[50]:


eli_total_pairgenes_df_remove_fre = eli_total_pairgenes_df_remove_fre.sort_values(by=['Num. of occurrence where pair of sugars processing (without observation)'], 
                                                                      ascending = False)
eli_total_pairgenes_df_remove_fre = eli_total_pairgenes_df_remove_fre.reset_index()
eli_total_pairgenes_df_remove_fre = eli_total_pairgenes_df_remove_fre.drop(['index'], axis = 1)

eli_total_pairgenes_df_remove_fre


# In[ ]:




