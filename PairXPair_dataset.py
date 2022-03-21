#!/usr/bin/env python
# coding: utf-8

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

import json
import os


path = os.path.relpath('./serotypes.json')
f = open(path, 'r')
serotypes = json.load(f)



simply_genes = [] 
lin_sugars = [] 

wholePairs_data = dict() 

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
    
    genes_sugars_pairs = dict()
    for x in range(len(simply_genes)-1):
        genes_sugars_pairs[(simply_genes[x], simply_genes[x+1])] = sugars_pair
        
    wholePairs_data[key] = genes_sugars_pairs


wholePairs_data



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


pairgene_cols = []
pairgene_cols = list(wholePairs.keys())

pairsugar_rows = []

for pairgene in wholePairs:
    for pairsugar in wholePairs[pairgene]:
        if tuple(pairsugar) in pairsugar_rows:
            continue
        else:
            pairsugar_rows.append(tuple(pairsugar))       


num_rows = len(pairsugar_rows)
num_cols = len(pairgene_cols)
data = np.zeros(shape=(num_rows, num_cols), dtype=np.int32)

show_compare(pairgene_cols, pairsugar_rows, data)



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
pairs_fre = dict()

def rules_frequency(pairgene):
    if pairgene in pairs_fre:
        pairs_fre[pairgene] += 1
    else:
        pairs_fre[pairgene] = 1

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



pairgene_cols = []
pairgene_cols = list(wholePairs.keys())

pairsugar_rows = []

for pairgene in wholePairs:
    for pairsugar in wholePairs[pairgene]:
        if tuple(pairsugar) in pairsugar_rows:
            continue
        else:
            pairsugar_rows.append(tuple(pairsugar))    



num_rows = len(pairsugar_rows)
num_cols = len(pairgene_cols)
data = np.zeros(shape=(num_rows, num_cols), dtype=np.int32)


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


for n in range(num_rows):
    for m in range(num_cols):
        if data[n][m] == 0:
            continue
        else:
            data[n][m] += -1
            
show_compare(pairgene_cols, pairsugar_rows, data)


li_pairs_genes = list(wholePairs.keys()) #Pairs of genes


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
    
    p_genes = [] # Store genes as pair for test
    for i in range(len(simply_genes)-1):
        p_genes.append((simply_genes[i], simply_genes[i+1]))

    
    for j in range(len(p_genes)):
        for pairgenes in li_pairs_genes:
            if not pairgenes == p_genes[j]:
                continue
            else: #pairgene == p_genes[k]
                fre_pairgenes[p_genes[j]] += 1
               
            
            
li_total_pair_genes_vals = list(fre_pairgenes.values())



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


pair_occur = [0]*1194

for genes in fre_pairgenes:
    i = li_store_genes.index(genes)
    
    pair_occur[i] = fre_pairgenes[genes]



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



li_pair_occur = list(pair_occur)
li_pair_processing = list(pair_processing)    


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
                       
                        if ser_pair_processing[p]==['']:
                            ser_pair_processing[p]=key
                            print(ser_pair_processing[p])
                        else:
                            toy_list=[]
                            toy_list.append(key)                            
                            ser_pair_processing[p] = list(ser_pair_processing[p]) + toy_list


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


total_pairgenes_df_remove_fre = pd.DataFrame({
    'Pair of genes': li_store_genes, 
    'Pair of sugars': li_store_sugars,
    #'Num. of occurrence of pair of gene (Written in front of each pair of genes)': li_pair_occur,
    'Num. of occurrence where pair of sugars processing (included observation)': li_pair_processing,
    'Num. of occurrence where pair of sugars processing (without observation)': li_pair_processing_without_one,
    'Serotype(s) where processing pair of genes and sugars': li_ser_pair_processing
})




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



eli_total_pairgenes_df_remove_fre = pd.DataFrame({
    'Pair of genes': eli_li_store_genes, 
    'Pair of sugars': eli_li_store_sugars,
    'Num. of occurrence where pair of sugars processing (included observation)': eli_li_pair_processing,
    'Num. of occurrence where pair of sugars processing (without observation)': eli_li_pair_processing_without_one,
    'Serotype(s) where processing pair of genes and sugars': eli_li_ser_pair_processing
})

eli_total_pairgenes_df_remove_fre




eli_total_pairgenes_df_remove_fre.sort_values(by=['Num. of occurrence where pair of sugars processing (without observation)'], 
                                                ascending = False)


condition = eli_total_pairgenes_df_remove_fre['Num. of occurrence where pair of sugars processing (included observation)'] > 4

eli_total_pairgenes_df = eli_total_pairgenes_df_remove_fre.sort_values(by=['Num. of occurrence where pair of sugars processing (without observation)'], 
                                                ascending = False)[condition]



eli_total_pairgenes_df 





