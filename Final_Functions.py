#!/usr/bin/env python
# coding: utf-8

# # Functions for final year project

# In[1]:


# from Serotype_Data import *

# from Bio.Seq import Seq
# from Bio.pairwise2 import format_alignment
# from Bio import pairwise2

import numpy as np
import tabulate


# In[2]:


# processing functions
def remove_brown_gene(gene_list):
    new_list = gene_list[:]
    for gene in gene_list:
        if '-' in gene:
            new_list.remove(gene)
    return new_list

def simplify_genes(gene_list):
    new_list = gene_list[:]
    housekeeping = ['wzg', 'wzh', 'wzd', 'wze'] #Sequence start with ['wzg', 'wzh', 'wzd', 'wze']
    ignore = ['tnp'] #Ignore gene 'tnp'
    finisher = ['aliA'] #Sequence finish with 'aliA'
    polymerase = ['wzy']
    flippase = ['wzx']
    
    to_remove = housekeeping + ignore + finisher + polymerase + flippase #['wzg', 'wzh', 'wzd', 'wze', 'tnp', 'aliA', 'wzy', 'wzx']
    
    #new_list = remove_brown_gene(new_list) #Remove brown genes
    
    for gene in to_remove:
        while gene in new_list: #Remove all of genes belong to to_remove
            new_list.remove(gene)
    
    return new_list


# ##### Pairwise sequence alignment
# https://readiab.org/pairwise-alignment.html
# 
# Fundamental problems is determining how similar a pair of biological sequences are.
# 
# We’ll explore why determining biological sequence similarity is harder than it might initially seem, and learn about pairwise sequence alignment, the standard approach for determining sequence similarity.

# In[3]:


def show_compare(h_sequence, v_sequence, data, hide_zeros=False, nonzero_val=None):
    rows = []
    col_headers = [c for c in h_sequence]
    row_headers = [c for c in v_sequence]
    pad_headers = data.shape == (len(row_headers) + 1, len(col_headers) + 1)
    if pad_headers:
        row_headers = [" "] + row_headers
        col_headers = [" "] + col_headers
    for h, d in zip(row_headers, data):
        current_row = [h]
        for e in d:
            if e == 0:
                if hide_zeros:
                    current_row.append('')
                else:
                    current_row.append(0)
            else:
                if nonzero_val is not None:
                    current_row.append(nonzero_val)
                else:
                    current_row.append(e)
        rows.append(current_row)
    return tabulate.tabulate(rows, headers=col_headers, tablefmt='html')


# In[ ]:





# In[ ]:




