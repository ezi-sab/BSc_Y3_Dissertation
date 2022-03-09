#!/usr/bin/env python
# coding: utf-8

# # Functions for finding the rules
# 

# In[15]:


from Serotype_Data import *

from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from Bio import pairwise2

import numpy as np
import tabulate


# In[16]:


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


# In[18]:


def compare_genes(list):
    alignments = pairwise2.align.globalms(simplify_genes(list[0]['genes']), simplify_genes(list[1]['genes']), 2, -1, -5, -.1, gap_char=["-"], one_alignment_only=True, penalize_end_gaps=True) # identical character, mismatched character, gap open, gap extend

    for a in alignments: 
        print(format_alignment(*a, full_sequences=True))

def compare_sugars(list):
    try: 
        list[0]['sugars'] and list[1]['sugars']
        alignments = pairwise2.align.globalms(list[0]['sugars'], list[1]['sugars'], 2, -1, -2, -.1, gap_char=["-"], penalize_end_gaps=True)

        for a in alignments: 
            print(format_alignment(*a, full_sequences=True))
    except KeyError:
        print("** No sugar structure for one or both serotypes **")

def compare_bonds(list):
    alignments = pairwise2.align.globalms(list[0]['bonds'], list[1]['bonds'], 2, -1, -2, -.1, gap_char=["-"], penalize_end_gaps=True)

    for a in alignments: 
        print(format_alignment(*a, full_sequences=True))

def compare_side_branches(list):
    alignments = pairwise2.align.globalms(list[0]['side branches'], list[1]['side branches'], 2, -1, -2, -.1, gap_char=["-"], one_alignment_only=True, penalize_end_gaps=True)

    for a in alignments: 
        print(format_alignment(*a, full_sequences=True))

def compare_lists(list):
    alignments = pairwise2.align.globalms(list[0], list[1], 2, -1, -2, -.1, gap_char=["-"], penalize_end_gaps=True, one_alignment_only=True)

    for a in alignments: 
        print(format_alignment(*a, full_sequences=True))


# In[ ]:


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


# In[19]:




