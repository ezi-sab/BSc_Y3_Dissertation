#!/usr/bin/env python
# coding: utf-8

# # Serotype data

# In[ ]:





# In[4]:


#Encode
#(nameOfSerotype) = {
#    'genes': [-],
    # we read the sugar sequences from right to left
#    'sugars': [n],
#    'side branches': [n], //ignore bond of side branches
#    'modifications': [n],
#    'bonds': [n-1]
#}


# In[5]:





# In[6]:


# serotype 1
ser_1 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchB', 'wchC', 'wchD', 'wzy', 'wzx', 'gla', 'ugd', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['2-acetamido-4-amino-2,4,6-trideoxygalactose', 'galacturonic acid', 'galacturonic acid'],
    'modifications': ['', '', ''],
    'bonds': ['a1-3', 'a1-3']
}

# serotype 2
ser_2 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wchG', 'wchH', 'wzy', 'wchI', 'wzx', 'ugd', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'rhamnose', 'rhamnose', 'glucose', 'glucuronic acid'],
    'modifications': ['', '', '', '','',''],
    'bonds': ['b1-4', 'a1-3', 'a1-3', 'a1-2', 'a1-6']
}

# serotype 3
ser_3 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'tnp', 'ugd', 'wchE', 'galU', 'pgm'],
    'sugars': ['glucose', 'glucuronic acid'],
    'modifications': ['', ''],
    'bonds': ['b1-4']
}

# serotype 4
ser_4 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'HG261', 'wciI', 'wciJ', 'wciK', 'wciL', 'wzy', 'wciM', 'wzx', 'mnaA', 'fnlA', 'fnlB', 'fnlC', 'tnp', 'aliA'],
    'sugars': ['N-acetylgalactosamine', 'N-acetylfucosamine', 'N-acetylmannosamine', 'galactose', 'pyruvate'],
    'modifications': ['', '', '', '',''],
    'bonds': ['a1-3', 'b1-3', 'a1-3', '2,3(S)']
}

# serotype 5
ser_5 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wciI', 'wciJ', 'wzy', 'wzx', 'whaC', 'whaD', 'whaE', 'fnlA', 'ugd', 'fnlB', 'fnlC', 'tnp', 'tnp', 'tnp', 'aliA'],
    'sugars': ['4-keto-N-acetyl-D-quinovosamine','N-acetylfucosamine','glucose'],
    'side branches': ['', 'P-glucuronic acid-PrepNAc', ''],
    'modifications': ['', '', ''],
    'bonds': ['a1-3', 'b1-4']
}

# serotype 6A

ser_6A = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciN', 'HG262', 'wciO', 'wciP', 'wzy', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'], #missing rmlB
    'sugars': ['glucose', 'galactose', 'ribitol', 'rhamnose'],
    'modifications': ['', '', '', ''],
    'bonds': ['a1-3', '5-P-2', 'a1-3']
}

# serotype 6B

ser_6B = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciN', 'HG263', 'wciO', 'wciP', 'wzy', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'galactose', 'ribitol', 'rhamnose'],
    'modifications': ['', '', '', ''],
    'bonds': ['a1-3', '5-P-2', 'a1-4']
}

# serotype 7F
ser_7F = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wcwA', 'wcwC', 'wcwD', 'HG140', 'wcwF', 'wcwG', 'wcwH', 'wzy', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'tnp', 'aliA'], #under tnp
    'sugars': ['glucose', 'rhamnose', 'galactose', 'N-acetylgalactosamine', 'rhamnose', 'N-acetylglucosamine'],
    'side branches': ['', '', 'P-galactose', '', '', ''],
    'modifications': ['', '2Ac', '', '', '', ''],
    'bonds': ['b1-4', 'a1-3', 'b1-6', 'a1-4', 'a1-2']
}

# serotype 7A
ser_7A = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wcwA', 'wcwC', 'wcwD', 'HG140', 'wcwF', 'wcwG', 'wcwH', 'wzy', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'tnp', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'galactose', 'N-acetylgalactosamine', 'rhamnose', 'N-acetylglucosamine'],
    'modifications': ['', '2Ac', '', '', '', ''],
    'bonds': ['b1-4', 'a1-3', 'b1-6', 'a1-4', 'a1-2']
}

# serotype 7B
ser_7B = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wcwI', 'wcwL', 'wcwK', 'wcxU', 'wzy', 'rsbF', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'rhamnose', 'N-acetylglucosamine', 'glucose', 'rhamnose', 'ribose-f'],
    'modifications': ['', '', '', '', '', '', ''],
    'bonds': ['b1-4', 'a1-2', 'a1-2', 'a1P-6', 'a1-3', 'b1-4']
}

# serotype 7C - No structure
ser_7C = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wcwI', 'wcwL', 'wcwK', 'wcxU', 'wzy', 'rbsF', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 8
ser_8 = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wciQ', 'wciR', 'wciS', 'wzx', 'wciT', 'wzy', 'ugd', 'HG265', 'HG266', 'aliA'],
    'sugars': ['glucose', 'glucuronic acid', 'galactose', 'glucose'],
    'modifications': ['', '', '', ''],
    'bonds': ['b1-4', 'a1-4', 'a1-4']
}

# serotype 9A
ser_9A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchO', 'wcjA', 'mnaA', 'wzy', 'wcjB', 'wzx', 'wcjC', 'wcjD', 'ugd', 'wcjE', 'aliA'],
    'sugars': ['glucose', 'N-acetylmannosamine', 'galactose', 'glucuronic acid', 'glucose'],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'a1-3', 'a1-3', 'a1-4']
}

# serotype 9L
ser_9L = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchO', 'wcjA', 'mnaA', 'wzy', 'wcjB', 'wzx', 'wcjC', 'ugd', 'tnp', 'wcjE', 'aliA'],
    'sugars': ['glucose', 'N-acetylmannosamine', 'galactose', 'glucuronic acid', 'N-acetylglucosamine'],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'a1-3', 'a1-3', 'a1-4']
}

# serotype 9N
ser_9N = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchO', 'wcjA', 'mnaA', 'wzy', 'wcjB', 'wzx', 'wcjC', 'ugd', 'tnp', 'wcjE', 'aliA'],
    'sugars': ['glucose', 'N-acetylmannosamine', 'glucose', 'glucuronic acid', 'N-acetylglucosamine'],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'a1-3', 'a1-3', 'a1-4']
}

# serotype 9V
ser_9V = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchO', 'wcjA', 'mnaA', 'wzy', 'wcjB', 'wzx', 'wcjC', 'wcjD', 'ugd', 'wcjE', 'aliA'],
    'sugars': ['glucose', 'N-acetylmannosamine', 'galactose', 'glucuronic acid', 'glucose'],
    'modifications': ['', '6Ac4Ac', '', '3Ac2Ac', '3Ac2Ac'],
    'bonds': ['b1-4', 'a1-3', 'a1-3', 'a1-4']
}

# serotype 10F
ser_10F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wcjG', 'wciB', 'wzy', 'wcrB', 'wcrC', 'wcrD', 'wciF', 'wzx', 'wciG', 'wcrH', 'aliA'],
    'sugars': ['galactose', 'galactose', 'ribitol', 'galactose', 'N-acetylgalactosamine'],
    'side branches': ['', '', '', 'F-galactose', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-3', '5-P-6', 'a1-2', 'b1-3']
}

# serotype 10A
ser_10A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wcjG', 'wciB', 'wzy', 'wcrB', 'wcrC', 'wcrD', 'wciF', 'wcrG', 'wzx', 'glf', 'aliA'],
    'sugars': ['galactose', 'galactose', 'ribitol', 'galactose', 'N-acetylgalactosamine', 'galactose'],
    'side branches': ['', '', '', '', 'P-galactose', ''],
    'modifications': ['', '', '', '', '', ''],
    'bonds': ['b1-3', '5-P-5', 'a1-2', 'b1-3', 'b1-3']
}

# serotype 10B - No structure
ser_10B = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wcjG', 'wciB', 'wzy', 'wcrB', 'wcrC', 'wcrD', 'wciF', 'wcrG', 'wzx', 'glf', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 10C - No structure
ser_10C = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wcjG', 'wciB', 'wzy', 'wcrB', 'wcrC', 'wcrD', 'wciF', 'wzx', 'wciG', 'glf', 'wcrH', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 11F
ser_11F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wcyK', 'wcwC', 'wcrL', 'wzy', 'wcwT', 'wcwU', 'wzx', 'gct', 'wcjE', 'aliA'],
    'sugars': ['glucose', 'galactose', 'galactose', 'N-acetylglucosamine', 'ribitol'],
    'modifications': ['', '', '2Ac', '3Ac0.5', ''],
    'bonds': ['b1-4', 'a1-3', 'a1-4', '1-P-4']
}

# serotype 11A
ser_11A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wcyK', 'wcwC', 'wcrL', 'wzy', 'wcwT', 'wcwU', 'wzx', 'gct', 'wcjE', 'aliA', 'aliA'],
    'sugars': ['glucose', 'galactose', 'galactose', 'glucose', 'glycerol'],
    'modifications': ['', '', '', '2/3Ac', ''],
    'bonds': ['b1-4', 'a1-3', 'a1-4', '1-P-4']
}

# serotype 11B
ser_11B = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wcyK', 'wcwR', 'wcrL', 'wzy', 'wcwT', 'wcwU', 'wzx', 'gct', 'aliA'],
    'sugars': ['glucose', 'galactose', 'galactose', 'N-acetylglucosamine', 'ribitol'],
    'modifications': ['', '', '', '3Ac0.9', ''],
    'bonds': ['b1-4', 'a1-3', 'a1-4', '1-P-4']
}

# serotype 11C
ser_11C = {
    'genes': ['tnp', 'wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wcyK', 'wcwR', 'wcrL', 'wzy', 'wcwT', 'wcwU', 'wzx', 'gct', 'aliA'],
    'sugars': ['glucose', 'galactose', 'galactose', 'N-acetylglucosamine', 'glycerol'],
    'modifications': ['', '', '', '3Ac', ''],
    'bonds': ['b1-4', 'a1-3', 'a1-4', '1-P-4']
}

# serotype 11D - No structure
ser_11D = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wcyK', 'wcwC', 'wcrL', 'wzy', 'wcwT', 'wcwU', 'wzx', 'gct', 'wcjE', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 12F
ser_12F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wciI', 'wciJ', 'wcxB', 'wzy', 'wcxD', 'wcxE', 'wcxF', 'wzx', 'mnaB', 'mnaA', 'fnlA', 'fnlB', 'fnlC', 'tnp', 'tnp', 'aliA'],
    'sugars': ['N-acetylgalactosamine', 'N-acetylfucosamine', 'N-acetylmannosaminuronic acid', 'glucose', 'glucose'],
    'side branches': ['', 'P-galactose', '', '', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['a1-3', 'b1-4', 'a1-3', 'a1-2']
}

# serotype 12A
ser_12A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'tnp', 'wciI', 'wciJ', 'wcxB', 'wzy', 'wcxD', 'wcxE', 'wcxF', 'wzx', 'mnaB', 'mnaA', 'fnlA', 'fnlB', 'fnlC', 'tnp', 'tnp', 'aliA'],
    'sugars': ['N-acetylgalactosamine', 'N-acetylfucosamine', 'N-acetylmannosaminuronic acid', 'glucose', 'glucose'],
    'side branches': ['', 'P-N-acetylgalactosamine', '', '', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['a1-3', 'b1-4', 'a1-3', 'a1-2']
}

# serotype 12B - No structure
ser_12B = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wciI', 'wciJ', 'wcxB', 'wzy', 'wcxD', 'wcxE', 'wcxF', 'wzx', 'mnaB', 'mnaA', 'fnlA', 'fnlB', 'fnlC', 'tnp', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 13
ser_13 = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'whaG', 'abp1', 'abp2', 'wciF', 'wcrD', 'wzy', 'wzx', 'wciG', 'glf', 'aliA'],
    'sugars': ['glucose', 'galactose', 'ribitol', 'N-acetylglucosamine', 'galactose'],
    'modifications': ['2/3Ac', '', '', '', ''],
    'bonds': ['b1-4', '5-P-4', 'b1-4', 'b1-4']
}

# serotype 14
ser_14 = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wzy', 'wchL', 'wchM', 'wchN', 'wzx', 'wciY', 'lrp', 'aliA'],
    'sugars': ['glucose', 'galactose', 'N-acetylglucosamine', 'galactose'],
    'modifications': ['', '', '', ''],
    'bonds': ['b1-4', 'b1-3', 'b1-4']
}

# serotype 15F
ser_15F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wzy', 'wchL', 'wchM', 'wchN', 'wzx', 'wciZ', 'wchX', 'gtp1', 'gtp2', 'gtp3', 'rmlB', 'rmlD', 'glf', 'wcjE', 'aliA'],
    'sugars': ['glucose', 'galactose', 'N-acetylglucosamine', 'galactose', 'galactose'],
    'side branches': ['', '', '', 'choline0.2', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'b1-3', 'b1-4', 'a1-2']
}

# serotype 15A
ser_15A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wzy', 'wchL', 'wchM', 'wchN', 'wzx', 'wciZ', 'wchX', 'gtp1', 'gtp2', 'gtp3', 'tnp', 'aliA'],
    'sugars': ['glucose', 'galactose', 'N-acetylglucosamine', 'galactose', 'galactose'],
    'side branches': ['', '', '', 'glycerol0.7', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'b1-3', 'b1-4', 'a1-2']
}

# serotype 15B
ser_15B = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wzy', 'wchL', 'wchM', 'wchN', 'wzx', 'wciZ', 'wchX', 'gtp1', 'gtp2', 'gtp3', 'tnp', 'aliA'],
    'sugars': ['glucose', 'galactose', 'N-acetylglucosamine', 'galactose', 'galactose'],
    'side branches': ['', '', '', 'choline0.2', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'b1-3', 'b1-4', 'a1-2']
}

# serotype 15C
ser_15C = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchJ', 'wchK', 'wzy', 'wchL', 'wchM', 'wchN', 'wzx', 'wciZ', 'wchX', 'gtp1', 'gtp2', 'gtp3', 'tnp', 'aliA'],
    'sugars': ['glucose', 'galactose', 'N-acetylglucosamine', 'galactose', 'galactose'],
    'side branches': ['', '', '', 'choline', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'b1-3', 'b1-4', 'a1-2']
}

# serotype 16F - No structure
ser_16F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wciU', 'wcxM', 'wcxN', 'HG191', 'wzy', 'wcxP', 'wzx', 'wcxQ', 'gct', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 16A - No structure
ser_16A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'gct', 'wcxR', 'wchJ', 'wchK', 'wcyK', 'wcxS', 'wzy', 'wcxT', 'wciB', 'wzx', 'wciG', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 17F
ser_17F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wcxG', 'abp1', 'abp2', 'wciP', 'wcrT', 'wcrU', 'wzy', 'wcrV', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'arabinitol', 'rhamnose', 'rhamnose', 'galactose'],
    'side branches': ['', '', '', '', 'P-galactose', ''],
    'modifications': ['', '', '', '', '2Ac', ''],
    'bonds': ['b1-4', '1-P-3', 'a1-2', 'b1-4', 'a1-3']
}

# serotype 17A
ser_17A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wciB', 'wcrP', 'wcrQ', 'wcrR', 'wcrT', 'wcrU', 'wcrV', 'wzy', 'wzx', 'ugd', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'galactose', 'glucuronic acid', ' rhamnose', 'rhamnose', 'galactose'],
    'side branches': ['', 'P-glucose', '', '', 'P-galactose', ''],
    'modifications': ['', '', '', '', '2Ac', ''],
    'bonds': ['b1-3', 'b1-3', 'a1-4', 'b1-4', 'a1-3']
}

# serotype 18F
ser_18F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wciU', 'wcxM', 'wciV', 'wciW', 'wzx', 'wzy', 'wciX', 'wciY', 'gct', 'HG94', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'glucose', 'galactose', 'glucose'],
    'side branches': ['', '', '', 'glycerol', ''],
    'modifications': ['', '2Ac', '', '', '6Ac'],
    'bonds': ['b1-4', 'a1-3', 'b1-4', 'a1-2']
}

# serotype 18A
ser_18A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wciU', 'wciV', 'wciW', 'wzx', 'wzy', 'wciY', 'gct', 'HG94', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'N-acetylglucosamine', 'galactose', 'glucose'],
    'side branches': ['', '', '', 'D-glycerol', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'a1-3', 'b1-4', 'a1-2']
}

# serotype 18B
ser_18B = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wciU', 'wciV', 'wciW', 'wzx', 'wzy', 'wciX', 'wciY', 'gct', 'HG94', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'glucose', 'galactose', 'glucose'],
    'side branches': ['', '', '', 'D-glycerol', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'a1-3', 'b1-4', 'a1-2']
}

# serotype 18C
ser_18C = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wciU', 'wciV', 'wciW', 'wzx', 'wzy', 'wciX', 'wciY', 'gct', 'HG94', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'glucose', 'galactose', 'glucose'],
    'side branches': ['', '', '', 'glycerol', ''],
    'modifications': ['', '', '', '', '6Ac0.3'],
    'bonds': ['b1-4', 'a1-3', 'b1-4', 'a1-2']
}

# serotype 19F
ser_19F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchO', 'wchP', 'wchQ', 'wzy', 'wzx', 'mnaA', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'N-acetylmannosamine', 'rhamnose'],
    'modifications': ['', '', ''],
    'bonds': ['b1-4', '1-P-4']
}

# serotype 19A
# Doesn't include aliB at the begining
ser_19A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchO', 'wchP', 'wchQ', 'wzy', 'wzx', 'mnaA', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'N-acetylmannosamine', 'rhamnose'],
    'modifications': ['', '', ''],
    'bonds': ['b1-4', '1-P-4']
}

# serotype 19B
ser_19B = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchO', 'wchP', 'wchQ', 'wchR', 'wzy', 'wchS', 'rbsF', 'wzx', 'mnaA', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'N-acetylmannosamine', 'rhamnose', 'N-acetylmannosamine', 'rhamnose', 'D-Ribf'],
    'modifications': ['', '', '', '', '', ''],
    'bonds': ['b1-4', '1-P-4', 'b1-4', 'a1-3', 'b1-4']
}

# serotype 19C
ser_19C = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchO', 'wchP', 'wchQ', 'wchR', 'wzy', 'wchS', 'rbsF', 'wzx', 'mnaA', 'wchU', 'HG264',  'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'N-acetylmannosamine', 'rhamnose', 'N-acetylmannosamine', 'rhamnose', 'D-Ribf'],
    'side branches': ['', '', '', 'P-glucose', '', ''],
    'modifications': ['', '', '', '', '', ''],
    'bonds': ['b1-4', '1-P-4', 'b1-4', 'a1-3', 'b1-4']
}

# serotype 20
ser_20 = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wciB', 'wzy', 'whaJ', 'wciL', 'wcwK', 'wciD', 'whaF', 'wzx', 'wciG', 'glf', 'wcjE', 'tnp', 'tnp', 'aliA'],
    'sugars': ['glucose', 'galactose', 'glucose', 'glucose', 'N-acerylglucosamine', 'galactose'],
    'modifications': ['', '5,6Ac2', '', '', '', ''],
    'bonds': ['b1-3', 'b1-3', 'a1-6', 'a1-P-6', 'b1-4']
}

# serotype 21 - No structure
ser_21 = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wcwA', 'wcwK', 'wcyT', 'wcyU', 'wzy', 'wzx', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'wcyO', 'tnp', 'aliA'] #,
#    'sugars': ['', '', '', ''],
#    'modifications': ['', '', '', ''],
#    'bonds': ['', '', '', '']
}

# serotype 22F
ser_22F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wcwA', 'wcwC', 'ugd', 'wcwV', 'whaB', 'wzy', 'wcwX', 'wzx', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'glucuronic acid', 'rhamnose', 'galactose'],
    'side branches': ['', 'P-glucose', '', '', ''],
    'modifications': ['', '2AC0.8', '', '', ''],
    'bonds': ['b1-4', 'b1-4', 'a1-4', 'a1-2']
}

# serotype 22A - No structure
ser_22A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wcwA', 'wcwC', 'ugd', 'wcwV', 'whaB', 'wzy', 'wcwX', 'wzx', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 23F
ser_23F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wzy', 'wchV', 'wchW', 'wzx', 'wchX', 'gtp1', 'gtp2', 'gtp3', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'galactose', 'rhamnose'],
    'side branches': ['', '', 'glycerol', ''],
    'modifications': ['', '', '', ''],
    'bonds': ['b1-4', 'b1-4', 'a1-2']
}

# serotype 23A - No structure
ser_23A = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wzy', 'wchV', 'wchW', 'wzx', 'wchX', 'gtp1', 'gtp2', 'gtp3', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 23B - No structure
ser_23B = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wzy', 'wchV', 'wchW', 'wzx', 'wchX', 'gtp1', 'gtp2', 'gtp3', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 24F - No structure
ser_24F = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wcxG', 'abp1', 'abp2', 'HG267', 'wzy', 'wcxI', 'wcxJ', 'wcxK', 'wzy', 'rbsF', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'],
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 24A - No structure
ser_24A = {
    'genes': ['tnp', 'wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wcxG', 'abp1', 'abp2', 'HG268', 'wzy', 'wcxI', 'wcxJ', 'wcxK', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 24B - No structure
ser_24B = {
    'genes': ['wzg', 'wzh', 'wzd', 'wze', 'wchA', 'wchF', 'wcxG', 'abp1', 'abp2', 'HG269', 'wzy', 'wcxI', 'wcxJ', 'wcxK', 'wzy', 'rbsF', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 25F - No structure
ser_25F = {
    'genes': ['tnp', 'glf', 'wzd', 'wze', 'tnp', 'wzg', 'wzh', 'wciI', 'wcyA', 'wzy', 'wcyB', 'wcyC', 'wcyD', 'wcyE', 'wzx', 'wcyF', 'gla', 'ugd', 'tnp', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 25A - No structure
ser_25A = {
    'genes': ['tnp', 'glf', 'wzd', 'wze', 'tnp', 'wzg', 'wzh', 'wciI', 'wcyA', 'wzy', 'wcyB', 'wcyC', 'wcyD', 'wcyE', 'wzx', 'wcyF', 'gla', 'ugd', 'tnp', 'tnp', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 27
ser_27 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'whaK', 'wzy', 'whaL', 'wzx', 'wcyS', 'wcrN', 'HG270', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'galactose', 'N-acetylglucosamine', 'pyruvate'],
    'side branches': ['', 'choline', '', '', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-4', 'a1-4', 'b1-3', '4,6(S)']
}

# serotype 28F - No structure
ser_28F = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wciU', 'wcxM', 'wcxN', 'wzy', 'wcxP', 'wzx', 'wcxQ', 'gtp1', 'gtp2', 'gtp3', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 28A - No structure
ser_28A = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wciU', 'wcxM', 'wcxN', 'wzy', 'wcxP', 'wzx', 'wcxQ', 'gtp1', 'gtp2', 'gtp3', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 29
ser_29 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wcjH', 'wciB', 'wzy', 'wcrJ', 'wcrM', 'wcrH', 'wzx', 'glf', 'aliA'],
    'sugars': ['galactose', 'galactose', 'N-acetylgalactosamine', 'ribitol', 'galactose'],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-3', 'b1-6', '5-P-4', 'b1-1']
}

# serotype 31
ser_31 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wcjG', 'wciB', 'wcrP', 'wcrR', 'wzy', 'wcrW', 'wzx', 'wcrX', 'ugd', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'wcjE', 'aliA'],
    'sugars': ['galactose', 'glucuronic acid', 'rhamnose', 'galactose', 'rhamnose'],
    'modifications': ['', '', '', '', ''],
    'bonds': ['b1-3', 'b1-4', 'b1-3', 'b1-3']
}

# serotype 32F
ser_32F = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wcyH', 'wzy', 'wcyI', 'wchQ', 'wzx', 'wcyS', 'wcrN', 'HG272', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'glucose', 'rhamnose'],
    'side branches': ['', 'choline', '', ''],
    'modifications': ['', '2Ac', '', ''],
    'bonds': ['b1-4', 'a1-4', 'a1-P-2']
}

# serotype 32A
ser_32A = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wcyH', 'wzy', 'wcyI', 'wchQ', 'wzx', 'wcyS', 'wcrN', 'HG271', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['glucose', 'rhamnose', 'glucose', 'rhamnose'],
    'side branches': ['', 'choline', '', ''],
    'modifications': ['', '2Ac', '4Ac', ''],
    'bonds': ['b1-4', 'a1-4', 'a1-P-2']
}

# serotype 33F
ser_33F = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wciC', 'wciD', 'wciE', 'wciF', 'wzy', 'wzx', 'wciG*', 'glf', 'wcjE', 'aliA'],
    'sugars': ['glucose', 'galactose', 'galactose', 'galactose', 'galactose'],
    'side branches': ['', '', 'P-galactose', '', ''],
    'modifications': ['', '2Ac0.4', '', '', ''],
    'bonds': ['b1-3', 'a1-3', 'b1-3', 'b1-3']
}

# serotype 33A - No structure
ser_33A = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wciC', 'wciD', 'wciE', 'wciF', 'wzy', 'wzx', 'wciG', 'glf', 'wcjE', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 33B
ser_33B = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciN', 'wciO', 'wcrC', 'wciD', 'wciE', 'wciF', 'wzy', 'wzx', 'wciG', 'glf', 'aliA'],
    'sugars': ['glucose', 'ribitol', 'galactose', 'N-acetylgalactosamine', 'galactose'],
    'side branches': ['', '', 'P-galactose', '', ''],
    'modifications': ['', '', '', '', ''],
    'bonds': ['5-P-6', 'a1-2', 'b1-4', 'b1-3']
}

# serotype 33C - No structure
ser_33C = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wcjG', 'wciN', 'wcrO', 'wcrC', 'wcrD', 'wciF', 'wzy', 'wzx', 'glf', 'wcyO', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 33D - No structure
ser_33D = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciN', 'wciO', 'wcrC', 'wciD', 'wciF', 'wzy', 'wzx', 'wciG', 'glf', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 34
ser_34 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wzy', 'wcrO', 'wcrC', 'wcrD', 'wzx', 'glf', 'wcyO', 'aliA'],
    'sugars': ['glucose', 'galactose', 'ribitol', 'galactose', 'galactose'],
    'modifications': ['', '', '', '', '6Ac0.5'],
    'bonds': ['b1-3', '5-P-3', 'a1-2', 'b1-3']
}

# serotype 35F - No structure
ser_35F = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wcjH', 'wciB', 'wzy', 'wcrO', 'wcrC', 'wcrD', 'wzx', 'wciG', 'glf', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 35A
ser_35A = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wzy', 'wcrI', 'wcrJ', 'wcrK', 'mnp1', 'wcrH', 'mnp2', 'wzx', 'wciG', 'glf', 'wcjE', 'tnp', 'tnp', 'aliA'],
    'sugars': ['glucose', 'galactose', 'galactose', 'mannirol', 'galactose'],
    'modifications': ['', '5,6Ac2', '', '', '2Ac'],
    'bonds': ['b1-3', 'b1-3', '6-P-3', 'b1-1']
}

# serotype 35B
ser_35B = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wzy', 'wcrJ', 'wcrM', 'wcrH', 'tnp', 'wzx', 'wciG', 'glf', 'aliA'],
    'sugars': ['glucose', 'galactose', 'N-acetylgalactosamine', 'ribitol', 'galactose'],
    'modifications': ['', '', '', '', '2Ac0.7'],
    'bonds': ['b1-3', 'b1-6', '5-P-4', 'b1-1']
}

# serotype 35C - No structure
ser_35C = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wzy', 'wcrI', 'wcrJ', 'wcrK', 'mnp1', 'wcrH', 'mnp2', 'wzx', 'wciG', 'glf', 'wcjE', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 36 - No structure
ser_36 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wcrO', 'wchO', 'wcjA', 'mnaA', 'wciF', 'wzx', 'wzy', 'glf', 'wcrH', 'tnp', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 37
ser_37 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wciC', 'wciD', 'wciE', 'wciF', 'wzy', 'wzx', 'wciG', 'glf', 'wcjE', 'aliA'],
    'sugars': ['glucose', 'glucose'],
    'modifications': ['', ''],
    'bonds': ['b1-2']
}

# serotype 38 - No structure
ser_38 = {
    'genes': ['tnp', 'glf', 'wzd', 'wze', 'tnp', 'wzg','wzh', 'wciI', 'wcyA', 'wzy', 'wcyB', 'wcyC', 'wcyD', 'wcyV', 'wzx', 'wcyF', 'gla', 'ugd', 'tnp', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 39 - No structure
ser_39 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wcjH', 'wciB', 'whaI', 'wciE', 'wcrC', 'wcrD', 'wciF', 'wzy', 'wcrG', 'wzx', 'glf', 'wcyO', 'aliA'] #,
#    'sugars': ['', '', '', ''],
#    'modifications': ['', '', '', ''],
#    'bonds': ['', '', '', '']
}

# serotype 40 - No structure
ser_40 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wcwI', 'wcwL', 'wcwK', 'wcxU', 'wzy', 'rbsF', 'wzx', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 41F - No structure
ser_41F = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wcrP', 'wcrQ', 'wcrR', 'wzy', 'wcrW', 'wzx', 'wcrX', 'ugd', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'glf', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': [n]
}

# serotype 41A - No structure
ser_41A = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wcrP', 'wcrQ', 'wcrR', 'wzy', 'wcrW', 'wzx', 'wcrX', 'ugd', 'glf', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 42 - No structure
ser_42 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wciB', 'wzy', 'wcrI', 'wcrJ', 'wcrK', 'mnp1', 'wcrH', 'mnp2', 'wzx', 'wciG', 'glf', 'wcjE', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 43 - No structure
ser_43 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wcjH', 'wciB', 'whaI', 'wciE', 'wcrC', 'wcyN', 'wcrH', 'wcyO', 'wcyP', 'wzx', 'glf', 'wcjE', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 44 - No structure
ser_44 = { 
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wciI', 'wciJ', 'wcxB', 'wzy', 'wcxD', 'wcxE', 'wcxF', 'wzx', 'mnaB', 'mnaA', 'fnlA', 'fnlB', 'fnlC', 'tnp', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 45
ser_45 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wciI', 'wciJ', 'wcxB', 'wciL', 'wcyQ', 'wcxS', 'wzx', 'wzy', 'wcyR', 'gct', 'fnlA', 'fnlB', 'fnlC', 'HG273', 'tnp', 'rmlA', 'rmlC', 'rmlB', 'rmlD', 'tnp', 'aliA'],
    'sugars': ['N-acetylgalactosamine', 'N-acetylfucosamine', 'galactose', 'rhamnose'],
    'side branches': ['', 'P-N-acetylglucosamine', 'P-galactose', ''],
    'modifications': ['', '', '', ''],
    'bonds': ['a1-3', 'a1-3', 'a1-3']
}

# serotype 46 - No structure
ser_46 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'tnp', 'wciI', 'wciJ', 'wcxB', 'wzy', 'wcxD', 'wcxE', 'wcxF', 'wzx', 'mnaB', 'mnaA', 'fnlA', 'fnlB', 'fnlC', 'tnp', 'tnp', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 47F - No structure
ser_47F = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wcjH', 'wciB', 'wzy', 'whaI', 'wcrC', 'wcrD', 'wzx', 'wciG', 'glf', 'wcjE', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 47A - No structure
ser_47A = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wcjG', 'wciB', 'whaI', 'wcrC', 'wcyM', 'wzy', 'wcyN', 'whaM', 'wzx', 'glf', 'wcjE', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}

# serotype 48 - No structure
ser_48 = {
    'genes': ['wzg', 'wzh', 'wzd','wze', 'wchA', 'wchF', 'wcxG', 'abp1', 'abp2', 'wcyS', 'wcwY', 'wzy', 'wzx', 'glf', 'rmlA', 'rmlC','rmlB','rmlD', 'aliA'] #,
#    'sugars': [''],
#    'modifications': [''],
#    'bonds': ['']
}


# In[ ]:





# In[ ]:



