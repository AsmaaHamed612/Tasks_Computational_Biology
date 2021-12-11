#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Task number 44444444 "Digestion" "Human Example"


# In[2]:


from pyopenms import *

digestion = ProteaseDigestion()
digestion.getEnzymeName() 
seq = "".join([l.strip() for l in open("uniprot-yourlist_M2021121192C7BAECDB1C5C413EE0E0348724B68235CE10N.fasta").readlines()[1:]])
seq = AASequence.fromString(seq)

result = []
digestion.digest(seq, result)
print(result[4].toString())
len(result) 


# In[ ]:




