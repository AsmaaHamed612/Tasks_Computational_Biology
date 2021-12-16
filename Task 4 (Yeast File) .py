#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Task number 4


# In[5]:


from pyopenms import *

digestion = ProteaseDigestion()
digestion.getEnzymeName() 
seq = "".join([l.strip() for l in open("uniprot-yourlist_M2021112192C7BAECDB1C5C413EE0E0348724B68230CFBAX.fasta").readlines()[1:]])
seq = AASequence.fromString(seq)

result = []
digestion.digest(seq, result)
print(result[4].toString())
len(result) 


# In[ ]:





# In[ ]:




