#!/usr/bin/env python
# coding: utf-8

# In[1]:



import pyopenms as ms
sum = 0
seq=ms.AASequence.fromString("VAKA")
for l in seq:
    sum +=l.getMonoWeight()
sum


# In[2]:


seq=ms.AASequence.fromString("VAKA")
total=seq.getMonoWeight()
total


# In[3]:


total==sum


# In[ ]:




