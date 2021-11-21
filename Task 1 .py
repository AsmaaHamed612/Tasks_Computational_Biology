#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Task number 11111111111


# In[2]:


import pyopenms as ms
sum = 0
seq=ms.AASequence.fromString("VAKA")
for l in seq:
    sum +=l.getMonoWeight()
sum
    


# In[3]:


seq=ms.AASequence.fromString("VAKA")
total=seq.getMonoWeight()
total


# In[4]:


total==sum


# In[ ]:




