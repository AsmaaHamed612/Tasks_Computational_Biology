#!/usr/bin/env python
# coding: utf-8

# In[2]:


#                                            [1] Amino Acid Sequences

from pyopenms import *
seq = AASequence.fromString("DFPIANGER") 
prefix = seq.getPrefix(4) 
suffix = seq.getSuffix(5) 
concat = seq + seq
 
print("Sequence:", seq)
print("Prefix:", prefix)
print("Suffix:", suffix)
print("Concatenated:", concat)

mfull = seq.getMonoWeight() 
mprecursor = seq.getMonoWeight(Residue.ResidueType.Full, 2) 

mz = seq.getMonoWeight(Residue.ResidueType.Full, 2) / 2.0 
mz = seq.getMZ(2) 

print()
print("Monoisotopic mass of peptide [M] is", mfull)
print("Monoisotopic mass of peptide precursor [M+2H]2+ is", mprecursor)
print("Monoisotopic m/z of [M+2H]2+ is", mz)


# In[3]:


seq = AASequence.fromString("DFPIANGER")

print("The peptide", str(seq), "consists of the following amino acids:")
for aa in seq:
    print(aa.getName(), ":", aa.getMonoWeight())


# In[5]:


#                                                   [2] N- and C-Terminus

seq = AASequence.fromString("C[143]PKCK(Label:13C(6)15N(2))CR")
if seq.hasNTerminalModification():
    print("N-Term Modification: ", seq.getNTerminalModification().getFullId())
if seq.hasCTerminalModification():
    print("C-Term Modification: ", seq.getCTerminalModification().getFullId())
for aa in seq:
    if (aa.isModified()):
        print(aa.getName(), ":", aa.getMonoWeight(), ":", aa.getModificationName())
    else:
        print(aa.getName(), ":", aa.getMonoWeight())


# In[6]:


#                                        [3] Molecular formula

seq = AASequence.fromString("DFPIANGER")
seq_formula = seq.getFormula()
print("Peptide", seq, "has molecular formula", seq_formula)


# In[8]:


#                                         [4]Isotope patterns

coarse_isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(6) )
for iso in coarse_isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")


# In[9]:


fine_isotopes = seq_formula.getIsotopeDistribution( FineIsotopePatternGenerator(0.01) )
for iso in fine_isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")


# In[10]:


#                              [5] isotop distribution coarse and fine_structure

import math
from matplotlib import pyplot as plt

def plotIsotopeDistribution(isotope_distribution, title="Isotope distribution"):
    plt.title(title)
    distribution = {"mass": [], "abundance": []}
    for iso in isotope_distribution.getContainer():
        distribution["mass"].append(iso.getMZ())
        distribution["abundance"].append(iso.getIntensity() * 100)

    bars = plt.bar(distribution["mass"], distribution["abundance"], width=0.01, snap=False) 

    plt.ylim([0, 110])
    plt.xticks(range(math.ceil(distribution["mass"][0]) - 2,
                     math.ceil(distribution["mass"][-1]) + 2))
    plt.xlabel("Atomic mass (u)")
    plt.ylabel("Relative abundance (%)")

plt.figure(figsize=(10,7))
plt.subplot(1,2,1)
plotIsotopeDistribution(coarse_isotopes, "Isotope distribution - coarse")
plt.subplot(1,2,2)
plotIsotopeDistribution(fine_isotopes, "Isotope distribution - fine structure")
plt.show()


# In[11]:


#                                             [6]Fragment ions

suffix = seq.getSuffix(3) 
print("="*35)
print("y3 ion sequence:", suffix)
y3_formula = suffix.getFormula(Residue.ResidueType.YIon, 2) 
suffix.getMonoWeight(Residue.ResidueType.XIon, 2) / 2.0 
suffix.getMonoWeight(Residue.ResidueType.BIon, 2) / 2.0

print("y3 mz:", suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 )
print("y3 molecular formula:", y3_formula)


# In[13]:


#                                               [7]Modified Sequences

seq = AASequence.fromString("PEPTIDESEKUEM(Oxidation)CER")
print(seq.toUnmodifiedString())
print(seq.toString())
print(seq.toUniModString())
print(seq.toBracketString())
print(seq.toBracketString(False))

print(AASequence.fromString("DFPIAM(UniMod:35)GER"))
print(AASequence.fromString("DFPIAM[+16]GER"))
print(AASequence.fromString("DFPIAM[+15.99]GER"))
print(AASequence.fromString("DFPIAM[147]GER"))
print(AASequence.fromString("DFPIAM[147.035405]GER"))


# In[19]:


#                                              [8]Proteins and FASTA files

bsa = FASTAEntry() 
bsa.sequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
bsa.description = "BSA Bovine Albumin (partial sequence)"
bsa.identifier = "BSA"
alb = FASTAEntry()
alb.sequence = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGE"
alb.description = "ALB Human Albumin (partial sequence)"
alb.identifier = "ALB"
entries = [bsa, alb]

f = FASTAFile()
f.store("example.fasta", entries)


# In[21]:


entries = []
f = FASTAFile()
f.load("example.fasta", entries)
print( len(entries) )
for e in entries:
    print (e.identifier, e.sequence)


# In[ ]:




