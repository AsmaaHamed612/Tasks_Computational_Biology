#!/usr/bin/env python
# coding: utf-8

# In[1]:


# code number 111111    CONSTANTS

import pyopenms


# In[2]:


help(pyopenms.Constants)


# In[3]:


print ("Avogadro's number is", pyopenms.Constants.AVOGADRO)


# In[8]:


# code number 222       ELEMENTS

from pyopenms import *


# In[9]:


edb = ElementDB()
edb.hasElement("O")
edb.hasElement("S")
oxygen = edb.getElement("O")
print(oxygen.getName())
print(oxygen.getSymbol())
print(oxygen.getMonoWeight())
print(oxygen.getAverageWeight())


# In[10]:


sulfur = edb.getElement("S")
print(sulfur.getName())
print(sulfur.getSymbol())
print(sulfur.getMonoWeight())
print(sulfur.getAverageWeight())
isotopes = sulfur.getIsotopeDistribution()


# In[11]:


print ("One mole of oxygen weighs", 2*oxygen.getAverageWeight(), "grams")
print ("One mole of 16O2 weighs", 2*oxygen.getMonoWeight(), "grams")


# In[12]:


#  code number 333333    Isotopes


edb = ElementDB()
oxygen_isoDist = {"mass": [], "abundance": []}
sulfur_isoDist = {"mass": [], "abundance": []}


# In[13]:


oxygen = edb.getElement("O")
isotopes = oxygen.getIsotopeDistribution()


# In[14]:


for iso in isotopes.getContainer():
    print ("Oxygen isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    oxygen_isoDist["mass"].append(iso.getMZ())
    oxygen_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[15]:


sulfur = edb.getElement("S")
isotopes = sulfur.getIsotopeDistribution()
for iso in isotopes.getContainer():
    print ("Sulfur isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    sulfur_isoDist["mass"].append(iso.getMZ())
    sulfur_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[16]:


# code number 444444    Amino Acids

lys = ResidueDB().getResidue("Lysine")


# In[17]:


print(lys.getName())


# In[18]:


print(lys.getThreeLetterCode())


# In[19]:


print(lys.getAverageWeight())


# In[20]:


print(lys.getMonoWeight())


# In[21]:


print(lys.getPka())
print(lys.getFormula().toString())


# In[22]:


# code number 555555   Molecular Formulae

methanol = EmpiricalFormula("CH3OH")


# In[23]:


water = EmpiricalFormula("H2O")


# In[24]:


ethanol = EmpiricalFormula("CH2") + methanol


# In[25]:


print("Ethanol chemical formula:", ethanol.toString())
print("Ethanol composition:", ethanol.getElementalComposition())
print("Ethanol has", ethanol.getElementalComposition()[b"H"], "hydrogen atoms")


# In[26]:


# code number 666666 Isotopic Distributions

methanol = EmpiricalFormula("CH3OH")
ethanol = EmpiricalFormula("CH2") + methanol


# In[27]:


methanol_isoDist = {"mass": [], "abundance": []}
ethanol_isoDist = {"mass": [], "abundance": []}


# In[28]:


print("Coarse Isotope Distribution:")
isotopes = ethanol.getIsotopeDistribution( CoarseIsotopePatternGenerator(4) )
prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
print("This covers", prob_sum, "probability")


# In[29]:


for iso in isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    methanol_isoDist["mass"].append(iso.getMZ())
    methanol_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[30]:


print("Fine Isotope Distribution:")
isotopes = ethanol.getIsotopeDistribution( FineIsotopePatternGenerator(1e-3) )
prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
print("This covers", prob_sum, "probability")


# In[31]:


for iso in isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    ethanol_isoDist["mass"].append(iso.getMZ())
    ethanol_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[35]:


# code number 77777   Amino Acid Sequences

seq = AASequence.fromString("DFPIANGER") 
prefix = seq.getPrefix(4)
suffix = seq.getSuffix(5) 
concat = seq + seq


# In[36]:


print("Sequence:", seq)
print("Prefix:", prefix)
print("Suffix:", suffix)
print("Concatenated:", concat)


# In[37]:


mfull = seq.getMonoWeight() 


# In[38]:


mprecursor = seq.getMonoWeight(Residue.ResidueType.Full, 2)


# In[39]:


mz = seq.getMZ(2)


# In[40]:


print()
print("Monoisotopic mass of peptide [M] is", mfull)
print("Monoisotopic mass of peptide precursor [M+2H]2+ is", mprecursor)
print("Monoisotopic m/z of [M+2H]2+ is", mz)


# In[41]:


seq = AASequence.fromString("DFPIANGER")

print("The peptide", str(seq), "consists of the following amino acids:")
for aa in seq:
    print(aa.getName(), ":", aa.getMonoWeight())


# In[42]:


# code number 88888  N- or C-terminal modification

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


# In[43]:


# code number 99999 Molecular formula

seq = AASequence.fromString("DFPIANGER")
seq_formula = seq.getFormula()
print("Peptide", seq, "has molecular formula", seq_formula)


# In[44]:


# code number 10 10 10   Fragment ions




suffix = seq.getSuffix(3)
print("="*35)
print("y3 ion sequence:", suffix)
y3_formula = suffix.getFormula(Residue.ResidueType.YIon, 2) 
suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 
suffix.getMonoWeight(Residue.ResidueType.XIon, 2) / 2.0 
suffix.getMonoWeight(Residue.ResidueType.BIon, 2) / 2.0 

print("y3 mz:", suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 )
print("y3 molecular formula:", y3_formula)


# In[45]:


# code number 11 11 11     Modified Sequences


seq = AASequence.fromString("PEPTIDESEKUEM(Oxidation)CER")


# In[47]:


print(seq.toUnmodifiedString())
print(seq.toString())
print(seq.toUniModString())
print(seq.toBracketString())
print(seq.toBracketString(False))


# In[48]:


print(AASequence.fromString("DFPIAM(UniMod:35)GER"))
print(AASequence.fromString("DFPIAM[+16]GER"))
print(AASequence.fromString("DFPIAM[+15.99]GER"))
print(AASequence.fromString("DFPIAM[147]GER"))
print(AASequence.fromString("DFPIAM[147.035405]GER"))


# In[3]:


# code number 12 12 12 Spectrum

from pyopenms import *


# In[4]:


spectrum = MSSpectrum()

mz = range(1500, 500, -100)

i = [0 for mass in mz]

spectrum.set_peaks([mz, i])


# In[5]:


spectrum.sortByPosition()


# In[6]:


for p in spectrum:

  print(p.getMZ(), p.getIntensity())


# In[7]:


for mz, i in zip(*spectrum.get_peaks()):

  print(mz, i)


# In[8]:


print(spectrum[2].getMZ(), spectrum[2].getIntensity())


# In[9]:


help(MSSpectrum)


# In[ ]:




# code number 14 14 14 14 Fasta File


# In[32]:


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

entries = []
f = FASTAFile()
f.load("example.fasta", entries)
print( len(entries) )
for e in entries:
    print (e.identifier, e.sequence)


# In[30]:





# In[ ]:





# In[ ]:





# In[ ]:




