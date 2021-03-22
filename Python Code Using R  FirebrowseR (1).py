#!/usr/bin/env python
# coding: utf-8

# In[134]:


pip install rpy2


# In[135]:


import rpy2


# In[136]:


import rpy2.robjects as robjects


# In[137]:


print(rpy2.__version__)


# In[138]:


from rpy2.robjects.packages import importr
# import R's "base" package
base = importr('base')

# import R's "utils" package
utils = importr('utils')


# In[139]:


import rpy2.robjects.packages as rpackages

# import R's utility package
utils = rpackages.importr('utils')

# select a mirror for R packages
utils.chooseCRANmirror(ind=1) # select the first mirror in the list


# In[141]:


import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
packageNames = ('afex', 'emmeans', 'ggplot2', 'hexbin', 'FirebrowseR', 'ggplot2', 'ggthemes')
utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]

# Running R in Python example installing packages:
if len(packnames_to_install) > 0:
    utils.install_packages(StrVector(packnames_to_install))


# In[142]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
get_ipython().run_line_magic('Rdevice', 'png')


# In[143]:


get_ipython().run_cell_magic('R', '', '#barplot(c(1,3,2,5,4), ylab="value")\ndirectORY_ <- getwd()\ndirectORY_\nrequire (FirebrowseR)\nrequire(ggplot2)\nrequire(ggthemes)\nmRNA_Exp_PDL1 = Samples.mRNASeq(format = "csv", cohort = "",  gene = c("B7-H", "B7H1", "PD-L1", "PDL1", "B7-H1", "CD274"), protocol = "RSEM", sort_by = "cohort", page = "1", page_size = "500")\nmRNA_Exp_PD1 = Samples.mRNASeq(format = "csv", cohort = "",  gene = c("PDCD1", "CD279", "hSLE1", "PD-1", "PD1"), protocol = "RSEM", sort_by = "cohort", page = "1", page_size = "500")\n\nmRNA_Exp_PDL1\nmRNA_Exp_PD1')


# In[144]:


get_ipython().run_cell_magic('R', '', "listCommon <- intersect(mRNA_Exp_PD1$tcga_participant_barcode,mRNA_Exp_PDL1$tcga_participant_barcode)\nhead(listCommon)\nwrite.table(listCommon, file='table_listCommon.tsv', quote=FALSE, sep='\\t')\nlistCommon")


# In[ ]:


get_ipython().run_cell_magic('R', '', 'listDF <- as.data.frame(listCommon)\ncountex <- nrow(listDF)\nmRNA_Exp_PD1_barcode <- 0\nfor (i in 1:countex){\n  bindingAgent1 = Samples.mRNASeq(format = "csv",  tcga_participant_barcode = listCommon[i], gene = c("PDCD1", "CD279", "hSLE1", "PD-1", "PD1"), protocol = "RSEM")\n  mRNA_Exp_PD1_barcode = rbind(mRNA_Exp_PD1_barcode, bindingAgent1)\n}\nmRNA_Exp_PD1_barcode = mRNA_Exp_PD1_barcode[2:countex,]\nhead(mRNA_Exp_PD1_barcode)')


# In[ ]:


get_ipython().run_cell_magic('R', '', 'mRNA_Exp_PDL1_barcode <- 0\nfor (i in 1:countex){\n  bindingAgent2 = Samples.mRNASeq(format = "csv",  tcga_participant_barcode = listCommon[i], gene = c("B7-H", "B7H1", "PD-L1", "PDL1", "B7-H1", "CD274"), protocol = "RSEM")\n  mRNA_Exp_PDL1_barcode = rbind(mRNA_Exp_PDL1_barcode, bindingAgent2)\n}\nmRNA_Exp_PDL1_barcode = mRNA_Exp_PDL1_barcode[2:countex,]\nhead(mRNA_Exp_PDL1_barcode)')


# In[ ]:


get_ipython().run_cell_magic('R', '', 'totalCohorts <- nrow(mRNA_Exp_PD1_barcode) #Counting total number of participant based on barcodes\nfdata = factor(mRNA_Exp_PD1_barcode$cohort) #Identifying the common cohort from the list\nidentifiedDisease <- levels(fdata)\n\ncohortsList = Metadata.Cohorts(format = "csv")\nnamesOfDiseases <- cohortsList[cohortsList$cohort == identifiedDisease,]$description\nnamesOfDiseases\nwriteReport <- as.data.frame(table(mRNA_Exp_PD1_barcode$cohort))\nfigureThis <- qplot(as.numeric(mRNA_Exp_PD1_barcode$expression_log2), as.numeric(mRNA_Exp_PDL1_barcode$expression_log2), xlab="PD1 expression (Log Scale)", ylab="PDL1 expression (Log Scale)", colour = mRNA_Exp_PDL1_barcode$cohort)\nfigureThis + labs(colour = \'Identified Cohorts =\') + ggtitle("Median RSEM normalized expression / cohort") + theme_stata() + scale_color_stata() + geom_smooth(method = "lm", se = FALSE, linetype="dotted")\n')


# In[158]:


namesOfDiseases


# In[159]:


get_ipython().run_line_magic('R', '-o mRNA_Exp_PDL1 mRNA_Exp_PDL1_barcode')
get_ipython().run_line_magic('R', '-o mRNA_Exp_PD1 mRNA_Exp_PD1_barcode')


# In[152]:



#%R -o yaxisValue mRNA_Exp_PD1_barcode$expression_log2
#yaxisValue
get_ipython().run_line_magic('R', 'mRNA_Exp_PDL1_barcode$expression_log2')
get_ipython().run_line_magic('R', '-o x_Value mRNA_Exp_PDL1_barcode$expression_log2')
x_Value


# In[93]:


mRNA_Exp_PD1


# In[94]:


mRNA_Exp_PD1.info(verbose=True)


# In[95]:


import pandas as pd


# In[116]:


test = mRNA_Exp_PD1.sort_values('expression_log2', ascending=True, na_position='last')


# In[117]:


test


# In[82]:


mRNA_Exp_PD1.info(verbose=True)


# In[112]:


mRNA_Exp_PD1.sort_values(by=['expression_log2'])


# In[ ]:




