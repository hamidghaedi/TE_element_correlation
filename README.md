### 1) Raw count to vst normalized values in R :

```R
# read repeat masker data
# getting TEs
rmk <- readRDS("~/LINE1-BLCA/rmsk_annotation.RDS")
intrestedElelements <- rmk$repName[rmk$repClass %in% c("LINE", "SINE", "LTR","DNA", "Retroposon")]

# save a subset for later 
rmkSub <- rmk[rmk$repClass %in% c("LINE", "SINE", "LTR","DNA", "Retroposon"),]

# Reading data 
# clinical data
clinicalU <- data.frame(data.table::fread("~/LINE1-BLCA/UROMOL_final/clinU.csv"))

# all loci
teExpU <- readRDS("~/LINE1-BLCA/UROMOL_final/RE_all_1_raw_counts.RDS")

# reading expression data for each loci

# intergenic
IntGenTeExpU <- readRDS("~/LINE1-BLCA/UROMOL_final/RE_intergenic_1_raw_counts.RDS")
IntGenTeExpU <- IntGenTeExpU$counts
IntGenTeExpU <- IntGenTeExpU[rownames(IntGenTeExpU) %in% intrestedElelements,]

# exonic
exTeExpU <- readRDS("~/LINE1-BLCA/UROMOL_final/RE_exon_1_raw_counts.RDS")
exTeExpU <- exTeExpU$counts
exTeExpU <- exTeExpU[rownames(exTeExpU) %in% intrestedElelements,]

# intronic
IntronTeExpU <- readRDS("~/LINE1-BLCA/UROMOL_final/RE_intron_1_raw_counts.RDS")
IntronTeExpU <- IntronTeExpU$counts
IntronTeExpU <- IntronTeExpU[rownames(IntronTeExpU) %in% intrestedElelements,]

# vst normalization
library(DESeq2)

# list of count matrices
exp <- list(teExpU, exTeExpU, IntGenTeExpU, IntronTeExpU)
# empty list to be populated in loop
vst_res <- list()

for(i in 1:length(exp)){
  # preparing expression sets
  expDatU = exp[[i]]
  clinU = clinicalU
  rownames(clinU) <- clinU$UROMOL.ID

  # to see whether all rows of clinical are present in expression datset
  #all(rownames(clinU) %in% colnames(expDatU))
  # whether they are in the same order:
  all(rownames(clinU) == colnames(expDatU))
  # reorder expression df columns based on sample orders in clinical file
  expDatU <- expDatU[, rownames(clinU)]
  #all(rownames(clinU) == colnames(expDatU))
  
  # Making DESeqDataSet object which stores all experiment data
  dds <- DESeqDataSetFromMatrix(countData = round(expDatU),
                                colData = clinU,
                                design = ~ 1)
  
  # prefilteration: it is not necessary but recommended to filter out low expressed genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  # data tranfromation
  vsdU <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
  vst_res[[i]] <- assay(vsdU)
}

# saveing results into csv files:
write.csv(vst_res[[1]], "teExpU_vst.csv")
write.csv(vst_res[[2]], "exTeExpU_vst.csv")
write.csv(vst_res[[3]], "IntGenTeExpU_vst.csv")
write.csv(vst_res[[4]], "IntronTeExpU_vst.csv")

# exporting repeat masker ddata
d <- rmkSub[!duplicated(rmkSub$repName), ]
write.csv(d[c(3,5,4)], "rmk.csv")

```

### 2) Correlation analysis in Python:

```python
import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


# rading files
os.chdir("C:\\Users\\qaedi\\OneDrive - Queen's University\\Documents\\")

all_loc = pd.read_csv("teExpU_vst.csv", index_col=0)
exon = pd.read_csv("exTeExpU_vst.csv", index_col=0)
intron = pd.read_csv("IntronTeExpU_vst.csv", index_col=0)
intergenic = pd.read_csv("IntGenTeExpU_vst.csv", index_col=0)
rmk = pd.read_csv("rmk.csv")



# to get LINE elements
line = rmk.loc[rmk['repClass'] == 'LINE']['repName'].tolist()
# to get SINE elements
sine = rmk.loc[rmk['repClass'] == 'SINE']['repName'].tolist()
# to get DNA elements
dna = rmk.loc[rmk['repClass'] == 'DNA']['repName'].tolist()
# to get retroposon elements
retroposon = rmk.loc[rmk['repClass'] == 'Retroposon']['repName'].tolist()
# to get LTR elements
ltr = rmk.loc[rmk['repClass'] == 'LTR']['repName'].tolist()


# to create a df with five columns LINE, SINE, DNA, Retroposon and LTR 
def convert_to_TEclass (df):
    #line
    filtered_df = df[df.index.isin(line)] 
    l = list(filtered_df.mean())
    #sine
    filtered_df = df[df.index.isin(sine)] 
    s = list(filtered_df.mean())
    #dna
    filtered_df = df[df.index.isin(dna)] 
    d = list(filtered_df.mean())
    #ret
    filtered_df = df[df.index.isin(retroposon)] 
    r = list(filtered_df.mean())
    #ltr
    filtered_df = df[df.index.isin(ltr)] 
    lt = list(filtered_df.mean())
    co = list(filtered_df.columns)
    res = pd.DataFrame({'LINE': l,
                            'SINE' : s,
                            'DNA' : d,
                            'Retroposon' : r,
                            'LTR' : lt}, index = co)
    return(res)

# converting dfs to class wise expression and calculating correlation
all_loc_class = convert_to_TEclass(all_loc)
exon_class = convert_to_TEclass(exon)
intron_class = convert_to_TEclass(intron)
intergenic_class = convert_to_TEclass(intergenic)
```

### 3) Visualization

```python
# all loci

# Select only the numerical columns
df_num = all_loc_class.select_dtypes(include=[np.number])
# Compute the correlation matrix
corr = df_num.corr()
# Plot the correlation matrix using a heatmap from seaborn
palette = sns.color_palette("RdBu", n_colors=256)
plt.figure(figsize=(10,10))
sns.heatmap(corr, annot=True, cmap=palette)
plt.show()
```
![alt text](https://github.com/hamidghaedi/TE_element_correlation/blob/main/image/all_loci_in_one_file.png)

```python
# exonic

df_num = exon_class.select_dtypes(include=[np.number])
# Compute the correlation matrix
corr = df_num.corr()
# Plot the correlation matrix using a heatmap from seaborn
palette = sns.color_palette("RdBu", n_colors=256)
plt.figure(figsize=(10,10))
sns.heatmap(corr, annot=True, cmap=palette)
plt.show()
```
![alt text](https://github.com/hamidghaedi/TE_element_correlation/blob/main/image/exon.png)

```python
# intronic

# Select only the numerical columns
df_num = intron_class.select_dtypes(include=[np.number])
# Compute the correlation matrix
corr = df_num.corr()
# Plot the correlation matrix using a heatmap from seaborn
palette = sns.color_palette("RdBu", n_colors=256)
plt.figure(figsize=(10,10))
sns.heatmap(corr, annot=True, cmap=palette)
plt.show()

```
![alt text](https://github.com/hamidghaedi/TE_element_correlation/blob/main/image/intron.png)

```python
#intergenic

# Select only the numerical columns
df_num = intergenic_class.select_dtypes(include=[np.number])
# Compute the correlation matrix
corr = df_num.corr()
# Plot the correlation matrix using a heatmap from seaborn
palette = sns.color_palette("RdBu", n_colors=256)
plt.figure(figsize=(10,10))
sns.heatmap(corr, annot=True, cmap=palette)
plt.show()
```
![alt text](https://github.com/hamidghaedi/TE_element_correlation/blob/main/image/intergenic.png)


### 4) Summary table


|            |            | exon   | intron | intergenic | all_in_one_file |
| ---------- | ---------- | ------ | ------ | ---------- | --------------- |
| LINE       | SINE       | \-0.55 | \-0.19 | \-0.44     | \-0.32          |
| LINE       | DNA        | 0.56   | 0.08   | \-0.09     | \-0.06          |
| LINE       | Retroposon | \-0.3  | \-0.45 | \-0.53     | \-0.4           |
| LINE       | LTR        | 0.38   | \-0.5  | \-0.15     | \-0.35          |
| SINE       | DNA        | \-0.46 | \-0.03 | \-0.24     | \-0.22          |
| SINE       | Retroposon | 0.45   | 0.4    | 0.66       | 0.55            |
| SINE       | LTR        | \-0.12 | 0.17   | 0.16       | 0.19            |
| DNA        | Retroposon | \-0.34 | \-0.67 | \-0.45     | \-0.59          |
| DNA        | LTR        | 0.54   | \-0.27 | 0.23       | \-0.4           |
| Retroposon | LTR        | 0.16   | 0.54   | 0.14       | 0.37            |
