# Data Collection

I collected data from 4 databases:
1. [AqSolDB](https://www.nature.com/articles/s41597-019-0151-1) dataset. It contains 9982 unique compound. The dataset contained a SMILES column
2. [Delaney-processed](https://github.com/deepchem/deepchem/blob/master/datasets/delaney-processed.csv) dataset. It contains 1128 unique compound. The dataset contained a SMILES column
3. [T. J. Hou, K. Xia, W. Zhang, and X. J. Xu 2004](https://pubmed.ncbi.nlm.nih.gov/14741036/) dataset. It contains 1290 unique compound. The dataset did not contain a SMILES coulumn just CAS nuber column, so I converted the CAS number to the crossponding SMILES representation.
4. [Junmei Wang, George Krudy, Tingjun Hou, Wei Zhang, George Holland, and Xiaojie Xu 2007](https://pubs.acs.org/doi/full/10.1021/ci700096r?casa_token=JQmRIoLinscAAAAA%3AMc43u7KVHSyUOa_fW6BfMvKWnoQ3HnPFpgvKqkb4E1daBmN83DPibzG67zAEDf1-GkG9GExjCEYkbvE-) dataset. It contains 1708 unique compound. For some reason the SMILES provided in this dataset were incorrect and RDKit did not read it, but I found a [sdf file](http://modem.ucsd.edu/adme/databases/databases_logS.htm) which I used to generated RDKit elements then valid SMILES.

# Data Curation
- I organized the data inside the 4 datasets to have the same forum which consists of 19 columns. 18 columns are the numerical data and a column for SMILES.

| SMILES | Solubility | MolWt | MolLogP | MolMR | HeavyAtomCount | NumHAcceptors | NumHDonors | NumHeteroatoms | NumRotatableBonds | NumValenceElectrons | NumAromaticRings | NumSaturatedRings | NumAliphaticRings | RingCount | TPSA | LabuteASA | BalabanJ | BertzCT |
|--------|------------|------:|--------:|------:|---------------:|--------------:|-----------:|---------------:|------------------:|--------------------:|-----------------:|------------------:|------------------:|----------:|-----:|----------:|---------:|--------:|

- I then merged the 4 datasets together into one big dataset. It contained 13934 molecule, then I removed duplicated moleules based on SMILES column to have 11714 compound. 
- But there isone problem with SMILES that some compounds can have many valied SMILES, so it good to change SMILES to canonical smiles then droup duplicates, after this I got 10379 unique SMILES which are then saved to the mega dataset that will be used then for the model training.
