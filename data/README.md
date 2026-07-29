# Data setup

Data files are intentionally excluded from version control. GitHub
blocks ordinary Git files larger than 100 MiB, and the RD file required
here is approximately 430 MiB.

## Gene expression

Expected local path:

```text
data/processed/GSE9101_gene_expression_matrix.csv
```

The matrix used in the analysis has 20,989 genes (rows) and 12 samples
(columns), with the 3 unstimulated controls first. The source experiment
is NCBI GEO accession GSE9101:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9101

The processed matrix is available at:
https://github.com/wyling01/TwoSampleEPB-Paper/tree/main/Real%20Dataset/GSE9101/GSE9101%20Matrix

## Regression discontinuity

Expected local path:

```text
data/raw/rd/data-AER-3.dta
```

Obtain the replication package for:

Pop-Eleches and Urquiola (2013), "Going to a Better School: Effects and
Behavioral Responses," American Economic Review 103(4), 1289-1324.

https://pubs.aeaweb.org/doi/10.1257/aer.103.4.1289

From the replication package, copy `data-AER-3.dta` to the expected
path.
