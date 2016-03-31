# MLST from SRST2

This tool uses a table of compiled results from [SRST2](https://github.com/katholt/srst2) to create an MLST-like scheme.  Each sample in the file will be assigned a sequence type number based on its combination of alleles for a set of genes.

## Usage:
```
mlst_from_srst2.py [-h] --srst2_table SRST2_TABLE
                   [--existing_scheme EXISTING_SCHEME]
                   [--mlst_genes MLST_GENES] --out_types OUT_TYPES
                   --out_scheme OUT_SCHEME
```

### Input

This script requires a compiled table of results from SRST2.  The first column should be sample names and remaining columns are gene clusters where the cells contain specific alleles.

Example SRST2 input table:

Sample | abcA | abcB | abcC | abcD
:---: | :---: | :---: | :---: | :---:
sample1 | abcA_1 | abcB_1 | abcC_1 | abcD_1
sample2 | abcA_2 | abcB_1 | abcC_1 | abcD_1
sample3 | abcA_2 | abcB_1 | abcC_1 | abcD_3

#### Creating a new scheme

If you do not already have a scheme for a set of genes, you can create one by using the `--mlst_genes` argument: a comma-delimited list of the genes (i.e. SRST2 clusters) you want to use in your scheme.

Example:

`mlst_from_srst2.py --srst2_table input_table.txt --mlst_genes abcA,abcB,abcC --out_types sample_results.txt --out_scheme mlst_scheme.txt`

#### Expanding an existing scheme

If you already have a scheme for a set of genes (produced by a previous run of this script), it is given with the `--mlst_genes` argument.  This scheme will be expanded upon if your samples contain new allele combinations.

Example:

`mlst_from_srst2.py --srst2_table input_table.txt --existing_scheme old_scheme.txt --out_types sample_results.txt --out_scheme new_scheme.txt`

### Output

This script produces two output files: the MLST scheme file and the sequence type assignments for each sample.

#### Example MLST scheme output

ST | abcA | abcB | abcC | abcD
:---: | :---: | :---: | :---: | :---:
1 | abcA_1 | abcB_1 | abcC_1 | abcD_1
2 | abcA_2 | abcB_1 | abcC_1 | abcD_1
3 | abcA_1 | abcB_3 | abcC_1 | abcD_1
4 | abcA_2 | abcB_1 | abcC_1 | abcD_3

#### Example sequence type assignment output

Sample | ST | alleles
:---: | :---: | :---:
sample1 | 1 | abcA_1,abcB_1,abcC_1,abcD_1
sample2 | 2 | abcA_2,abcB_1,abcC_1,abcD_1
sample3 | 4 | abcA_2,abcB_1,abcC_1,abcD_3

Samples will not be given a sequence type if any of the following are true:
* The allele assignment in SRST2 is uncertain or imperfect (contains '?' or '*')
* The allele was not found in SRST2 (is '-')
* The gene was not found in SRST2

## License

GNU General Public License, version 3
