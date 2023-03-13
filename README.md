# VCFtoProteogenomics

VCFtoProteogenomics is a Python script that generates a protein FASTA file out of a variant call format (VCF) file. The script makes assumptions about the VCF file format and requires certain annotation files to function properly.

## VCF Assumptions:
```
- The Indels in the VCf are reported as following:
      Insertion: 1    125423    A    AGCC,AT
      Insertion: 1    125423    A    AGCC
      Deletion:  1    125423    ACC    A
```
## How to use:

### 1. Prepare the annotation files:

Download the following annotation files from https://www.ensembl.org/info/data/ftp/index.html, using any version of ENSEMBL you prefer:

    - CDS (FASTA)
    - Gene sets (GTF)

### 2. Prepare the python environment
```
# use conda to create the environment with the required dependencies
conda env create -f ../envs/VCFtoProteogenomics.yml
```
Activate the environement
```
conda activate VCFtoProteogenomics
```

### 3. Launch the script for conversion
To generate a protein FASTA file, use the following command:
```
python src/VCFtoProteogenomics.py \
    -input_gtf /path/to/input_gtf \
    -input_cds /path/to/input_cds \
    -input_vcf /path/to/input_vcf \
    -output /path/to/output

```

Note that you should replace **/path/to/input_gtf**, **/path/to/input_cds**, **/path/to/input_vcf**, and **/path/to/output** with the actual paths to the input and output files.

#### Example Use Case
Here's an example command using variables:

```
gtf="/path/to/gtf_file.gtf"
cds="/path/to/cds_file.fasta"
vcf="/path/to/vcf_file.vcf"
output="/path/to/output.fasta"

conda activate VCFtoProteogenomics
python src/VCFtoProteogenomics.py \
    -input_gtf $gtf \
    -input_cds $cds \
    -input_vcf $vcf \
    -output $output
```

Author: 

    - Georges BEDRAN, gbadran_90@live.com
