vcf_to_prot is for generating a protein fasta file out of a variant call format (vcf) file.

## VCF Assumptions:
```
- The Indels in the VCf are reported as following:
      Insertion: 1    125423    A    AGCC,AT
      Insertion: 1    125423    A    AGCC
      Deletion:  1    125423    ACC    A
```
## How to use:

### 1. Prepare the annotation files:

please download:

    - CDS (FASTA)
    - Gene sets (GTF)

from <https://www.ensembl.org/info/data/ftp/index.html> using any version of ENSEMBL you prefer


### 2. Prepare the python environment
```
# use conda to create the environment with the required dependencies
conda env create -f ../envs/vcf_to_prot.yml

# to activate the environement
conda activate vcf_to_prot
```

### 3. Launch the script for conversion
```
conda activate vcf_to_prot
$ python vcf_to_prot.py --help
usage: vcf_to_prot.py [-h] -input_gtf INPUT_GTF -input_cds INPUT_CDS
                      -input_vcf INPUT_VCF -permutation {0,1} -output OUTPUT

generates a protein fasta file out of a variant call format (vcf) file

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  -input_gtf INPUT_GTF  Path to input GTF fasta
  -input_cds INPUT_CDS  Path to input cds fasta
  -input_vcf INPUT_VCF  Path to input vcf variant file
  -permutation {0,1}    generate aa permutations for novel stop codons
  -output OUTPUT        Path to output fasta


NOTE: -perumutation argument:
    translate with termination at the second in-frame stop codon.
    Sequences with novel stop codons will be further cleaved to 21 amino acids centering the stop sites which will be substituted with the 20 known amino acids.
    Next, sequences will be unspecifically cleaved using the pyteomics library and the cleaved peptides covering the substituted novel stop codon site within the 8 to 12 amino acids range will only be considered.
    otherwise if you need a regular variant fasta file just ignore it

```
#### usecase example

```
$gtf="/path/to/gtf_file.gtf"
$cds="/path/to/cds_file.fasta"
$vcf="/path/to/vcf_file.vcf"
$output="/path/to/output.fasta"

conda activate vcf_to_prot
python vcf_to_prot.py \
    -input_gtf $gtf \
    -input_cds $cds \
    -input_vcf $vcf \
    -output $ output
```

Author: 

    - Georges BEDRAN, gbadran_90@live.com
