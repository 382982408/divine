======================================================================
DIVINE: [DI]sease-causing genes/[V]ariant prioritization in cl[IN]ical whole [E]xome sequencing data

Author: Changjin Hong, Ph.D (changjin.hong@gmail.com)
Version: 0.1.1
Website: https://github.com/cjhong
License: GNU GENERAL PUBLIC LICENSE
======================================================================

1. Input:
-Human Phenotype Ontology (HPO) IDs describing patient clinical features
-or/and VCF file

2. Objectives:
-Divine is designed to make molecular diagnosis with clinical whole exome sequencing data more effective. Using both patient phenotypic information and genetic variants, Divine applies a machine learning technique to systematic multi-omic data integration such that it can prioritize a gene list causing a patient disease.

3. Output:
-Divine output is deterministic.
-If only HPO IDs are given, Divine generates a prioritized gene list and inferred diseases.
-If VCF file (or with HPO IDs), Divine generates an annotated variant table with a ranking score.

4. Developer note:
This is the README file for the DIVINE program. DIVINE is in active development.
Please check the above website for newer available versions. Please contact the
developers at changjin.hong@gmail.com to report any problems or for additional help.

5. Dependency:
-Prerequisite: Linux, python 2.7, pip, grep, awk, sort, rsync, and internet connection
-Python modules: Divine requires python modules. The setup script will install the following modules if necessary.

*fastsemsim-0.9.4: https://sites.google.com/site/fastsemsim
*hgvs: https://github.com/counsyl/hgvs
*hpo_similarity: https://github.com/jeremymcrae/hpo_similarity
ConfigParser
backports-abc
backports.ssl-match-hostname
certifi
decorator
matplotlib
networkx
nose
numpy
pandas
pygr
pyparsing
pysam
python-dateutil
pytz
scikit-learn
singledispatch
six
tornado
xlwt

[*] the module is already included in the package.

6. Installation

6.1. github
	Download divine from www.github.com/cjhong/divine

6.2. run the following installation script
	setup_divine.py --install --update_db
	
	*note that it will take very long to download/extract all resource files (16GB)
	
6.3 setup shell configuration after setup

	Get environment variables (DIVINE,GCN,GCN_DATA_DIR,GCN_DB_DIR,GCN_LOGFILE,PATH,PYTHONPATH) at the end of the installation setup message

	Add/modify the environment variables into your shell configuration file ($HOME/.bash_profile or $HOME/.profile).

6.4 apply the new shell configuration
	Logoff/login to your account to make sure that the new configuration is on
	or run "source $HOME/.bashrc"
	or run "source $HOME/.bash_profile"
	or run "source $HOME/.profile"

7. FAQ

7.1.Q:
From a patient (e.g, say, patient ID is P0001), I have both patient phenotypic description (clinical features) and VCF file (e.g, P0001.vcf).

7.1.A: 
First, visit either "http://compbio.charite.de/phenomizer" or "https://mseqdr.org/search_phenotype.php". Enter the patient phenotype terms, description, or specific keywords you think important. Get the best matching HPO IDs. Paste the HPO IDs in the format of HP:XXXXXXX (e.g., HP:0002307) into a text file line by line and save it (e.g., P0001.hpo). Refer to a sample HPO input file in $DIVINE/gcn/bin/prioritize/example. Then, run this.

$DIVINE/gcn/bin/prioritize/divine.py -q dir_to_the_hpo/P0001.hpo -v dir_to_the_vcf/P0001.vcf -o dir_to_the_output/P0001
----------------
7.2.Q:
I have both phenotype file (e.g., P0002.hpo) containing HPO IDs and VCF file (e.g., P0002.vcf).

7.2.A:
$DIVINE/gcn/bin/prioritize/divine.py -q dir_to_the_hpo/P0002.hpo -v dir_to_the_vcf/P0002.vcf -o dir_to_the_output/P0002
----------------
7.3.Q:
Patient phenotypic information is unaviable and I have only VCF file (e.g., P0003.vcf).

7.3.A:
$DIVINE/gcn/bin/prioritize/divine.py -v dir_to_the_vcf/P0003.vcf -o dir_to_the_output/P0003
----------------
7.4.Q:
I have patient phenotypic information but I don't have a VCF file. Or, I only need a set of genes to work with. Refer to 7.1.A to obtain HPO IDs (e.g., P0004.hpo).

7.4.A:
$DIVINE/gcn/bin/prioritize/divine.py -q dir_to_the_hpo/P0004.hpo -o dir_to_the_output/P0004
----------------
7.5.Q:
I have a VCF file (e.g., P0005.vcf) generated from either WGS or WES dataset or the VCF file is too big! I think that Divine is slow and I want to have a result as soon as possible.

7.5.A:
We are actively working on Divine to speed up the process. However, you can make the process faster by masking only NCBI RefSeq genes with +/-20bp flanking from each exon boundary. Be aware that it may not detect a pathogenic variant in UTR or intergenic region.

$DIVINE/gcn/bin/prioritize/divine.py -q dir_to_the_hpo/P0005.hpo -v dir_to_the_vcf/P0005.vcf --mask -o dir_to_the_output/P0005
----------------
7.6.Q:
We purchase HGMD professional license and we want to use the database.

7.6.A:
Contact to changjin.hong@gmail.com
----------------
7.7.Q
I don't like a default filtering scheme (e.g., Allowing both PASS and Low in Filter; Not to use ExAC; 0.03 for MAF cutoff) used in Divine and I want my own filtering strategy.

7.7.A
open $DIVINE/gcn/config/filterconf.txt and edit the configuration file
----------------
7.8.Q
Previously, I ran Divine on a patient sample dataset. I modify the annotation configuration file and I want to repeat the same process. I am sure that I didn't update Divine/GCN database. I want to bypass the annotation process to make my run faster.

7.8.A
$DIVINE/gcn/bin/prioritize/divine.py -q dir_to_the_hpo/P0005.hpo -v dir_to_the_vcf/P0005.vcf --reuse_varant -d rev2 -o dir_to_the_output/P0005
----------------
7.9.Q
I want to use Divine to clinical sequencing dataset. For the same input and configuration, I want Divine to generate the same output. Is Divine a deterministic software?

7.9.A
Yes, Divine is deterministic. There is no randomness in the analysis. Also, Divine generates all log files so that you can audit/trace the previous experiment results and database maintenance any time.
----------------
7.10.Q
I think database or annotation Divine provides is outdate. I would like to keep all the latest database so that I can improve my diagnosis result at best.

7.10.A
Annotation framework Divine currently uses is "varant". We will frequently update the database as each dataset is archieved. Furthermore, we plan to develop a stream pipeline to update annotation database.
----------------
7.11.Q
I know Divine is designed for germline or constitutional disease samples. Can I use Divine for somatic mutation (cancer sample) analysis?

7.11.A
We will work on this feature as well.
----------------
7.12.Q
I have a patient dataset but I am concerning if Divine collects the patient dataset without my permission and send it somewhere.

7.12.A
This standalone package never collect any user information or input dataset the user works. But, we need your feedback and bug report!

----------------
7.13.Q
What are the output files? How I can interpret them?

7.13.A
Depending on the options user provides, Divine generates 3 to 6 output files that can help you on the diagnosis procedure.

-hpo_to_diseases.tsv: From an input HPO file, Divine prioritize which disease the patient likely has. The output format is
#omim	genes	score[BMA]
OMIM:101600	FGFR1,FGFR2	0.000911782
OMIM:101200	FGFR2	0.000674322
:
:

-rank_pheno_gene.tsv: More than a single gene are associated with more than two diseases. This output reformat "hpo_to_diseases.tsv" by sorting the phenotype matching score by gene. The output format is
#gene   phenotypic_score
FGFR1   0.000860439
FGFR2   0.000860439
FGFR3   0.000628328
TWIST1  0.000605784

-divine.vcf: This is a VCF file annotated by Varant where we improve the original method significantly. Refer to Varnat website (http://compbio.berkeley.edu/proj/varant) and our manuscript for more detail.

-vfilter.vcf: This is a subset of the VCF file above (i.e., divine.vcf). This VCF file contains variants after filtering out all variants that is either located in an intergenic region or its MAF is high or in repeated region.

-rank.tsv: Ranking score assigned to each gene appeared into vfilter.vcf

-run_vcf2xls.xls: This is an excel table by generated by Varant. Divine annotates a ranking score per gene so that you can sort the variant list by the ranking score. 

----------------
7.14.Q
Does Divine analyze trio samples?

7.14.A
Yes, if VCF file contains multiple samples and trio analysis can be done in run_vcf2xls.xls file generated by Varant.
----------------
7.15.Q
How long does it take to run Divine? What about memory usage?

7.15.A
It depends on how many variants are contained in an input VCF file. For example, it may take 25 to 30 minute to handle a VCF with 370,000 variants. In the Varant annotation step, memory usage peaks around 3 to 4GB.
----------------

8. Divine Help
usage: divine.py [-h] [-q HPO_QUERY_FN] [-v VCF] [-o OUT_DIR] [-d EXP_TAG]
                 [-I INDEL_MODE] [-r TOP_SEED_RATE] [--wes_mask] [--no_cadd]
                 [--hgmd] [-k VKNOWN] [--reuse_varant] [-c CAPKIT]

Divine (v0.1.1) [author:changjin.hong@gmail.com]

optional arguments:
  -h, --help        show this help message and exit
  -q HPO_QUERY_FN   Input patient hpo file. A file contains HPO IDs (e.g.,
                    HP:0002307), one entry per line. Refer to
                    http://compbio.charite.de/phenomizer or
                    https://mseqdr.org/search_phenotype.php
  -v VCF            input vcf file
  -o OUT_DIR        output directory without white space. If not exist, it
                    will create the directory for you.
  -d EXP_TAG        specify experiment tag without white space. The tag will
                    be contained in the output file name.
  -I INDEL_MODE     the level of fidelity of indell call in VCF [1]:low (e.g.,
                    samtools), 2:high (GATK haplotype caller)
  -r TOP_SEED_RATE  the rate of choosing matched diseases from top for gene
                    enrichment [0.0015]; set to 0. to disable
  --wes_mask        to make the annotation process faster; the annotation
                    process only runs on RefSeq coding regions [False]
  --no_cadd         disable CADD [False]
  --hgmd            enable HGMD (requires a license) [False]
  -k VKNOWN         apply variant-level pathogenic annotation (e.g., either
                    ClinVar or HGMD) to prioritization strategy [1:Yes], 0:No
  --reuse_varant    Reuse previous annotation file (varant.vcf) if it is
                    available [False]
  -c CAPKIT         capture kit symbol[SureSelect_V6]



9. Example

cd $DIVINE/gcn/bin/prioritize/examples
./runme_pfeisffer.sh
./runme_millerSyndrome.sh
./runme_angelman.sh

10. Change Log
v.0.1.1 (June 15 2016)
- Original release

11. Contanct
Changjin Hong, Ph.D (changjin.hong@gmail.com)
