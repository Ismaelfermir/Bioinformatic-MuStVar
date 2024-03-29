# Bioinformatic-MuStVar

Tool in development...

**Bioinformatic-MuStVar** is a tools development in the Instituto de Investigación Sanitario Puerta de Hierro for somatic SNVs and Indels detection by **Mutect2**, **Strelka2** and **VarScan2** from cfDNA samples.

The tool consists of 3 parts:

### 1. MuStVar
The MuStVar.sh file contains the script to call variants from the aligned files (GRCh38-hg38) with the three tools (**Mutect2**, **Strelka2** and **VarScan2**).

First, make sure that your normal and tumor alignmet files are named as:
- Sample1.**normal.bam**
- Sample1.**tumor.bam**

Needs:
1. Computer threads
2. Work directory
3. Reference fasta
4. Panel bed file

Note that the following databases must be downloaded for the ANNOVAR annotation tool:
- refGene
- cytoBand
- 1000g2015aug_all
- 1000g2015aug_eas
- 1000g2015aug_eur
- avsnp150
- cosmic70
- clinvar_20190305

### 2. Filter
The MuStVar.filter.R file contains the script to filter the SNVs and Indels from the three tools by these information:
- Allele frequency > 0.5%
- Altered region: exonic, splicing and UTRs
- Read depth => 500
- Read depth of altered allele: => 30

### 3. Compare
The compare.callers.R file will compare the filtered SNVs and Indels and select just the alterations that are shared by at least two tools. Also it will annotate each alteration with the data base for detecting PCR and panel errors. The script it will report the mutations present in just one variant calling tool.

Note that the final result of alterations may not be 100% correct. A manual check for each mutation is needed. For example: alterations annotated by dbSNP as SNPs will be reported for manual check.

