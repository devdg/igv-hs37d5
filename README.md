### hs37d5: Tools and IGV enhancements

IGV natively supports hg19 and hg38, but not hs37d5 which is a fairly popular reference assembly. This repo is to assist those working with this assembly. 

If you follow the below the following will work:
1. Cytobands
2. Reference Sequence
3. RefSeq Genes with annotation and protein sequences
4. dbSNP annonations 

First you will need to get teh FASTA for hs37d5. Illumina has this here: https://ilmn-dragen-giab-samples.s3.amazonaws.com/FASTA/hs37d5.fa other options can be found here: https://knowledge.illumina.com/software/on-premises-software/software-on-premises-software-troubleshooting-list/000007409

For dbSNP to work you need to run the following:

```
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz \                                                       INT ✘  53m 22s    13:41:30 
| gunzip -c \
| awk 'BEGIN{OFS="\t"} {c=$2; sub(/^chr/,"",c); print c,$3,$4,$5,$6,$7}' \
| bgzip -c > snp151_hs37d5.bed.gz

tabix -p bed snp151_hs37d5.bed.gz
```

Ensure you have all the files from here in the same directory. Then load the JSON in IGV using "Genomes"->"Load Genome from file..."
