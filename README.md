### hs37d5: Tools and IGV enhancements

IGV natively supports hg19 and hg38, but not hs37d5 which is a fairly popular reference assembly. This repo is to assist those working with this assembly. 

For dbSNP to work you need to run the following:

```
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz \                                                       INT ✘  53m 22s    13:41:30 
| gunzip -c \
| awk 'BEGIN{OFS="\t"} {c=$2; sub(/^chr/,"",c); print c,$3,$4,$5,$6,$7}' \
| bgzip -c > snp151_hs37d5.bed.gz

tabix -p bed snp151_hs37d5.bed.gz
```
