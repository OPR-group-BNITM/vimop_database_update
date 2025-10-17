# Host database for ViMOP documentation

The host data base is manually created.
Add all the genomes that you want to use as filters as gzipped fasta files into a directory and create a yaml file with the following format:

```yaml
filters:
  reagent: "reagent-db.fasta.gz"
  human_rna: "GCF_000001405.39_GRCh38.p13_rna.fna.gz"
  human_dna: "GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
  mouse: "GCA_000001635.8_GRCm38.p6_genomic.fna.gz"
  mastomys: "GCF_008632895.1_UCSF_Mcou_1_genomic.fna.gz"
version: 1.0
description: "Human (GRCh38), mouse (8_GRCm38), mastomys and contaminant filter set"
```
