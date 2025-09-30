# Genomic characterization of human respiratory syncytial virus circulating in Islamabad, Pakistan, during an outbreak in 2022 to 2023

Reproducible code and workflow that mirrors the analysis in the paper.  
DOI: https://doi.org/10.1007/s00705-024-06036-0

## Program summary
End to end analysis for RSV-A and RSV-B as in the study. All steps are automated or scripted and documented here so a reviewer can reproduce the flow without guesswork.

1) Inputs
   - Paired-end FASTQ files from RSV-positive samples.

2) Quality control and trimming
   - FastQC for initial QC.
   - Trimmomatic for adapter and quality trimming.

3) Reference-based assembly and masked consensus
   - Fetch RSV-A NC_038235.1 or RSV-B NC_001781.1 by accession.
   - Index references with BWA, samtools, Picard.
   - Map reads with BWA-MEM, sort, mark duplicates.
   - Compute depth and mask low coverage positions to N using a user-set threshold default 10.
   - Call variants with bcftools and generate a masked consensus with bcftools consensus.

4) De novo assembly and contig validation
   - Assemble with SPAdes.
   - Contig QC: length, percent Ns, GC, N50. Expect a principal contig near 15.2 kb.
   - Identify closest matches using BLAST remote against NCBI nt.
   - Fetch selected GenBank sequences for context.

5) Multiple sequence alignment and phylogeny
   - MAFFT alignment of your consensus plus context sequences.
   - IQ-TREE with configurable model. Use GTR+G for strict parity with the paper or MFP to enable ModelFinder automatic model selection.
   - 1000 ultrafast bootstraps by default. SH-aLRT can be enabled if required.

6) Clade and genotype assignment
   - Nextclade datasets for RSV-A and RSV-B to assign genotypes and clades.

7) Mutational analysis
   - Amino acid differences in G and F compared to RSV-A NC_038235.1 or RSV-B NC_001781.1. Results exported as TSV.

8) Outputs
   - Per sample consensus: results/consensus/<sample>.fa
   - Combined consensus: results/consensus/all_consensus.fasta
   - Alignment: results/aln/consensus_alignment.fasta
   - Maximum likelihood tree: results/iqtree/consensus.treefile
   - QC tables and mutation reports under results/

9) Repro and compliance
   - CITATION.cff, MIT LICENSE, CI check for quick verification.
   - Never commit restricted data. Set NCBI_EMAIL once for Entrez-based steps.

## Requirements
- Python 3.11 or newer
- Option A: pip and virtualenv
- Option B: conda or mamba
- Snakemake if you want to run the full pipeline
- Optional: Docker

### NCBI usage note
Set a contact email once per shell for E-utilities steps. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

## Quick verification
This confirms the repo executes locally and produces an image.
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -r env/requirements.txt
python analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
```

## One-command end to end run
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate rsv-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit config/config.yaml. Example:
```yaml
pairs:
  - sample: RSV_A_001
    subtype: A
    r1: data-private/RSV_A_001_R1.fastq.gz
    r2: data-private/RSV_A_001_R2.fastq.gz
  - sample: RSV_B_001
    subtype: B
    r1: data-private/RSV_B_001_R1.fastq.gz
    r2: data-private/RSV_B_001_R2.fastq.gz

ref_accessions:
  A: NC_038235.1   # RSV-A reference
  B: NC_001781.1   # RSV-B reference

tree_refs_acc:
  - OP890340.1
  - OR143212.1
  - OQ024128.1

params:
  threads: 4
  trim_adapters: env/TruSeq3-PE.fa
  min_depth_consensus: 10
  min_qual: 20
  iqtree_model: GTR+G    # or MFP for ModelFinder auto selection
  bootstrap: 1000
```

You can override parameters at runtime if needed:
```bash
snakemake -s workflow/Snakefile -c 4 --config params:iqtree_model=MFP params:bootstrap=1000
```

## Reference-based assembly and consensus details
Per sample:
1) Fetch reference by subtype and index it.
2) Trim adapters with Trimmomatic.
3) bwa mem mapping, samtools sort.
4) Picard MarkDuplicates then samtools index.
5) Depth and mask bed: positions with depth below min_depth_consensus are masked to N.
6) bcftools mpileup and bcftools call, then bcftools consensus with the mask to produce a high-quality consensus.

Primary deliverable: results/consensus/<sample>.fa

## De novo assembly and contig QC
**SPAdes and contig QC**
```bash
spades.py -1 results/<sample>_R1.trim.fastq.gz -2 results/<sample>_R2.trim.fastq.gz -o results/spades/<sample>
python analysis/scripts/contig_qc.py --in results/spades/<sample>/contigs.fasta --out_tsv results/spades/<sample>/contigs.qc.tsv
# Inspect the .summary.txt. Expect one principal contig near 15200 bases.
```

**Closest matches and context set**
```bash
python analysis/scripts/blast_closest.py --in results/spades/<sample>/contigs.fasta --out_tsv results/blast/<sample>.nt.tsv --max_hits 100
cut -f2 results/blast/<sample>.nt.tsv | head -n 50 > config/refs.txt
python analysis/scripts/fetch_genbank.py --acc config/refs.txt --out_fasta results/refs.fasta
```

**Alignment and phylogeny**
```bash
mafft --auto results/refs.fasta > results/aln/context.fasta
iqtree -s results/aln/context.fasta -m GTR+G -bb 1000 -nt 4 -pre results/iqtree/context
```

## Clade and genotype assignment
```bash
nextclade dataset get --name rsv-a --output-dir datasets/rsv-a
nextclade dataset get --name rsv-b --output-dir datasets/rsv-b

# Run on your consensus
nextclade run --input-dataset datasets/rsv-a --output-all results/nextclade/a results/consensus_rsv_a.fasta
nextclade run --input-dataset datasets/rsv-b --output-all results/nextclade/b results/consensus_rsv_b.fasta
```

## Mutational analysis in G and F
Report amino acid differences vs the appropriate reference.
```bash
# RSV-A
python analysis/scripts/mutations.py --subtype A --genomes results/consensus_rsv_a.fasta --out_tsv results/mutations_rsv_a.tsv
# RSV-B
python analysis/scripts/mutations.py --subtype B --genomes results/consensus_rsv_b.fasta --out_tsv results/mutations_rsv_b.tsv
```
These scripts use Biopython Entrez to fetch reference annotations. Set NCBI_EMAIL once per shell or pass --email.

## Outputs checklist
- results/consensus/<sample>.fa and results/consensus/all_consensus.fasta
- results/aln/consensus_alignment.fasta
- results/iqtree/consensus.treefile and associated IQ-TREE files
- results/spades/<sample>/contigs.qc.tsv and .summary.txt if you run the de novo path
- results/mutations_rsv_*.tsv for G and F differences

## Project layout
```
analysis/           scripts
config/             configuration
env/                environment files
workflow/           Snakemake workflow
data-example/       example dataset for quick verification
results-example/    example outputs after quick verification
docs/               SOP and notes
.github/workflows/  CI
```

## Data access and compliance
- Do not upload raw reads or restricted data.
- GISAID terms apply for sequences obtained through GISAID.
- For any E-utilities requests, set NCBI_EMAIL as shown earlier. Optional NCBI_API_KEY can improve the rate limit.

## Methods and tool references
- FastQC v0.11.9, Trimmomatic v0.39, Picard MarkDuplicates, SPAdes v3.15.5, BWA, MAFFT, IQ-TREE with GTR+G and 1000 bootstraps, Nextclade CLI.
- SPAdes: Bankevich A et al., J Comput Biol 2012, 19(5):455-477.
- Geneious: Kearse M et al., Bioinformatics 2012, 28(12):1647-1649.
- RSV genotype definition: Goya S et al., Influenza Other Respir Viruses 2020, 14(3):274-285.
- IQ-TREE: Nguyen LT et al., Mol Biol Evol 2015, 32(1):268-274.

## How to cite
- Cite the paper using the DOI above.
- Cite this repository. CITATION.cff is included so GitHub can display the cite box.
