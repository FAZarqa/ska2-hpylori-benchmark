# H. pylori SKA2 transmission genomics — analysis scripts

Custom scripts supporting the evaluation of SKA2 as a reference-free
distance estimator for *Helicobacter pylori* transmission epidemiology,
benchmarked against DNAdiff SNP distances and a recombination-masked
IQ-TREE phylogeny.

## Pipeline overview

1. Genome QC and representative genome selection (BUSCO, FastANI, QUAST, CheckM2)
2. Panaroo core-gene alignment, re-aligned per gene with MAFFT
3. fastGEAR recombination detection on each gene alignment
4. Masking of recombinant regions and concatenation into a core-genome alignment
5. IQ-TREE ML phylogeny on the masked alignment
6. InSilicoSeq read simulation (30x/100x) and Shovill/SKESA reassembly
7. SKA2 distance estimation (k = 31, 17, 15) on complete genomes, simulated
   reads, and reassemblies
8. Comparison of SKA2 distances against DNAdiff and the masked ML phylogeny

## Scripts

| Script | Purpose |
|---|---|
| `concat_masked_alignments.py` | Concatenates per-gene alignments (post fastGEAR masking) into one core-genome alignment, gap-filling any genome missing from a given gene |
| `fastgear_import.sh` | Converts a fastGEAR output directory into a single importation status file (recent + ancestral events) |
| `mask_recombinant_regions.py` | Masks recombinant regions in a gene alignment with N, using the importation status file and a minimum length cutoff |
| `plot_ska2_vs_benchmark.R` | Generic function to scatter-plot SKA2 distances against a benchmark distance matrix (DNAdiff or ML), with Pearson r, R², and n annotated |


## Requirements

Python 3 with Biopython (for the masking script). R with `ggplot2`,
`dplyr`, `ape`, `phangorn`, `vegan`, and `patchwork` for the plotting
scripts.

## Usage

```bash
# generate importation status file from a fastGEAR output directory
bash fastgear_import.sh path/to/fastgear/output/gene1.aln/output

# mask recombinant regions >= 100 bp
python mask_recombinant_regions.py --align gene1.aln.fas \
    --imp gene1.aln/output/FastGear.importation_status.txt --min 100

# concatenate all masked gene alignments into one core-genome alignment
python concat_masked_alignments.py --alignments masked_alignments/ \
    --genomes genome_list.txt
```

```r
# compare SKA2 distances against DNAdiff or ML distances
source("plot_ska2_vs_benchmark.R")
plot_ska2_vs_benchmark(
  ska_file = "SKA_distances/hpylori_distance_k31.txt",
  benchmark_file = "dnadiff_cleaned.txt",
  benchmark_col = "SNP_distance",
  k_label = "k31",
  benchmark_label = "DNAdiff",
  out_dir = "figures"
)
```

## Data availability

Genome accessions are listed in Supplementary Table S1. 
