FUNGI RNA-SEQ DATA PROCESSING AND ANALYSIS
==========================================
` velazqks `

Fungi have diverse and unique lifestyles. For instance, the digestion of their food takes place outside of the fungal body before consuming it. This mechanism is aided by the secretion of enzymes and acids into the environment ([Brakhage, 2013](https://doi.org/10.1038/nrmicro2916), [Keller, 2019](https://doi.org/10.1038/s41579-018-0121-1)). The ability to synthesize a variety of molecules that support the survival in hostile environments has probably been one of the bases for the evolution of some of the most pathogenic organisms we know. Fungi in the Aspergillus genus are especially recognized as the most pathogenic, and A. fumigatus is the most pathogenic of these species, responsible for 90% of the infections ([Lin, 2001](https://doi.org/10.1086/318483)). Living inside a host organism poses many challenges to pathogens. To evade the immune system, reactive oxidative species (ROS) and heat shock stress, mechanisms to surpass these stress responses are essential ([Grützmann et al., 2014](https://doi.org/10.1093/dnares/dst038)). Many genes involved in pathogenicity in A. fumigatus have been reported ([Abad et al., 2010](https://doi.org/10.1016/j.riam.2010.10.003)). The most relevant groups are involved in thermotolerance, cell wall composition, immune response avoidance, toxin production, nutrient acquisition, and response to stress.

Alternative splicing generates higher variation in the expression in the genome through the transcription of various of isoforms from a single gene. This process has been associated with higher biological complexity and multicellularity. A research study suggested that alternative splicing not only is involved in enhancing the evolutionary dynamics in fungi, but also in regulating pathogenicity and virulence ([Grützmann et al., 2014](https://doi.org/10.1093/dnares/dst038)). 

Here, I processed RNA-Seq data from 4 wild-type A fumigatus samples ([Lind et al., 2016](https://doi.org/10.1534/g3.116.033084)): two samples grown at 30⁰ and two samples grown at 37⁰, as a proxy to body temperature. The reads were aligned to the FingiDB reference genome and annotations using HISAT2. The transcripts were assembled and quantified with StringTie. Transcript counts were generated using StringTie's prepDE.py. Transcript-level differential expression analysis was done using Ballgown.

`bin/pipeline_rnaseq.sh` generates the raw data.

`data/` already contains the Ballgown data to run differential expression analyses.

<br>

**NEXT STEPS**

**Perform gene-level analysis using the RNA-Seq data isoforms.**

1. Quantify transcript expression levels using StringTie to get FPKM or TPM values.
2. Aggregate the expression values transcripts per gene to get the gene-level expression values for each sample.
3. Perform differential gene expression (DE) analysis by comparing gene expression levels between different temperatures to identify genes that are differentially expressed. Use volcano plots, heatmaps, or bar plots to visualize DE (DESeq2 or edgeR).
4. Examine the biological processes and pathways that are affected by the DE genes: Functional enrichment analysis.


**Quantify alternative splicing (AS) using RNA-Seq data.**

1. Map the RNA-Seq reads to the reference genome (HISAT2 or STAR) and quantify transcript abundance (StringTie or featureCounts).
2. Identify AS events in each gene (e.g., exon inclusion or skipping, alternative splice-site selection, mutually exclusive exons, intron retention). Use rMATS, SUPPA, or Whippet.
3. For each AS event, calculate the exon-inclusion ratio, or percent spliced-in (PSI).
4. Perform differential AS analysis to identify splicing events that change between temperatures. Use rMATS, MISO, or DEXSeq.
5. Visualize the PSI values for each condition/sample. Create plots showing PSI values for different splicing events to highlight changes in AS.
6. Examine the biological consequences of the changes in AS: Functional enrichment analysis.


