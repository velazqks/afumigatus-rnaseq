#!/bin/bash

# Directories
fungidb="https://fungidb.org/common/downloads/release-2.0/Afumigatus_AF293B"
data_dir="../data"
ref_genome="$data_dir/Afumigatus_AF293B_Genome.fasta"
annotation="$data_dir/a_fumigatus_Af293"
index_prefix="$data_dir/af_index"
# Create a list of sample names and gtf file directories. Used to generate DE counts
sample_gtf_list="$data_dir/de_samples.txt"
# Clear the file if it exists
> "$sample_gtf_list"

# Create data directory if it doesn't exist
#mkdir -p "$data_dir"

# Function to download reference files
download_files() {
    # Download reference genome and annotation
    wget -q "$fungidb/fasta/data/FungiDB-CURRENT_Afumigatus_AF293B_Genome.fasta" -O "$ref_genome"
    wget -q "$fungidb/gff/data/a_fumigatus%20Af293.gff" -O "$annotation.0.gff"
    
    # Format ids similar to the ids in the reference genome
    sed 's/apidb/JVCI/g' "$annotation.0.gff" > "$annotation.gff"
    
    # Convert file type from gff to gtf
    agat_convert_sp_gff2gtf.pl --gff "$annotation.gff" -O "$annotation.gtf"
    
    # Index the reference genome
    hisat2-build -p 8 "$ref_genome" "$index_prefix"
}

# Function to process SRA data
process_sra() {
    sra="$1"
    sra_dir="$data_dir/$sra"
    sra_sra="$sra_dir/$sra.sra"
    
    # Download SRA data
    #prefetch "$sra" --output-file "$sra_sra"
    
    # Convert SRA to FASTQ
    #fastq-dump "$sra_sra" --outdir "$sra_dir"
    #gzip "$sra_dir/$sra.fastq"
    
    # Map reads to reference genome
    hisat2 -p 8 --dta -x "$index_prefix" -U "$sra_dir/$sra.fastq.gz" -S "$sra_dir/$sra.sam"
    
    # Convert SAM to BAM and sort
    samtools sort -@ 8 "$sra_dir/$sra.sam" -o "$sra_dir/$sra.bam"
    samtools index "$sra_dir/$sra.bam"
    
    # Assemble and quantify transcripts
    stringtie -p 8 "$sra_dir/$sra.bam" -G "$annotation.gtf" -o "$sra_dir/$sra.gtf"
    
    # Append sample name and GTF file directory to the list file
        echo "$sra $sra_dir.gtf" >> "$sample_gtf_list"

    
    # Clean up
    #rm "$sra_sra"
    #rm "$sra_dir.sam"
}

# Function to estimate transcript abundances and create table counts for Ballgown
transcript_abundances() {
    sra="$1"
    sra_dir="$data_dir/$sra"
    stringtie -e -B -p 8 "$sra_dir/$sra.bam" -G "$data_dir/stringtie_merged.gtf" -o "$data_dir/ballgown/$sra/$sra.gtf"
}


# Download genome reference files
download_files

# Process SRA data
process_sra "SRX2000912"
process_sra "SRX2000913"
process_sra "SRX2000916"
process_sra "SRX2000917"

# Merge transcripts from all samples
stringtie --merge -p 8 -G "$annotation.gtf" -o "$data_dir/stringtie_merged.gtf" "$data_dir/mergelist.txt"

# Compare the transcripts with the reference annotation
gffcompare -r "$annotation.gtf" -G "$data_dir/stringtie_merged.gtf" -o "$data_dir/merged"

# Estimate transcript abundances and create table counts for Ballgown
transcript_abundances "SRX2000912"
transcript_abundances "SRX2000913"
transcript_abundances "SRX2000916"
transcript_abundances "SRX2000917"


# Generate transcript counts from 1) a list of gtf samples and paths, or 2) a ballgown path
prepDE.py -i "$data_dir/ballgown/" -t ""$data_dir/transcript_count_matrix.csv"
