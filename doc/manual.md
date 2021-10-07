---
project: [[legacy-pipeline]]
people:
  - [[Ana Duggan]]
  - [[Gabriel Renaud]]
  - [[Brian Golding]]
  - [[Hendrik Poinar]]
  - [[Katherine Eaton]]
  - [[Jessica Hider]]
authors:
  - name: [[Ana Duggan]]
    orcid: 
    affiliations:
      [	
        "[[McMaster Ancient DNA Center]]",
        "[[Department of Anthropology]], [[McMaster University]]"
      ]
lang: en-US
repo: ktmeaton/legacy-pipelinee
tags: 🧨
status: priority
title: Legacy Pipeline
subtitle: The legacy pipeline for the McMaster Ancient DNA Centre.
type:
  - [[Project]]
  - [[Tool]]
  - [[Note]]
numberSections: False
autoSectionLabels: True
sectionsDepth: 3
tblPrefix: Table
figPrefix: Figure
secPrefix: Section
---


## Introduction

The Legacy Pipeline uses Make to pre-processes raw data from FASTQ inputs, align reads to a reference genome, and create input fasta files for metagenomic analysis.

## Pipeline Steps


- Create reference genome indices for mapping (bwa, samtools)
- Organize and rename input FASTQ sequences.
- Format conversion (fastq -> bam)
- Sequencing adapter removal and for paired end data merging (leeHom)
- Read mapping to reference using (bwa aln)
- Post-mapping filtering (samtools)
- PCR duplicate removal (bam-rmdup, retrieveMappedProperlyPaired...)
- Ancient DNA C-to-T damage pattern visualisation (mapDamage)
- Metagenomics filtering (awk)

## Preparation

Remember that the purpose of keeping all the references in a single location is first and foremost consistency across users and projects.

1. Upload your reference.fasta to a communal directory.
    - Example: `/home/username/Reference_sequences`
2. Index the reference for pipeline compatibility\*.
    - `bwa index reference.fasta`
3. Index the reference for mapDamage compatibility
    - `samtools faidx reference.fasta`
4. Make sure the read and execute permissions are open to all
    - `chmod +rx reference.fasta*`

\* network-aware-bwa/bwa

## Running The Pipeline

Updated pipeline usage:
- As per November 26, 2019 lab meeting discussion, this updated pipeline produces new files meant to ease data entry into metagenomic classifiers such as BLAST or Kraken, creates a text file summarising the exhaustion estimate for each library , and automatically adds a map quality filter in the final step.

Remember to save your .fastq.gz files pre-organised into folders such as `Sample_ID`.

    ```bash
    for f in `ls *R1_001.*fastq.gz`; 
    do 
        mkdir Sample_${f%%_*} | mv ${f%%_*}* ./Sample_${f%%_*}; 
    done
    ```

Use legacy.pl to generate a Makefile which contains all the commands necessary to process your data.

    ```bash
    legacy.pl /path/to/data/folder /home/poinarlab Reference_sequences/reference.fasta > Makefile
    ```
Launch make to process data\*.
- `make -j N`

\* Number of CPUs/samples to process in parallel.
