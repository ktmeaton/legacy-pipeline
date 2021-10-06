#!/bin/bash

INDIR=$1
OUTDIR=$2
ID_SEP=$3
ID_SEP="${ID_SEP:=_}"

for old_sample_dir in `ls -d $INDIR/*`; do

  # Check if directory
  if [[ ! -d $old_sample_dir ]]; then continue; fi;

  sample_id=`basename $old_sample_dir | cut -d "$ID_SEP" -f 1 `;
  new_sample_dir="${OUTDIR}/Sample_${sample_id}";

  # Create new directory
  echo -e "\nCreating sample directory: $new_sample_dir";
  mkdir -p $new_sample_dir;

  # Symlink files
  for file in `ls $old_sample_dir/*.fastq.gz`; do
    filename=`basename $file`;
    target="${new_sample_dir}/${filename}";
    echo -e "\t$target";
    ln -s $file $target;
  done;

  break

done
