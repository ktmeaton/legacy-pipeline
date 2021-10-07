#!/bin/bash

args=( "$@" );
OUTDIR=$1
INDIR=$2
ID_SEP="_"

# Check number of args
if [[ ${#args[@]} < 2 ]]; then
    echo ""
    echo -e "Description : Prepare a legacy project from a BaseSpace PROJECT."
    echo -e "Version     : v0.1.0";
    echo -e "Date        : 2021-10-07";
    echo -e "Usage       : prep_legacy_input_project.sh output_project input_project";
    echo ""
    exit;
fi

# Check if INDIR exists
if [[ ! -e $INDIR ]]; then
  echo "The input project directory '$INDIR' does not exist.";
  exit 1;
fi;



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
    if [[ ! -L $target ]]; then
      echo -e "\t$target";
      ln -s $file $target;
    fi;
  done;

done
