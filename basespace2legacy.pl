#!/bin/bash

INDIR=$1

for old_sample_dir in `ls -d $INDIR/*`; do
  sample_id=`basename $old_sample_dir | cut -d "-" -f 1 `;
  new_sample_dir="${INDIR}Sample_${sample_id}";
  echo -e "Renaming:\n\t$old_sample_dir to\n\t$new_sample_dir";
  mv $old_sample_dir $new_sample_dir;  
done
