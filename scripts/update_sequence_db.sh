#!/bin/bash

CONFIG=$1
CONFIG="${CONFIG:=secondary}"
max_projects=$2
max_projects="${max_projects:=100}"
args=("$@");

if [[ ${#args[@]} < 1 ]]; then
    echo ""
    echo -e "Description : Update the Illumina Sequencing Database from BaseSpace."
    echo -e "Version     : v0.1.0";
    echo -e "Date        : 2021-10-06";
    echo -e "Usage       : update_sequence_db.sh (config) (max_num_projects)";
    echo -e "Notes       : max_num_projects\tINT [1]";
    echo -e "            : config STR [secondary]"
    echo ""
    exit;
fi

project_names=(`
  bs list projects -c $CONFIG -F Name | \
  grep -v "\-\+\-" | \
  tail -n+2 | \
  sed 's/|//g' | \
  awk '{$1=$1;print}' | \
  tail -n $max_projects | \
  sed 's/ /_/g' | \
  tac`); 

project_ids=(`
  bs list projects -c $CONFIG -F Id | \
  grep -v "\-\+\-" | \
  tail -n+2 | \
  sed 's/|//g' | \
  awk '{$1=$1;print}' | \
  tail -n $max_projects | \
  tac`); 

for i in "${!project_names[@]}"; do
  # TBD: Catch spaces!!
  name=${project_names[$i]};
  id=${project_ids[$i]};
  
  progress=`expr $i + 1`;
  echo -e "\nDownloading project ($progress/$max_projects): $name\t$id";
  bs download project -c $CONFIG --quiet -i $id -o $name;
done
