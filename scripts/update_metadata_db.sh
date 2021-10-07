#!/bin/bash

CONFIG=$1

bs list datasets -c $CONFIG \
  -F Name \
  -F Project.Name \
  -F Id \
  -F Project.Id \
  -F TotalSize \
  -F DataSetType.Id \
  -F DateCreated \
  -F DateModified
