#!/bin/bash

bs list datasets \
  -F Name \
  -F Project.Name \
  -F Id \
  -F Project.Id \
  -F TotalSize \
  -F DataSetType.Id \
  -F DateCreated \
  -F DateModified
