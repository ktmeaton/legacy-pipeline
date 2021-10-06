#!/bin/bash

bs list projects \
  -F Name \
  -F Id \
  -F Project.TotalSize \
  -F DateCreated \
  -F DateModified
