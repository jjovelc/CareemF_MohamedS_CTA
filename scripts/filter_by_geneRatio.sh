#!/usr/bin/bash

INFILE=$1

awk -F'\t' 'NR == 1 || ($3 ~ "/" && split($3, ratio, "/") && (ratio[1] / ratio[2] >= 0.02))' "$INFILE" > "${INFILE/.tsv/_filtered.tsv}"
