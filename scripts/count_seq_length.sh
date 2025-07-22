#!/bin/bash

if [ ! -f "$1" ]; then
    echo "Error: file '$1' not found."
    exit 1
fi



awk '/^>/{if (seq) print seq_id, "\t| Length:", length(seq); seq_id=$0; seq=""; next} \
    {seq = seq $0} END {if (seq) print seq_id, "\t| Length:", length(seq)}' "$1" | sed 's/>//'

