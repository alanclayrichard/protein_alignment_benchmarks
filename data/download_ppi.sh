#!/bin/sh
# download the human ppi gold-standard test set (figshare) + human swissprot sequences.
# the Intra*_rr.txt files list interacting/non-interacting uniprot accession pairs;
# human_swissprot_oneliner.fasta holds the sequences keyed by accession.
set -e
d="$(dirname "$0")/ppi"
mkdir -p "$d"
for pair in Intra0_neg:41270466 Intra0_pos:41270469 Intra1_neg:41270472 \
            Intra1_pos:41270475 Intra2_neg:41270478 Intra2_pos:41270481; do
  curl -fsSL "https://ndownloader.figshare.com/files/${pair##*:}" -o "$d/${pair%%:*}_rr.txt"
done
curl -fsSL https://ndownloader.figshare.com/files/52237748 -o "$d/human_swissprot_oneliner.fasta"
