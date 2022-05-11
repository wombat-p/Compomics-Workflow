#!/bin/bash -ue
echo "raw_file	exp_condition	biorep
OVEMB150205_14	Main	A1
OVEMB150205_12	Main	A2" > none
  cp "none" exp_design.txt
  mv "QuantifiedPeptides.tsv" q_input.txt
  mv "QuantifiedProteins.tsv" q_prot.txt
  Rscript /home/veits/devel/Bioinformatics/ELIXIR_EDAM/WOMBAT-P/Compomics-Workflow/Nextflow/scripts/runMSqRob.R
