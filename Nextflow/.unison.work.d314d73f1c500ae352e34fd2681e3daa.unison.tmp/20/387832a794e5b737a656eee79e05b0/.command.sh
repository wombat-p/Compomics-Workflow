#!/bin/bash -ue
peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
 peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in "./OVEMB150205_12.psdb" -out_reports "./" -reports "3,4,6,9"
mv "false_Default_PSM_Report_with_non-validated_matches.txt" "OVEMB150205_12.txt"
mv "false_Default_PSM_Report.txt" "OVEMB150205_12_filtered.txt"
mv "false_Default_Peptide_Report.txt" "OVEMB150205_12_peptides.txt"
mv "false_Default_Protein_Report.txt" "OVEMB150205_12_proteins.txt"
