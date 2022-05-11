#!/bin/bash -ue
# needed for Myrimatch, see https://github.com/compomics/searchgui/issues/245
LANG=/usr/lib/locale/en_US
export LC_ALL=C; unset LANGUAGE
 mkdir tmp
 mkdir log
 searchgui eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ./tmp -log ./log
 searchgui eu.isas.searchgui.cmd.SearchCLI -spectrum_files ./  -output_folder ./ -fasta_file "./ABRF_iPRG_2012_target_concatenated_target_decoy.fasta"  -id_params "./searchgui.par" -threads 2 \
     -xtandem 0 -msgf 0 -comet 0 -ms_amanda 0 -myrimatch 1
 mv searchgui_out.zip OVEMB150205_12.zip
