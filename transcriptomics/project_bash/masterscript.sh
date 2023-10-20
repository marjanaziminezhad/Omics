#!/bin/bash
set -x
set -e

# create Dir Structure, and links to raw bam files
bash createDirStructure_symLinkToRawBam.sh
# create links to run_toolX.sh
bash makeSymLinkToAllRunScripts.sh
# script that calls one after each other the for_loop_toolX.sh scripts
bash loop_script.sh

#  all bash files shoul be in the same dir
