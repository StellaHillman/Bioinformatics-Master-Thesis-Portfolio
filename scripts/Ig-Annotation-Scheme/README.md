Ig-Annotation-Sheme
======================

This directory contains files used to make the Ig numbering. The initial annotations were made using the refund_files.js in the "IgStrad-numbering-icn3d" directory as follows: 
for file in /path/directory-with-pdb-structures/*.pdb; do node refnum_file.js "$file" pdb; done

Then the output json files (not included here) were converted into csv and corrected in excel: numbered-FnIII-ttn.csv and numbered-Ig-ttn.csv. 

Then they were renumbered according to my own modfied numbering system using: "New-Numbering.R" -> FnIII-numbering_new-annotations.csv and Ig-numbering_new-annotations.csv

The excel contains spreadsheets where the numbering was visualised and explored. 