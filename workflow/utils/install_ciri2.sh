#!/bin/bash

################################################################################
## Script to manage CIRI2 dependencies
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: G.Molano, LA (gonmola@hotmail.es)
################################################################################

# 2. Installing CIRI2
wget https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip
unzip CIRI_v2.0.6.zip
rm -r CIRI_v2.0.6.zip
mv CIRI_v2.0.6/CIRI2.pl tools/
rm -r CIRI_v2.0.6/
rm -r __MACOSX/
