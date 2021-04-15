#!/bin/bash

################################################################################
## Dependencies Installation on Linux machines
################################################################################

# 2. Installing CIRI2
wget https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip
unzip CIRI_v2.0.6.zip
rm -r CIRI_v2.0.6.zip
mv CIRI_v2.0.6/CIRI2.pl libs/ciri2/
rm -r CIRI_v2.0.6/
rm -r __MACOSX/