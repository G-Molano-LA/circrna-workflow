
# http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format#The_opposite_direction.2C_GTF_to_GenePred
#

path=$1
pushd "$path/"

echo "Converting GTF annotationg file to refFlat GenePred table format.. "

gtfToGenePred -genePredExt  -geneNameAsName2 GRCh37_ann.gtf refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > GRCh37_refFlat.txt # refFlat Format

# awk -F "\t" '$3="chr"$3' GRCh37_refFlat.txt > GRCh37_refFlat_UCSC.txt # changing chr format

rm refFlat.tmp.txt
# rm GRCh37_refFlat.txt

echo "... DONE."

popd
