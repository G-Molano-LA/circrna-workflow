

pushd data/raw_data

echo "Converting GTF annotationg file to refFlat GenePred table format.. "

gtfToGenePred -genePredExt  GRCh38_ann.gtf refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > GRCh38_refFlat.txt
rm refFlat.tmp.txt

echo "... DONE."

popd
