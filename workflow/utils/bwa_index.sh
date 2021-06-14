
genome=$1
path=$2

if [[ $genome == "GRCh38" ]];
then
  echo "Downloading genome index files from NCBI repository (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz)..."
  wget -c  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz \
      -O "$genome/bwa/GRCh38.bwa_index.tar.gz"
  tar -zxvf "$genome/bwa/GRCh38.bwa_index.tar.gz" -C "$genome/bwa"

  pushd "$path/bwa/"
  rm GRCh38.bwa_index.tar.gz
  rename 's/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna/GRCh38/' *
  popd
elif [[ $genome == "GRCh37" ]];
then
  echo "Creating a genome index file from GRCh37.fna genome"
  bwa index -a bwtsw -p "$genome" "$path/$genome.fna"
else
  echo "Please specify a genome in the config.yaml file. Genome: hg38/GRCh38 or hg37/GRCh37/hg19."
fi
