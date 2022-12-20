## Download data from:
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

## get reference
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf

mkdir reference
mv chr22* reference/

## Decompress:
tar -xvf HBR_UHR_ERCC_ds_5pc.tar 

# Order folder
rm HBR_UHR_ERCC_ds_5pc.tar 
mkdir reads
mv *gz reads/

## create new_names
for i in `dir -l reads/ | awk '{print $NF}'`; do ori_name=$i; new_name=`echo $i | sed "s/_ERCC//g" | sed "s/UHR_//g" | sed "s/HBR_//" | sed "s/_Build37-ErccTranscripts-chr22//" | sed "s/.read/_R/"`; echo "mv reads/$i reads/$new_name"; done | grep 'gz' > rename.sh
sh rename.sh
