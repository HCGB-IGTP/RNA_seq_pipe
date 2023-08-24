echo "`date` ... Starting ..."
echo ""

mkdir logs

##########################################
## Test toy dataset
##########################################
echo "# ------------------------------ #"
echo "RSP prep -i toy_dataset/reads/ -o RSP_analysis"
echo "..."
echo ""
RSP prep -i toy_dataset/reads/ -o RSP_analysis | tee logs/RSP_analysis.prep.log
echo ""

echo "# ------------------------------ #"
echo "RSP QC -i RSP_analysis/"
echo "..."
echo ""
RSP QC -i RSP_analysis/ | tee logs/RSP_analysis.qc.log
echo ""

echo "# ------------------------------ #"
echo "RSP trim -s trimmomatic -i RSP_analysis/"
echo "..."
echo ""
RSP trim -s trimmomatic -i RSP_analysis/ | tee logs/RSP_analysis.trim.log
echo ""

echo "# ------------------------------ #"
echo "RSP map -s star -i RSP_analysis/ --ref_genome toy_dataset/reference/chr22_with_ERCC92.fa --ref_folder index_folder --ref_name chr22"
echo "..."
echo ""
RSP map -s star -i RSP_analysis/ --ref_genome toy_dataset/reference/chr22_with_ERCC92.fa --ref_folder index_folder --ref_name chr22 | tee logs/RSP_analysis.map.log
echo ""

echo "# ------------------------------ #"
echo "RSP count -s star -i RSP_analysis/ --ref_annot toy_dataset/reference/chr22_with_ERCC92.gtf"
echo "..."
echo ""
RSP count -s star -i RSP_analysis/ --ref_annot toy_dataset/reference/chr22_with_ERCC92.gtf | tee logs/RSP_analysis.counts.log
echo ""
