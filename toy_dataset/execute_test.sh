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