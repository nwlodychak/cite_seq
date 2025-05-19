git clone https://github.com/xjtu-omics/msisensor-rna.git
pip3 install .
#msisensor-rna -i <input> -o <output>
python MSIsensor-RNA.py predict -i your_expression_data.txt -o output_results


#-h, --help            show this help message and exit
#-i INPUT, --input INPUT
#                    The path of input file. e.g. xxx.csv [required]
#-o OUTPUT, --output OUTPUT
#                    The output file of gene information. e.g. xxx.csv [required]
#-thresh_t THREADS, --threads THREADS
#                    The threads used to run this program. [default=4]
#-thresh_cov THRESH_COV, --thresh_cov THRESH_COV
#                    Threshold for coefficient of variation of gene expression value of all samples (Mean/Std). [default=0.5]
#-thresh_p THRESH_P_RANKSUM, --thresh_p_ranksum THRESH_P_RANKSUM
#                    Threshold for Pvalue of rank sum test between MSI-H and MSS samples. [default=0.01]
#-thresh_auc THRESH_AUCSCORE, --thresh_AUCscore THRESH_AUCSCORE
#                    Threshold for AUC score: AUC score was calculating by the sklearn package. [default=0.65]
#-p POSITIVE_NUM, --positive_num POSITIVE_NUM
#                    The minimum positive sample of MSI for training. [default = 10]