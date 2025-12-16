echo "======================Bash output======================"

# Write line to put into the output file 
echo "a) List fasta files"
# command line
ls examples/*/seq_*_chain*.fasta | sort > 1_output_list.txt
# sorting file names alphabetically
# test output with 
cat 1_output_list.txt

echo "====================================================="

# Write sentence to put into the output file
echo "b) Count total number of pattern-matching fasta files"
# execution line
ls examples/*/seq_*_chain*.fasta | wc -l > 1_output_count.txt
# counting word count within each line “-l”
cat 1_output_count.txt

echo "====================================================="

# Write line to put into the output file
echo "c) Calculate total sequence length of all fasta files"
# execution line
sed '/^>/d' examples/*/seq_*_chain*.fasta | wc -m > 1_output_length.txt
# deleting header lines starting with “>”, counting words, sums the length of remaining lines
cat 1_output_length.txt
