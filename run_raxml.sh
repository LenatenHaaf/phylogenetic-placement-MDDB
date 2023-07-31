input_dir="/home/s2566818/standard-RAxML-master/removed_fasta_test3"
output_dir="/home/s2566818/standard-RAxML-master/test3"

for input_file in "$input_dir"/*.fasta; do
    file_name=$(basename "${input_file%.*}")
    ./raxmlHPC-PTHREADS-SSE3 -s "$input_file" -n "$file_name.out" -w "$output_dir" -m GTRGAMMA -p 12345 -T 8
done