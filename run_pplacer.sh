# Script to run pplacer for the first test
# #!/bin/bash

# # Path to the directory containing reference packages
# refpkg_dir="/home/s2566818/packages"

# # Path to the directory containing query sequences
# query_dir="/home/s2566818/queries1"

# output_dir="/home/s2566818/output_test1_bay" 

# # Loop over each query file in the query directory
# for query_file in "$query_dir"/*.fasta; do
#     # Get the query file name without extension
#     query_name=$(basename "$query_file" .fasta)

#     # Build the path to the corresponding reference package
#     refpkg_file="$refpkg_dir/${query_name}.refpkg"

#     # Run pplacer with the query file and reference package
#     pplacer -c "$refpkg_file" "$query_file" -o "$output_dir/$query_name.jplace" -p
# done

# Script to run pplacer for the second test

#!/bin/bash

# # Path to the directory containing reference packages
# refpkg_dir="/home/s2566818/packages_test2_best"

# # Path to the directory containing query sequences
# query_dir="/home/s2566818/queries_test2"

# output_dir="/home/s2566818/output_test2_bay" 

# # Loop over each query file in the query directory
# for query_file in "$query_dir"/*.fasta; do
#     # Get the query file name without extension
#     query_name=$(basename "$query_file" .fasta)

#     # Build the path to the corresponding reference package
#     refpkg_file="$refpkg_dir/${query_name}.refpkg"

#     # Run pplacer with the query file and reference package
#     pplacer -c "$refpkg_file" "$query_file" -o "$output_dir/$query_name.jplace" -p
# done

# Script to run pplacer for the third test

# Path to the directory containing reference packages
refpkg_dir="/home/s2566818/packages_test3_best"

# Path to the directory containing query sequences
query_dir="/home/s2566818/queries"

output_dir="/home/s2566818/output_test3_bay" 

# Loop over each query file in the query directory
for query_file in "$query_dir"/*.fasta; do
    # Get the query file name without extension
    query_name=$(basename "$query_file" .fasta)

    # Build the path to the corresponding reference package
    refpkg_file="$refpkg_dir/${query_name}.refpkg"

    # Run pplacer with the query file and reference package
    pplacer -c "$refpkg_file" "$query_file" -o "$output_dir/$query_name.jplace" -p
done