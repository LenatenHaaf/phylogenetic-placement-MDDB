# Script to run guppy tog on the output of pplacer
for file in /home/s2566818/output_test3_bay/*.jplace; do
    guppy tog -o "${file%.jplace}.tog.tre" "$file"
done