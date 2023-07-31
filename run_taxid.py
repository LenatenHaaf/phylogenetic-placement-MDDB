import subprocess
import os
import time
start_time = time.time()
'''
Script to run taxit create
Change paths to correct paths
'''
def run_taxit_create(input_file, package_name, aln_fasta, tree_stats, tree_file):
    command = [
        'taxit',
        'create',
        '-l', 'its',
        '-P', package_name + '.refpkg',
        '--aln-fasta', aln_fasta,
        '--tree-stats', tree_stats,
        '--tree-file', tree_file
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode == 0:
        print(f"Taxit create succeeded for {input_file}")
    else:
        print(f"Taxit create failed for {input_file}. Error: {stderr.decode('utf-8')}")

input_directory = 'Test3/removed_fasta'

for filename in os.listdir(input_directory):
    number = filename[8:11]
    package_name = 'Test3/packages_test3_timetest/timetest_' + number
    if filename.endswith('.fasta'):
        input_file = os.path.join(input_directory, filename)
        aln_fasta = input_file
        tree_stats = 'Test3/raxml_test3/test3/RAxML_info.' + filename[:11] + '.out'
        tree_file = 'Test3/raxml_test3/test3/RAxML_bestTree.' + filename[:11] + '.out'

        run_taxit_create(input_file, package_name, aln_fasta, tree_stats, tree_file)

end_time = time.time()
run_time = end_time - start_time
print(run_time)