diamond blastp --query rawdata/$1/$1.prokka/$1.faa -d all_seqs.dmnd -t dmnd_tmp -o rawdata/$1/$1.prokka/$1.m8 -f tab --max-target-seqs 10 --min-score 50 --threads 6
