#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
from Bio import SeqIO
import pandas as pd
from pandas import DataFrame, Series
import seaborn as sns
import matplotlib.pyplot as plt
import math

def parse_args():
	parser = argparse.ArgumentParser(description='''
For parsing contig files and producing statistics on N50, coverage, etc. Expects a single
fasta-formatted file with a SPAdes format header.
	''')
	parser.add_argument('fastafile', type=str,help='relative path to fasta file containing contigs')
	parser.add_argument('len_cutoff', type=int,help='cutoff for length of contigs')
	parser.add_argument('cov_cutoff', type=float,help='cutoff for coverage of contigs')
	parser.add_argument('prefix',type=str,help='prefix for output files')
	return parser.parse_args()

def parse_seqs(fastafile,prefix,len_cutoff,cov_cutoff):
	o = open(prefix+".filteredcontigs.fna",'w')
	leng = {}
	cov = {}
	for seq in SeqIO.parse(open(fastafile,'r'),'fasta'):
		vals = seq.id.split("_")
		contigname = "{}_{}".format(vals[0],vals[1])
		if float(vals[5]) > 1000:
			continue
		else:
			leng[contigname] = math.log(int(vals[3]),10)
			cov[contigname] = float(vals[5])
		if leng[contigname] > math.log(len_cutoff,10) and cov[contigname] > cov_cutoff:
			seq.id = "{}_contig{}".format(os.path.basename(prefix),vals[1])
			seq.description = ""
			SeqIO.write(seq,o,'fasta')
	d = {'log10length':pd.Series(leng),'coverage':pd.Series(cov)}
	o.close()
	return pd.DataFrame(d)

def main():
	args = parse_args()
	fastafile = os.path.abspath(args.fastafile)
	len_cutoff = args.len_cutoff
	cov_cutoff = args.cov_cutoff
	prefix = os.path.abspath(args.prefix)

	df = parse_seqs(fastafile,prefix,len_cutoff,cov_cutoff)
	fig = plt.figure()
	df.plot(x="coverage",y="log10length",kind="scatter")
	plt.axhline(math.log(len_cutoff, 10))
	plt.axvline(cov_cutoff)
	plt.savefig('{}.png'.format(prefix),dpi=300)


if __name__ == '__main__':
	main()
