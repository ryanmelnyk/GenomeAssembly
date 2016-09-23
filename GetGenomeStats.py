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
	parser.add_argument('prefix',type=str,help='prefix for output files')
	return parser.parse_args()

def parse_seqs(fastafile):
	contigcount = 0
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
	d = {'log10length':pd.Series(leng),'coverage':pd.Series(cov)}
	return pd.DataFrame(d)

def main():
	args = parse_args()
	fastafile = os.path.abspath(args.fastafile)
	prefix = args.prefix

	df = parse_seqs(fastafile)
	fig = plt.figure()
	df.plot(x="coverage",y="log10length",kind="scatter")
	plt.savefig('{}.png'.format(prefix),dpi=300)


if __name__ == '__main__':
	main()
