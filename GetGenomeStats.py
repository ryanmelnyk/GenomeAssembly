#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import argparse, os
from Bio import SeqIO
from Bio import SeqUtils
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
	parser.add_argument('outdir', type=str,help='location of SPAdes output - prokka folder should be within')
	parser.add_argument('len_cutoff', type=int,help='cutoff for length of contigs')
	parser.add_argument('cov_cutoff', type=float,help='cutoff for coverage of contigs')
	parser.add_argument('prefix',type=str,help='prefix for output files')
	return parser.parse_args()

def parse_seqs(outdir,prefix,len_cutoff,cov_cutoff):

	o = open(prefix+".filteredcontigs.fna",'w')
	leng = {}
	cov = {}
	GC_cont = {}

	hits = {}
	strain = os.path.basename(outdir)
	for line in open(os.path.join(outdir,"{}.prokka".format(strain),"{}.m8".format(strain)),'r'):
		vals = line.rstrip().split("\t")
		if vals[0] not in hits:
			hits[vals[0]] = (vals[1],float(vals[11]))
		elif float(vals[11]) > hits[vals[0]][1]:
			hits[vals[0]] = (vals[1],float(vals[11]))
		else:
			pass

	contigs = {}
	for line in open(os.path.join(outdir,"{}.prokka".format(strain),"{}.gff".format(strain)),'r'):
		if line.startswith("##"):
			continue
		elif line.startswith(">"):
			break
		else:
			vals = line.rstrip().split("\t")
			if vals[2] == "CDS":
				if vals[0] not in contigs:
					contigs[vals[0]] = [vals[8].split(";")[0].split("=")[1]]
				else:
					contigs[vals[0]].append(vals[8].split(";")[0].split("=")[1])

	annotations = {}
	for c in contigs:
		genera = {}
		name = c.split("_")
		contigname = "{}_{}".format(name[0],name[1])
		for tag in contigs[c]:
			if tag in hits:
				vals = hits[tag][0].split("_")
				if vals[0] == "":
					if vals[1] not in genera:
						genera[vals[1]] = 1
					else:
						genera[vals[1]] += 1
				else:
					if vals[0] not in genera:
						genera[vals[0]] = 1
					else:
						genera[vals[0]] += 1
		try:
			annotations[contigname] = sorted(genera.items(), reverse = True, key= lambda (k,v): v)[0][0]
		except IndexError:
			annotations[contigname] = "None"

	annot_counts = {}
	for a in annotations:
		if annotations[a] not in annot_counts:
			annot_counts[annotations[a]] = 1
		else:
			annot_counts[annotations[a]] += 1

	if len(annot_counts.keys()) > 5:
		top_annotations = [k[0] for k in sorted(annot_counts.items(), reverse = True, key = lambda (k,v): v)][0:5]
		for a in annotations:
			if annotations[a] not in top_annotations:
				if annotations[a] == "None":
					pass
				else:
					annotations[a] = "Other"

	for line in open(os.path.join(outdir,"{}.prokka".format(strain),"{}.phiX.hits".format(strain))):
		if line.startswith("#"):
			continue
		else:
			vals = line.rstrip().split()
			if float(vals[13]) > 100.0:
				name = vals[0].split("_")
				contigname = "{}_{}".format(name[0],name[1])
				annotations[contigname] = "phiX"

	contig_lengths = []
	for seq in SeqIO.parse(open(os.path.join(outdir,"contigs.fasta"),'r'),'fasta'):
		vals = seq.id.split("_")
		contigname = "{}_{}".format(vals[0],vals[1])
		if float(vals[5]) > 1000:
			continue
		else:
			GC_cont[contigname] = SeqUtils.GC(seq.seq)
			leng[contigname] = math.log(int(vals[3]),10)
			cov[contigname] = float(vals[5])
		if contigname not in annotations:
			annotations[contigname] = "None"
		if leng[contigname] > math.log(len_cutoff,10) and cov[contigname] > cov_cutoff:
			if annotations[contigname] == "phiX":
				continue
			else:
				contig_lengths.append(len(str(seq.seq)))
				seq.id = "{}_contig{}".format(os.path.basename(prefix),vals[1])
				seq.description = ""
				SeqIO.write(seq,o,'fasta')
	d = {'log10length':pd.Series(leng),'coverage':pd.Series(cov),'GC_content':pd.Series(GC_cont), 'annotations':pd.Series(annotations)}
	o.close()
	stats = open("{}.stats.txt".format(prefix),'w')
	stats.write("Total sequence length\t{}\n".format(str(sum(contig_lengths))))
	stats.write("Number of contigs\t{}\n".format(str(len(contig_lengths))))
	n50_ct = 0
	for i in sorted(contig_lengths,reverse=True):
		n50_ct += i
		if n50_ct > sum(contig_lengths)/2:
			stats.write("N50\t{}\n".format(str(i)))
			break
	stats.close()
	return pd.DataFrame(d)

def main():
	args = parse_args()
	outdir = os.path.abspath(args.outdir)
	len_cutoff = args.len_cutoff
	cov_cutoff = args.cov_cutoff
	prefix = os.path.abspath(args.prefix)

	df = parse_seqs(outdir,prefix,len_cutoff,cov_cutoff)
	fig = plt.figure()
	df.plot(x="coverage",y="log10length",c="GC_content",kind="scatter",colormap="viridis")
	plt.axhline(math.log(len_cutoff, 10))
	plt.axvline(cov_cutoff)
	plt.savefig('{}.png'.format(prefix),dpi=300)
	df.plot(x="coverage",y="log10length",c="GC_content",kind="scatter",xlim=(0,80),colormap="viridis")
	plt.axhline(math.log(len_cutoff, 10))
	plt.axvline(cov_cutoff)
	plt.savefig('{}.lowcov.png'.format(prefix),dpi=300)
	df.plot(x="coverage",y="GC_content",c="log10length",kind="scatter",colormap="viridis")
	plt.savefig('{}.GC.png'.format(prefix),dpi=300)
	df.plot(x="coverage",y="GC_content",c="log10length",kind="scatter",xlim=(0,80),colormap="viridis")
	plt.savefig('{}.GC.lowcov.png'.format(prefix),dpi=300)
	sns.lmplot(x="coverage",y="GC_content",hue="annotations",data=df,fit_reg=False)
	plt.savefig('{}.GC2.png'.format(prefix),dpi=300)
	sns.lmplot(x="coverage",y="GC_content",hue="annotations",data=df,fit_reg=False)
	plt.xlim(0,80)
	plt.savefig('{}.GC2.lowcov.png'.format(prefix),dpi=300)

if __name__ == '__main__':
	main()
