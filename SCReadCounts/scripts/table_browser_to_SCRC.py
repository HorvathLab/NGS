#!/usr/bin/python
import sys


def write_out_file(out_file, ch, pos, s):
	with open(out_file, "w") as o:
		bases = ['A', 'T', 'G', 'C']
		for p, val in s:
			for b in bases:
				if b != val:
					o.write("{}\t{}\t{}\t{}".format(ch, pos+p, val, b))


def main(fin, fout):
	with open(fin, 'r') as f:
		seq = ""
		for line in f:
			if line.startswith(">"):
				a = line.split()[1].split("=")[1]
				if seq != "":
					write_out_file(fout, chrom, exon_start, seq)
				chrom = a.split(":")[0]
				exon_start = a.split(":")[1].split("-")[0]
				# snp_end = a.split(":")[1].split("-")[1]
				seq = ""
			else:
				seq.join(line)


if __name__ == '__main__':
	main(fin=sys.argv[1], fout=sys.argv[2])
