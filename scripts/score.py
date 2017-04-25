#python3

import argparse
import csv
import re
from Bio import SeqIO

def process_reference():
	# setup reference
	reference_file = "ref.fasta"
	reference = {} # dict
	rid = re.compile(r'^genomeref\.(\d+)$')
	for record in SeqIO.parse(reference_file, "fasta"):
		# print(record.id)
		m = rid.match(record.id)
		id = m.group(1)
		length = len(record.seq)

		reference[id] = [0] * length # empty list of length `length`

	return reference
	

def process_fasta(query_file):
	reference = process_reference()
	rqid = re.compile(r'\.(\d+)$') # depends on input format

	# also may need to parse FASTA comment
	for record in SeqIO.parse(query_file, "fasta"):
		res = (r'\.(\d+)$', "\|(\d+)" )

		id = ""
		for r in res:
			m = re.search(r, record.id)
			if m != None:
				id = m.group(1)

		if id != "":
			# fast path
			if len(record.seq) == len(reference[id]):
				# this is good
				for i in range(0, len(reference[id])):
					reference[id][i] += 1
			else:
				# TODO: extract position information
				pass

	return reference


def output( reference):
	# compute sens and spez
	total_sites = 0
	correct_sites = 0
	missed_sites = 0
	extra_sites = 0
	for locus_id in reference:
		for site in reference[locus_id]:
			total_sites += 1
			if site == 0:
				missed_sites += 1
			elif site == 1:
				correct_sites += 1
			elif site > 1:
				extra_sites += site - 1

	print("total_sites:\t{}".format(total_sites))
	print("correct_sites:\t{} \t{}%".format(correct_sites, int(correct_sites/total_sites * 1000)/10))
	print("missed_sites:\t{} \t{}%".format(missed_sites, int(missed_sites/total_sites * 1000)/10))
	print("extra_sites:\t{} \t{}%".format(extra_sites, int(extra_sites/total_sites * 1000)/10))


def process_spine(path):
	reference = process_reference()
	# open file
	with open(path,'r') as f:
		next(f) # skip header
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			if len(row) != 6:
				continue

			m = re.search("\.(\d+)$", row[0])
			locus_id = m.group(1)
			start = int(row[2]) - 1
			end = int(row[3])
			for i in range(start, end):
				reference[locus_id][i] += 1

	return reference


def main():
	parser = argparse.ArgumentParser(description='Process some files.')
	parser.add_argument("-f", '--fasta', dest='fasta', 
		default="", help='the query file (fasta format)')
	parser.add_argument('--spine', dest='spine', 
		default="", help='the query file (spine coords format)')
	args = parser.parse_args()

	ref = None
	if args.fasta != "":
		ref = process_fasta(args.fasta)

	if args.spine != "":
		ref = process_spine(args.spine)

	output(ref)


if __name__ == '__main__':
	main()
