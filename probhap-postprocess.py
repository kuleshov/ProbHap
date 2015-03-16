#!/usr/bin/python
import argparse

#############################################################################

parser = argparse.ArgumentParser(description="")

parser.add_argument('--filtered-reads', help='path to filtered reads file', required=True)
parser.add_argument('--assignments', help='reads to the parent of origin', required=True)
parser.add_argument('--blocks', required=True)
parser.add_argument('--corrected-blocks', required=True)

args = parser.parse_args()

#############################################################################

# figure out the number of positions to phase:
with open(args.filtered_reads) as f:
    num_reads, num_positions = (int(x) for x in f.readline().split())
start, end = 1, num_positions

assignments = dict()
with open(args.assignments) as f:
	for line in f:
		fields = line.strip().split()
		assignments[fields[0]] = int(fields[1])

zero_support_set = {j:set() for j in xrange(start, end+1)}
one_support_set = {j:set() for j in xrange(start, end+1)}

with open(args.filtered_reads) as read_file:
	for line in read_file:
		fields = line.strip().split()
		if len(fields) <= 2: continue
		read = fields[1]
		start_j = int(fields[2])

		if read not in assignments:
			print 'WARNING: Could not find read', read
			continue

		for k, s in enumerate(fields[3]):
			if s == '-': continue
			if int(s) != assignments[read]:
				# s == 0 and ass = 0 or s == 1 and ass = 1
				zero_support_set[start_j+k].add(read)
			else:
				# s == 0 and ass = 1 or s == 1 and ass = 0
				one_support_set[start_j+k].add(read)

consensus = dict()
for j in xrange(start, end+1):
	if len(zero_support_set[j]) > len(one_support_set[j]):
		consensus[j] = 0
	elif len(zero_support_set[j]) < len(one_support_set[j]):	
		consensus[j] = 1

new_blocks = open(args.corrected_blocks, 'w')
with open(args.blocks) as f:
	for line in f:
		if line.startswith('BLOCK') or line.startswith('***'):
			new_blocks.write(line)
		else:
			fields = line.split()
			j = int(fields[0])
			if consensus.get(j, -1) == 0:
				new_blocks.write('%d\t0\t1\t%s\t%s\t%s\n' % (j, fields[3], fields[4], fields[5]))
			elif consensus.get(j, -1) == 1:
				new_blocks.write('%d\t1\t0\t%s\t%s\t%s\n' % (j, fields[3], fields[4], fields[5]))
			else:
				new_blocks.write('%d\t-\t-\t%s\t%s\t%s\n' % (j, fields[3], fields[4], fields[5]))

new_blocks.close()
