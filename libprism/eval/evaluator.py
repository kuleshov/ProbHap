#!/usr/bin/python
import os
import re
import argparse

from libprism.common.stats import N50

###############################################################################

CHR_LENGTHS={
 'chr10': 135534747,
 'chr11': 135006516,
 'chr12': 133851895,
 'chr13': 115169878,
 'chr14': 107349540,
 'chr15': 102531392,
 'chr16': 90354753,
 'chr17': 81195210,
 'chr18': 78077248,
 'chr19': 59128983,
 'chr1': 249250621,
 'chr20': 63025520,
 'chr21': 48129895,
 'chr22': 51304566,
 'chr23': 2639519,
 'chr24': 329516,
 'chr25': 152231523,
 'chr2': 243199373,
 'chr3': 198022430,
 'chr4': 191154276,
 'chr5': 180915260,
 'chr6': 171115067,
 'chr7': 159138663,
 'chr8': 146364022,
 'chr9': 141213431
}

###############################################################################

def evaluate(master_file, phase_file, post_thr, emi_thr, tr_thr):

	true_phase = dict()
	# load true phase
	with open(master_file) as master:
		for line in master:
			fields = line.strip().split()
			if len(fields) == 5 and fields[-1] == 'N':
			  continue
                        j = int(fields[0])
                        true_phase[j] = int(fields[2]), int(fields[3])

	def correct_pattern(pattern):
	    ones_at_ends_pattern=re.compile('^10|01$')
	    zeros_at_ends_pattern=re.compile('^01|10$')
	    pattern = ones_at_ends_pattern.sub('00', pattern)
	    pattern = zeros_at_ends_pattern.sub('11', pattern)

	    fives_pattern0 = re.compile('01010')
	    fives_pattern1 = re.compile('10101')
	    pattern = fives_pattern0.sub('00000', pattern)
	    pattern = fives_pattern1.sub('11111', pattern)

	    single_one_pattern = re.compile('010')
	    single_zero_pattern = re.compile('101')
	    pattern = single_one_pattern.sub('000', pattern)
	    pattern = single_zero_pattern.sub('111', pattern)

	    return pattern

	# evaluate algorithm positions:
	def evaluate_block(block, true_phase):
		j_indices = sorted([j for j in current_block.keys() if j in true_phase
							and block[j][0] != '-' and block[j][1] != '-'])

		if len(j_indices) < 2:
			return 0, 0, 0

		if int(block[j_indices[0]][0]) == true_phase[j_indices[0]][0]:
			sw = 0
			pattern = '0'
		else:
			sw = 1
			pattern = '1'

		for j in j_indices[1:]:
			if (int(block[j][0])+sw)%2 != true_phase[j][0]:
				sw = (sw+1)%2

			if int(block[j][0]) == true_phase[j][0]:
				pattern += '0'
			else:
				pattern += '1'

		num_sw = sum([1 for i in xrange(1,len(pattern)) if pattern[i] != pattern[i-1]])
		pattern2 = correct_pattern(pattern)
		num_long_sw = sum([1 for i in xrange(1,len(pattern2)) if pattern2[i] != pattern2[i-1]])
		num_pos = len(pattern) - 1

		return num_sw, num_long_sw, num_pos

	total_sw = 0
	total_long_sw = 0
	total_pos = 0
	phased_positions = set()
	block_lengths = list()

	with open(phase_file) as f:
		current_block = dict()
		for line in f:
			if line.startswith('BLOCK'):
				continue

			fields = line.strip().split()

			if line.startswith('****') or (tr_thr and float(fields[3]) < tr_thr and float(fields[4]) > 0.9):

				# if not line.startswith('****'):
				# 	print '>', line

				if len(current_block) >= 2:
					num_sw, num_long_sw, num_pos = evaluate_block(current_block, true_phase)
					total_sw += num_sw
					total_long_sw += num_long_sw
					total_pos += num_pos
					phased_positions.update(current_block.keys())
					block_lengths.append(len(current_block.keys()))

				current_block = dict()
				continue

			if fields[1] == '-' or fields[2] == '-':
				continue

			if post_thr and emi_thr:
				if float(fields[5]) < emi_thr or float(fields[4]) < post_thr:
			# # if float(fields[4]) < 0.6:
			# 	# print line
					continue
			pos = int(fields[0])
			current_block[pos] = fields[1], fields[2]

	return total_sw, total_long_sw, total_pos, phased_positions, block_lengths, pos

###############################################################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate simulated .bed file")

	parser.add_argument('--master', required=True)
	parser.add_argument('--alg-phase', required=True)
	parser.add_argument('--chr', required=True)
	parser.add_argument('--post-thr', type=float)
	parser.add_argument('--emi-thr', type=float)
	parser.add_argument('--tr-thr', type=float)

	args = parser.parse_args()

	total_sw, total_long_sw, total_pos, phased_positions, block_lengths, total_num_positions = \
		evaluate(args.master, args.alg_phase, args.post_thr, args.emi_thr, args.tr_thr)
	
	total_short_sw = total_sw - total_long_sw

	print "sw_acc:", (total_pos - total_sw) / float(total_pos), total_pos - total_sw, total_pos
	print "long_sw_acc:", (total_pos - total_long_sw) / float(total_pos), total_pos - total_long_sw, total_pos
	print "short_sw_acc:", (total_pos - total_short_sw) / float(total_pos), total_pos - total_short_sw, total_pos
	print "sw_per_mb:", total_sw / (CHR_LENGTHS[args.chr] / 10.0**6), total_sw, CHR_LENGTHS[args.chr] / 10.0**6
	print "long_sw_per_mb:", total_long_sw / (CHR_LENGTHS[args.chr] / 10.0**6), total_long_sw, CHR_LENGTHS[args.chr] / 10.0**6
	print "short_sw_per_mb:", total_short_sw / (CHR_LENGTHS[args.chr] / 10.0**6), total_short_sw, CHR_LENGTHS[args.chr] / 10.0**6
	print "percent_phased:", len(phased_positions) / float(total_num_positions), len(phased_positions), total_num_positions
	print "N50:", N50(block_lengths)
