#!/usr/bin/python
import argparse

from libprism.local import HMM
from libprism.local.prepare import clouds_from_refhap, merge_clouds, print_clouds

from math import log, exp

#############################################################################

parser = argparse.ArgumentParser(description="Single-individual haplotyping using probabilistic graphical models")

parser.add_argument('--reads', help='path to reads file', required=True)
parser.add_argument('--parsed-reads', help='path to parsed reads file', required=True)
parser.add_argument('--phase', help='output file', required=True)
parser.add_argument('--assignments', help='clouds to chromosomes', required=True)

args = parser.parse_args()

#############################################################################

# figure out the number of positions to phase:
with open(args.reads) as f:
    num_clouds, num_positions = (int(x) for x in f.readline().split())
data_start, data_end = 0, num_positions - 1

# load clouds
all_clouds, clouds_at_index = clouds_from_refhap(args.reads, data_start, data_end)
all_clouds, clouds_at_index = merge_clouds(all_clouds, clouds_at_index, data_start,
                                                       data_end, max_coverage=20)

all_clouds_list = list(all_clouds)

print_clouds(args.parsed_reads, all_clouds_list, data_start, data_end)

hmm = HMM.HMM(data_start,data_end, all_clouds_list, clouds_at_index)

M = hmm.matrix

out_file = open(args.phase, "w")
out_file.write("")
out_file.close()

assignments = open(args.assignments, 'w')
assignments.write("")
assignments.close()

print "Phasing blocks..."

# phase each locally phaseable block in turn

N_blocks = len(hmm.blocks)
for i, block in enumerate(hmm.blocks):
    if i % 100 == 0:
      print "Block %d/%d" % (i, N_blocks)
    start,end = block[0], block[-1]

    out_file = open(args.phase, "a")
    out_file.write("BLOCK: offset: %d len: %d positions: %d\n" % (start+1, end-start+1, len(block)))
    # print "BLOCK: offset: %d len: %d positions: %d" % (start+1, end-start+1, len(block))

    assignments = open(args.assignments, 'a')

    haplotype = hmm.run_viterbi(block)
    hmm.first_y = haplotype[0][0]

    right_segments = set()
    for (y,s) in haplotype: right_segments |= s

    left_clouds = set([M[i].name for j in block for i in hmm.segments_at_position[j]
                       if i not in right_segments and M[i][j] != -1])
    right_clouds = set([M[i].name for j in block for i in hmm.segments_at_position[j]
                        if i in right_segments and M[i][j] != -1])
    left_cloud_names = ','.join([cloud_name for cloud_name in left_clouds])
    right_cloud_names = ','.join([cloud_name for cloud_name in right_clouds])

    for c in left_clouds:
        assignments.write('%s\t0\n' % (c))

    for c in right_clouds:
        assignments.write('%s\t1\n' % (c))

    hmm.run_forwards(block)
    hmm.run_backwards(block)

    for k,(y,s) in enumerate(haplotype):
        j = block[k]
        out_file.write("%d\t" % (j))
        if j in hmm.uncovered_positions:
            out_file.write("-\t-\t")
        else:
            out_file.write("%d\t%d\t" % y)
        if k == 0:
            score = 0.5
        else:
            score = hmm.get_posterior_transition_prob_y_only(
                k, block, y, y_prev)
        emission_probability = exp(hmm.log_emission_prob(j, y, s))
        posterior_code = exp(hmm.get_posterior_code_y_only(k,block,y))
        out_file.write("%f\t%f\t%f\n" % (score,posterior_code,emission_probability))

        y_prev = y

    out_file.write("********\n")

    out_file.close()
    assignments.close()
