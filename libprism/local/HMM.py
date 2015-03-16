import os

from math import log, exp
from numpy import logaddexp
from libprism.local import cloud as cloud_lib
from libprism.common.math import logsumexp

import numpy as np

from scipy.weave import inline, converters
from libprism.local.prepare import verify_complexity, blocks_from_clouds, contigs_from_clouds

class HMM(object):
    def __init__(self,data_start,data_end, all_clouds, clouds_at_index):
        self.length = data_end - data_start + 1
        self.data_start = data_start
        self.data_end = data_end
        self.matrix = list()

        self.y_states = ((0,1),(1,0))
        self.fwd_values = list()
        self.bwd_values = list()

        self.dir = os.path.dirname(os.path.realpath(__file__))

        # parse input clouds
        contigs_at_index = contigs_from_clouds(all_clouds, self.data_start, self.data_end)

        for i,cloud in enumerate(all_clouds):
            self.matrix.append(cloud)

        # segments_at_position is an instance variable used by
        # the hmm algorihms; right now, we populate it with contigs
        self.segments_at_position = contigs_at_index

        # build blocks:
        self.blocks = blocks_from_clouds(all_clouds, clouds_at_index)

        # verify that blocks are disjoint:
        for i, block0 in enumerate(self.blocks):
            for j, block1 in enumerate(self.blocks):
                if i<=j: continue
                pos0 = set(block0)
                pos1 = set(block1)
                common_pos = pos0 & pos1
                if len(common_pos) > 0: keyboard()

        verify_complexity(contigs_at_index, self.data_start, self.data_end)

        uncovered_positions = set()
        # detect positions that have no coverage:
        for j in range(self.data_start, self.data_end+1):
            if len(clouds_at_index[j]) == 0:
                uncovered_positions.add(j)

        self.uncovered_positions = uncovered_positions

    #noinspection PyUnusedLocal
    def run_viterbi(self,block):
        start, end = block[0], block[-1]
        block_length = len(block)

        back_pointers = list()
        M, M_prev = np.array([]), np.array([])
        N, N_prev = 0, 0

        with open(self.dir + "/c/viterbi.c") as code_file:
            viterbi_code = code_file.read()

        # compute probabilities and pointers:
        for j, j_prev in ((block[k], block[(k-1) % block_length])
                          for k in xrange(len(block))):

#            if j>start:
#                n =  2**len(self.segments_at_position[j]) *\
#                     2**len(self.segments_at_position[j_prev] -\
#                            self.segments_at_position[j])
#                print "vit", j, n

            N = len(self.segments_at_position[j])
            N_prev = len(self.segments_at_position[j_prev])
            np_free_indices = self.get_free_indices(j, j_prev)
            N_free = len(np_free_indices)

            M_prev = M
            # print N
            M = np.zeros(2**(N+1))
            bpointers = np.zeros(2**(N+1))
            back_pointers.append(bpointers)

            np_prob0, np_prob1 = self.get_prob_vectors(j)
            np_current_to_prev_pos = self.get_position_map(j, j_prev)
            np_clouds_at_position = self.get_clouds_at_position(j)
            Nclouds = len(np_clouds_at_position)

            start_position = int(j == start)

            inline(viterbi_code,
                ['j', 'M', 'M_prev', 'bpointers',
                 'np_prob0', 'np_prob1',
                 'np_current_to_prev_pos', 'np_clouds_at_position',
                 'np_free_indices',
                 'N', 'N_prev', 'Nclouds', 'N_free', 'start_position'],
                headers = ["<math.h>",
                           "\"" + self.dir + "/c/common.h\""],
                type_converters = converters.blitz,
                compiler = "gcc")

        ptr = np.argmax(M)
        haplotype = list()
        haplotype.append(self.num_to_ys(end, ptr))

        k = len(block)-1
        while k > 0:
            ptr = int(back_pointers[k][ptr])
            k -= 1
            haplotype.append(self.num_to_ys(block[k], ptr))
        haplotype.reverse()

        return haplotype

    #noinspection PyUnusedLocal,PyTypeChecker
    def run_forwards(self,block):
        start, end = block[0], block[-1]
        block_length = len(block)

        fwd_values = dict()
        fwd, fwd_prev = np.array([]), np.array([])
        N, N_prev = 0, 0

        with open(self.dir + "/c/forwards.c") as code_file:
            forwards_code = code_file.read()

        # compute probabilities and pointers:
        for j, j_prev in ((block[k], block[(k-1) % block_length])
                          for k in xrange(len(block))):

            N = len(self.segments_at_position[j])
            N_prev = len(self.segments_at_position[j_prev])
            np_free_indices = self.get_free_indices(j, j_prev)
            N_free = len(np_free_indices)

            fwd_prev = fwd
            fwd = np.zeros(2**(N+1))

            np_prob0, np_prob1 = self.get_prob_vectors(j)
            np_current_to_prev_pos = self.get_position_map(j, j_prev)
            np_clouds_at_position = self.get_clouds_at_position(j)
            Nclouds = len(np_clouds_at_position)

            start_position = int(j == start)

            inline(forwards_code,
                ['j', 'fwd', 'fwd_prev',
                 'np_prob0', 'np_prob1',
                 'np_current_to_prev_pos', 'np_clouds_at_position',
                 'np_free_indices',
                 'N', 'N_prev', 'Nclouds', 'N_free', 'start_position'],
                headers = ["<math.h>",
                           "\"" + self.dir + "/c/common.h\""],
                type_converters = converters.blitz,
                compiler = "gcc")

            fwd_values[j] = fwd

        self.fwd_values = fwd_values
        self.log_data_prob = logsumexp(fwd)

    #noinspection PyUnusedLocal,PyTypeChecker
    def run_backwards(self,block):
        start, end = block[0], block[-1]
        block_length = len(block)

        bwd_values = dict()
        bwd, bwd_next = np.array([]), np.array([])
        N, N_next = 0, 0

        with open(self.dir + "/c/backwards.c") as code_file:
            backwards_code = code_file.read()

        # compute probabilities and pointers:
        for j, j_next in ((block[k], block[(k+1) % block_length])
                           for k in xrange(block_length-1,-1,-1)):
            N = len(self.segments_at_position[j])
            N_next = len(self.segments_at_position[j_next])
            np_free_indices = self.get_free_indices(j, j_next)
            N_free = len(np_free_indices)

            bwd_next = bwd
            bwd = np.zeros(2**(N+1))

            np_prob0, np_prob1 = self.get_prob_vectors(j_next)
            np_current_to_next_pos = self.get_position_map(j, j_next)
            np_clouds_at_position = self.get_clouds_at_position(j_next)
            Nclouds = len(np_clouds_at_position)

            end_position = int(j == end)

            inline(backwards_code,
                ['j', 'bwd', 'bwd_next',
                 'np_prob0', 'np_prob1',
                 'np_current_to_next_pos', 'np_clouds_at_position',
                 'np_free_indices',
                 'N', 'N_next', 'Nclouds', 'N_free', 'end_position'],
                headers = ["<math.h>",
                           "\"" + self.dir + "/c/common.h\""],
                type_converters = converters.blitz,
                compiler = "gcc")

            bwd_values[j] = bwd

        self.bwd_values = bwd_values

    def get_posterior_code(self,k,block,y,s):
#        print self.fwd_values[block[k]][self.ys_to_num(block[k], y, s)], \
#              self.bwd_values[block[k]][self.ys_to_num(block[k], y, s)]
        return exp(
            self.fwd_values[block[k]][self.ys_to_num(block[k], y, s)] +
            self.bwd_values[block[k]][self.ys_to_num(block[k], y, s)] -
            self.log_data_prob)

    def get_posterior_code_y_only(self,k,block,y):
        N = 2**len(self.segments_at_position[block[k]])
        if y == (0,1): slice = range(0,N)
        else: slice = range(N,2*N)

        return logsumexp(self.fwd_values[block[k]][slice] +
                         self.bwd_values[block[k]][slice]) - \
                         self.log_data_prob

    #noinspection PyUnusedLocal
    def get_posterior_transition_prob_y_only(self,k,block,y_state,
                                             y_state_prev):
        j, j_prev = block[k], block[k-1]
        fwd = self.fwd_values[j_prev]
        bwd = self.bwd_values[j]
        log_data_prob = self.log_data_prob

        np_free_indices = self.get_free_indices(j, j_prev)
        np_prob0, np_prob1 = self.get_prob_vectors(j)
        np_current_to_prev_pos = self.get_position_map(j, j_prev)
        np_clouds_at_position = self.get_clouds_at_position(j)

        N = len(self.segments_at_position[j])
        N_prev = len(self.segments_at_position[j_prev])
        N_free = len(np_free_indices)
        Nclouds = len(np_clouds_at_position)

        if y_state == (0,1): y = 0
        else: y = 1
        if y_state_prev == (0,1): y_prev = 0
        else: y_prev = 1

        with open(self.dir +
                  "/c/get_posterior_joint_probability_for_y.c") as code_file:
            code = code_file.read()

        joint_lik = inline(code,
            ['j', 'fwd', 'bwd', 'log_data_prob',
             'y', 'y_prev',
             'np_prob0', 'np_prob1',
             'np_current_to_prev_pos', 'np_clouds_at_position',
             'np_free_indices',
             'N', 'N_prev', 'Nclouds', 'N_free'],
            headers = ["<math.h>",
                       "\"" + self.dir + "/c/common.h\""],
            type_converters = converters.blitz,
            compiler = "gcc")

        return (exp(joint_lik - self.get_posterior_code_y_only(k-1,block,
                                                                y_state_prev)))


    # HELPER METHODS:

    def get_prob_vectors(self, j):
        M = self.matrix
        segments = sorted(list(self.segments_at_position[j]))

#        return np.array([M[s].lik(j,0) for s in segments]),\
#               np.array([M[s].lik(j,1) for s in segments])

        # this is for consistency with the previous (i.e. python) version of
        # the code:
        return np.array([log(max(M[s].prob(j,0),1e-14)) for s in segments]),\
               np.array([log(max(M[s].prob(j,1),1e-14)) for s in segments])

    def get_position_map(self, j, j_prev):
        segments = sorted(list(self.segments_at_position[j]))
        segments_prev = sorted(list(self.segments_at_position[j_prev]))

        seg_to_pos = dict()
        for i,s_prev in enumerate(segments_prev):
            seg_to_pos[s_prev] = i

        return np.array([seg_to_pos.get(s,-1) for s in segments])

    def get_clouds_at_position(self, j):
        M = self.matrix
        segments = sorted(list(self.segments_at_position[j]))

        return np.array([i for i,c in enumerate(segments) if M[c][j] != -1])

    def get_free_indices(self, j, j_prev):
        segments_prev = sorted(list(self.segments_at_position[j_prev]))
        return np.array([i for i,c in enumerate(segments_prev)
                           if c not in self.segments_at_position[j]])

    def num_to_ys(self, j, S):
        segments = sorted(list(self.segments_at_position[j]))
        N = 2**len(segments)

        if S / N == 0: y = ((0,1))
        else: y = ((1,0))

        s = frozenset([seg for i,seg in enumerate(segments)
                       if (S >> i) % 2 == 1])
        return y,s

    def ys_to_num(self, j, y, s):
        segments = sorted(list(self.segments_at_position[j]))
        N = 2**len(segments)

        S = 0
        for i, seg in enumerate(segments):
            if seg in s: S |= 1 << i

        if y == (1,0): S += N

        return S

    def log_emission_prob(self,j,y,s):
        matrix = self.matrix
        strand1_segments = set(s)
        strand0_segments = self.segments_at_position[j] - strand1_segments

        log_prob = 0
        for seg in strand0_segments:
            if matrix[seg][j] != -1:
                prob_match = max(matrix[seg].prob(j,matrix[seg][j]),1e-14)
                prob_mismatch = max(matrix[seg].prob(j,(matrix[seg][j]+1)%2),1e-14)
                if matrix[seg][j] == y[0]: 
                    log_prob += log(prob_match)
                else: 
                    log_prob += log(prob_mismatch)
        for seg in strand1_segments:
            if matrix[seg][j] != -1:
                prob_match = max(matrix[seg].prob(j,matrix[seg][j]),1e-14)
                prob_mismatch = max(matrix[seg].prob(j,(matrix[seg][j]+1)%2),1e-14)
                if matrix[seg][j] == y[1]: 
                    log_prob += log(prob_match)
                else: 
                    log_prob += log(prob_mismatch)

        return log_prob
