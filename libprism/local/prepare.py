#!/usr/bin/python

"""
Functions for loading and pre-processing clouds for the local HMM.

This module defines several types of constructs:
    Clouds : Objects defined in cloud.py that span a certain set of indices,
    and that contain alleles and probabilities.
    Segments : A contiguous interval of positions [start, end] associated with a cloud.
    Blocks : Set of positions that is assembled from clouds.
"""

from libprism.local.cloud import Cloud
import libprism.local.cloud as cloud_lib

###############################################################################
## LOAD CLOUDS

def parse_metadata(metadata):
    """
    Returns a dictionary that contains the metadata associated with a cloud.
    Metadata is assumed to come from the .bed file.
    """
    tags = dict()
    fields = metadata.split(':')
    for field in fields:
        key, value = field.split('=')
        tags[key] = value

    return tags

def clouds_from_cld(cloud_filename, data_start, data_end):
    cloud_file = open(cloud_filename)

    # we will organize clouds into the following:
    clouds_at_index = dict()
    all_clouds = set()

    phaseable_positions = set() # for debugging

    # initialize the above:
    for j in xrange(data_start,data_end+1):
        clouds_at_index[j] = set()

    print "Parsing cld file..."
    for i,line in enumerate(cloud_file):
        fields = line.strip().split()

        if len(fields) == 2: continue

        # extract data for the cloud:
        covered_positions = list()
        alleles = dict()
        qscores = dict()

        start = int(fields[2])
        allele_list = fields[3]
        qscore_list = fields[4]

        for j_index, j in enumerate(xrange(start + len(allele_list))):
            if j != '-':
                alleles[j] = [int(allele_list[j_index])]
                x = qscore_list[j_index]
                qscores[j] = [float(ord(x))-33]
                covered_positions.append(j)

        if len(covered_positions) == 0:
            continue

        covered_positions.sort()

        cloud_start = covered_positions[0]
        cloud_end = covered_positions[-1]
        cloud_length = cloud_end - cloud_start + 1

        # skip clouds outisde our interval:
        if cloud_start < data_start or cloud_end > data_end:
            continue

        # skip 1-snp clouds:
        if cloud_length == 1: continue

        # create the cloud object
        cloud = Cloud(cloud_start, cloud_end)
        cloud.name = fields[2]

        cloud.from_01_bed(covered_positions, alleles, qscores)

        phaseable_positions.update(set(cloud.positions))

        if len(cloud.positions) == 0:
            print "WARNING: All positions in cloud %s were filtered out" % cloud.name

        all_clouds.add(cloud)
        # place it into our sets:
        for j in xrange(cloud_start, cloud_end+1):
            if cloud[j] != -1:
                clouds_at_index[j].add(cloud)

    cloud_file.close()

    print 'Positions that can be phased:', len(phaseable_positions)

    return all_clouds, clouds_at_index

def clouds_from_refhap(cloud_filename, data_start, data_end):
    cloud_file = open(cloud_filename)

    # we will organize clouds into the following:
    clouds_at_index = dict()
    all_clouds = set()

    phaseable_positions = set() # for debugging

    # initialize the above:
    for j in xrange(data_start,data_end+1):
        clouds_at_index[j] = set()

    print "Parsing Refhap file..."
    for i,line in enumerate(cloud_file):
        fields = line.strip().split()

        if len(fields) == 2: continue

        # extract data for the cloud:
        covered_positions = list()
        alleles = dict()
        qscores = dict()

        qscore_list = list(fields[-1][::-1])
        num_pieces = int(fields[0])

        for k in xrange(num_pieces):
            start = int(fields[2+2*k])
            allele_list = fields[2+2*k+1]

            for j, allele in zip(xrange(start,start+len(allele_list)), allele_list):
                alleles[j] = [int(allele)]
                x = qscore_list.pop()
                qscores[j] = [float(ord(x)-33)]
                covered_positions.append(j)

        if len(covered_positions) == 0:
            continue

        covered_positions.sort()

        cloud_start = covered_positions[0]
        cloud_end = covered_positions[-1]
        cloud_length = cloud_end - cloud_start + 1

        # skip clouds outisde our interval:
        if cloud_start < data_start or cloud_end > data_end:
            continue

        # skip 1-snp clouds:
        if cloud_length == 1: continue

        # create the cloud object
        cloud = Cloud(cloud_start, cloud_end)
        cloud.name = fields[1]

        cloud.from_01_bed(covered_positions, alleles, qscores)

        phaseable_positions.update(set(cloud.positions))

        if len(cloud.positions) == 0:
            print "WARNING: All positions in cloud %s were filtered out" % cloud.name
            continue

        all_clouds.add(cloud)
        # place it into our sets:
        for j in xrange(cloud_start, cloud_end+1):
            if cloud[j] != -1:
                clouds_at_index[j].add(cloud)

    cloud_file.close()

    print 'Positions that can be phased:', len(phaseable_positions)

    return all_clouds, clouds_at_index
    return all_clouds, clouds_at_index

def print_clouds(f, clouds, data_start, data_end):
    f = open(f,"w")
    f.write('%d %d\n' % (len(clouds), data_end+1))
    for cloud in clouds:
        f.write("1 %s %d %s\n" % (cloud.name, cloud.start, cloud.to_str()))
    f.close()

###############################################################################
## PRE-PROCESS CLOUDS

def merge_clouds(all_clouds, clouds_at_index, data_start, data_end, max_coverage=10000):
    print "Merging clouds..."

    contigs_at_index = dict()
    for j in xrange(data_start, data_end+1):
        contigs_at_index[j] = set()
    for cloud in all_clouds:
        for j in xrange(cloud.start, cloud.end):
            contigs_at_index[j].add(cloud)

    for j in xrange(data_start+1, data_end+1):
        #if clouds_at_index[j].isdisjoint(clouds_at_index[j-1]): continue
        if len(contigs_at_index[j]) + len(contigs_at_index[j-1] - contigs_at_index[j]) > 16:
            # merge clouds that are similar enough:
            # define 3 stacks:
            comparing = list(contigs_at_index[j])
            compared = list()
            final = list()

            # this will mark clouds for deletion:
            clouds_to_delete = set()

            # keep merging clouds until all clouds have been compared:
            while len(comparing) + len(compared) > 0:
                # pick 1st cloud to which others will be compared
                current_cloud = comparing.pop()

                # while there are clouds left to be compared
                while len(comparing) > 0:
                    # pick a cloud to compare
                    potential_cloud = comparing.pop()

                    # if it is similar enough
                    if current_cloud[j] != -1 and current_cloud.matches(potential_cloud):
                        # merge them into a new cloud
                        new_cloud = cloud_lib.from_two_clouds(current_cloud,potential_cloud)

                        # tag new clouds for deletion:
                        clouds_to_delete.add(current_cloud)
                        clouds_to_delete.add(potential_cloud)

                        # make the new cloud the one we're comparing to
                        current_cloud = new_cloud
                    else:
                        compared.append(potential_cloud)

                # we will not be able to merge any more clouds,
                # so add current_cloud to the final set of clouds
                final.append(current_cloud)

                # switch the compared and comparing stacks
                comparing, compared = compared, comparing

            # print len(final)
            if len(final) > max_coverage:
                sorted_final = sorted(final, key=lambda x: len(x), reverse=True)
                final = set(sorted_final[0:max_coverage])
                clouds_to_delete |= set(sorted_final[max_coverage:])
            # print len(final)

            # remove the clouds marked for deletion:
            for c in clouds_to_delete:
                if c in all_clouds:
                    all_clouds.remove(c)
                    for k in xrange(c.start, c.start + len(c)):
                        if c in clouds_at_index[k]:
                            clouds_at_index[k].remove(c)
                        if c in contigs_at_index[k]:
                            contigs_at_index[k].remove(c)

            for cloud in final:
                all_clouds.add(cloud)
                for k in xrange(cloud.start, cloud.end+1):
                    if cloud[k] != -1:
                        clouds_at_index[k].add(cloud)
                    contigs_at_index[k].add(cloud)

        assert len(contigs_at_index[j]) <= max_coverage

    return all_clouds, clouds_at_index


def verify_complexity(clouds_at_position, data_start, data_end):
    first_warning = True

    for j in xrange(data_start+1,data_end+1):
        if clouds_at_position[j].isdisjoint(clouds_at_position[j-1]):
            continue

        n = len(clouds_at_position[j]) + len(clouds_at_position[j-1]\
                                             - clouds_at_position[j])
        if 2**n > 1000000:
            if first_warning:
                print "WARNING: The following positions have > 10e6 "\
                      "inner loop iterations:"
            print "\t%d:\t%d" % (j, n)

def contigs_from_clouds(all_clouds, data_start, data_end):
    # a contig is the [cloud_start, cloud_end] interval associated
    # with a cloud

    contigs_at_index = dict()
    for j in xrange(data_start,data_end+1):
        contigs_at_index[j] = set()

    for i,cloud in enumerate(all_clouds):
        cloud_start = cloud.start       # cloud positions start at 0
        cloud_length = len(cloud)

        # place it into our sets:
        for k in xrange(cloud_length):
            contigs_at_index[cloud_start+k].add(i)

    return contigs_at_index

###############################################################################
## CREATE BLOCKS

def blocks_from_contigs(contigs_at_position, data_start, data_end):
    print "Separating into blocks..."
    blocks = list()
    cur_start = data_start
    j = None
    for j in range(data_start, data_end):
        if contigs_at_position[j].isdisjoint(contigs_at_position[j + 1]):
            # if no segment as pos. j extends to pos. j+1, end block
            blocks.append((cur_start, j))
            cur_start = j + 1
            # add the very last block:
    if j:
        blocks.append((cur_start, j + 1))

    return blocks

def blocks_from_clouds(all_clouds, clouds_at_index):
    """
    Constructs blocks from clouds. A block is a list of positions such
    that no cloud over these positions leaves the block.
    """
    blocks = list()
    cloud_stack = set(all_clouds)
    all_visited_positions = set()
    positions_by_cloud = dict()
    processed_clouds = set()

    for cloud in all_clouds:
        cloud_pos = cloud.positions
        for j in cloud_pos:
            if cloud not in clouds_at_index[j]:
                exit("ERROR: clouds improperly placed in sets")


    while len(cloud_stack) > 0:
        initial_cloud = cloud_stack.pop()
        processed_clouds.add(initial_cloud)

        # extend initial_cloud to a block:
        new_positions = set(initial_cloud.positions)


        visited_positions = set()

        while len(new_positions) > 0:
            j = new_positions.pop()
            visited_positions.add(j)
            assert j not in all_visited_positions

            for c in clouds_at_index[j]:
                if c in cloud_stack:
                    cloud_stack.remove(c)
                    processed_clouds.add(c)
                    new_positions |= c.positions
                    for p in c.positions:
                        if p not in positions_by_cloud:
                            positions_by_cloud[p] = c

        block_positions = list(visited_positions)
        block_positions.sort()
        if len(visited_positions) == 0: exit("Block with no positions")
        blocks.append(block_positions)
        all_visited_positions.update(visited_positions)

    blocks.sort()
    return blocks
