from numpy import logaddexp
from math import  log

def from_two_clouds(c1, c2):
    start = min(c1.start,c2.start)
    end = max(c1.start + len(c1) - 1, c2.start + len(c2) - 1)

    cloud = Cloud(start, end)

    cloud.name = c1.name + "," + c2.name

    for j in xrange(start,end+1):
        if c1[j] == -1 and c2[j] == -1:
            cloud.seq.append(-1)
        elif c1[j] == -1:
            cloud.seq.append(c2[j])

            if c2[j] == 1: cloud.one_positions.add(j)
            elif c2[j] == 0: cloud.zero_positions.add(j)
        elif c2[j] == -1:
            cloud.seq.append(c1[j])

            if c1[j] == 1: cloud.one_positions.add(j)
            elif c1[j] == 0: cloud.zero_positions.add(j)
        else:
            if c1[j] == c2[j]:
                cloud.seq.append(c1[j])
                v = c1[j]     
            else:
                if c1.prob(j,c1[j]) > c2.prob(j,c2[j]):
                    cloud.seq.append(c1[j])
                    v = c1[j]
                elif c1.prob(j,c1[j]) < c2.prob(j,c2[j]):
                    cloud.seq.append(c2[j])
                    v = c2[j]
                else:
                    cloud.seq.append(-1)
                    v = -1

            if v == 1: cloud.one_positions.add(j)
            elif v == 0: cloud.zero_positions.add(j)

        cloud.accuracies.append(
            merge_probabilities(c1.get_prob(j), 
                                c2.get_prob(j)))

    cloud.positions = cloud.one_positions | cloud.zero_positions
    if len(cloud.positions) < len(c1.positions & c2.positions):
        print "WARNING: A position was discarded after obtaining cloud %s from merging." % cloud.name

    return cloud

def merge_probabilities(p1,p2):
    alleles = set()
    alleles.update(p1.keys(), p2.keys())
    p = dict()

    for a in alleles:
        p[a] = 1.0
        if a in p1:
            p[a] = p[a] * p1[a]
        if a in p2:
            p[a] = p[a] * p2[a]

    return p

class Cloud():
    def __init__(self,cloud_start, cloud_end):
        self.seq = list()               # 0/1 values of cloud
        self.start = cloud_start        # pos in snp # coordinates 
                                        # (e.g. 0 to 22k+ for chr. 22):
        self.end = cloud_end            # end pos in snp # coordinates
        self.coverage = list()          # how many reads cover each position
        self.accuracy = list()          # prob. of call being correct
        self.accuracies = list()
        self.one_positions = set()      # positions that contain ones
        self.zero_positions = set()     # positions that contain zeros
        self.positions = set()          # all positions covered by segment

    def __hash__(self):
        return hash(self.name)

    def __len__(self):
        return len(self.seq)

    def __eq__(x,y):
        return x.name == y.name

    def __getitem__(self,index):
        start, end = self.start, self.end
        # print self.start,end,index
        if index < self.start or index > end:
            return -1
        else:
            return self.seq[index-self.start]

    def get_coverage(self,index):
        start, end = self.start, self.end
        # print self.start,end,index
        if index < start or index > end:
            return 0
        else:
            return self.coverage[index-start]

    def get_prob(self,index):
        start, end = self.start, self.end
        if index < start or index > end:
            return dict()
        else:
            return self.accuracies[index-start]

    def set_coverage(self,index,n):
        if index >= self.start and index <= self.end:
            self.coverage[index - start] = n

    def from_str(self, sequence, coverage):
        """ Builds cloud form a sequence string, like 0011--01-1,
            and a coverage string, like 3$52--*?-4@ """

        for j,(s,c) in enumerate(zip(sequence,coverage)):
            if s == "-":
                self.seq.append(-1)
            else:
                self.seq.append(int(s))
                if int(s) == 1:
                    self.one_positions.add(j+self.start)
                else:
                    self.zero_positions.add(j+self.start)
            if c == "-":
                self.coverage.append(0)
            else:
                self.coverage.append(ord(c)-48)

        # positions covered by the segment are the union of 0,1 positions:
        self.positions = self.one_positions | self.zero_positions


    def from_bed(self, covered_positions, variants_at_pos, coverages_at_pos):
        seq = list()
        for j in xrange(cloud_start, cloud_end+1):
            if j in covered_positions:
                seq.append(variants_at_pos[j])
                coverage_str = coverages_at_pos[j]
                phred_score = (reduce(lambda x,y : x+y,
                                      [ord(x)-33 for x in coverage_str]))
                prob = 10**(-float(phred_score)/10.0)
                self.accuracy.append(prob)
            else:
                seq.append(-1)
                self.accuracy.append(-1)

        self.seq = seq
        cloud.name = ""

        # place it into our sets:
        for j in xrange(cloud_start, cloud_end):
            if cloud[j] != -1:
                clouds_at_position[cloud_start+k].add(cloud)

    def from_01_bed(self, covered_positions, alleles, qscores):
        start, end = self.start, self.end
        seq = list()
        accuracies = list()

        for j in xrange(start, end+1):
            if j in covered_positions and qscores[j][0] > 5:
                seq.append(alleles[j][0])
                if alleles[j][0] == 1:
                    self.one_positions.add(j)
                else:
                    self.zero_positions.add(j)

                p = dict()
                for a, qscore in zip(alleles[j], qscores[j]):
                    p[a] = 10**(-1*qscore/10)
                assert len(p.keys()) > 0
                accuracies.append(p)
            else:
                seq.append(-1)
                accuracies.append(dict())

        # if len(self.one_positions) + len(self.zero_positions) == 0: keyboard()
        self.seq = seq
        self.accuracies = accuracies

        # positions covered by the segment are the union of 0,1 positions:
        self.positions = self.one_positions | self.zero_positions

    def to_str(self):
        """ Returns a string with the entries and coverage of the cloud,
            for example: 0011--01-1 3$52--*?-4@. """

        subst = lambda x: x if x != "-1" else "-"
        string = "".join(subst(str(x)) for x in self.seq)
        string += " "

        for allele, accuracy in zip(self.seq, self.accuracies):
            if allele == -1: 
                string += "-"
            else: 
                qscore = int(-round(log(accuracy[allele], 10))*10)
                string += chr(min(33+qscore,126))

        return string

    def matches(self,c):
        """ Returns true if likelihood of not matching is < 1e-12 """

        threshold = 1e-12
        log_threshold = log(threshold)
        start = min(self.start, c.start)
        end = max(self.start + len(self) - 1, c.start + len(c) - 1)

        # Return false if segment doesn't match:
        # for j in xrange(start,end):
            # if self[j] != c[j] and self[j] != -1 and c[j] != -1: return False

        if not any((self[j] != -1 and c[j] != -1 
                    for j in xrange(start,end+1))): 
            return False

        equal_positions = [c[j] == self[j] for j in xrange(start, end+1)]
        if all(equal_positions): return True

        try:
            log_mismatch_prob = (reduce(lambda x,y: x+y,
                [logaddexp(min(self.lik(j,1)+c.lik(j,0), -1e-15),
                           min(self.lik(j,0)+c.lik(j,1), -1e-15))
                for j in xrange(start,end+1)
                if self[j] != -1 and c[j] != -1]))
            log_match_prob = (reduce(lambda x,y: x+y,
                [logaddexp(min(self.lik(j,0)+c.lik(j,0), -1e-15),
                    min(self.lik(j,1)+c.lik(j,1), -1e-15))
                 for j in xrange(start,end+1)
                 if self[j] != -1 and c[j] != -1]))
        except TypeError:
            keyboard()

        if log_mismatch_prob - logaddexp(log_match_prob, log_mismatch_prob) < \
           log_threshold:
            return True
        else:
            return False

    def prob(self, j, v):
        """ Returns prob(cloud[j]|v) """
        start = self.start
        cloud_v = self.seq[j-start]

        if self[j] == -1:
            return 0.0
        else:
            if len(self.accuracies[j-start]) == 1:
                p = self.accuracies[j-start][cloud_v]
                if cloud_v == v: return (1.0-p)
                else: return p
            elif len(self.accuracies[j-start]) == 2:
                p0 = self.accuracies[j-start][0]
                p1 = self.accuracies[j-start][1]

                if v == 1:
                    return (1-p1)*p0
                elif v == 0:
                    return (1-p0)*p1
                else:
                    print "Error!"
                    exit(-1)
            else:
                print "Error!"
                exit()

    def lik(self,j,v):
        """ Returns log prob(cloud[j]|v) """
        start = self.start
        cloud_v = self.seq[j-start]

        if self[j] == -1:
            return float("-inf")
        else:
            if len(self.accuracies[j-start]) == 1:
                p = self.accuracies[j-start][cloud_v]
                if cloud_v == v: return log(max(1.0-p, 1e-14))
                else:
                    if p == 0: return -1e20
                    else: return log(p)
            elif len(self.accuracies[j-start]) == 2:
                p0 = self.accuracies[j-start][0]
                p1 = self.accuracies[j-start][1]

                if v == 1:
                    p = (1-p1)*p0
                elif v == 0:
                    p = (1-p0)*p1

                if p == 0: return -1e20
                else: return log(p)
            else:
                print "Error!!"
                exit(1)
