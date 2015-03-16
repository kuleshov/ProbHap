#!/usr/bin/env python
import re
import argparse

parser = argparse.ArgumentParser(description="Combine .stats files")

parser.add_argument('--stat-files', required=True, nargs='+')
# parser.add_argument('--combined-stats', required=True)

args = parser.parse_args()

total_verifiable_pos = 0
total_het_pos = 0
total_phased_pos = 0
total_sw = 0
short_sw = 0
long_sw = 0
total_mb = 0.0

for sf in args.stat_files:
	with open(sf) as f:
		L = f.read()

		m = re.search('sw_acc: \d+.\d+ \d+ (\d+)\n', L)
		total_verifiable_pos += int(m.group(1))

		m = re.search('sw_per_mb: \d+.\d+ (\d+) (\d+.\d+)\n', L)
		total_sw += int(m.group(1))
		total_mb += float(m.group(2))

		m = re.search('long_sw_per_mb: \d+.\d+ (\d+) \d+.\d+\n', L)
		long_sw += int(m.group(1))

		m = re.search('short_sw_per_mb: \d+.\d+ (\d+) \d+.\d+\n', L)
		short_sw += int(m.group(1))

		m = re.search('percent_phased: \d+.\d+ (\d+) (\d+)\n', L)
		total_phased_pos += int(m.group(1))
		total_het_pos += int(m.group(2))


print "sw_acc:", (total_verifiable_pos - total_sw) / float(total_verifiable_pos), total_verifiable_pos - total_sw, total_verifiable_pos
print "long_sw_acc:", (total_verifiable_pos - long_sw) / float(total_verifiable_pos), total_verifiable_pos - long_sw, total_verifiable_pos
print "short_sw_acc:", (total_verifiable_pos - short_sw) / float(total_verifiable_pos), total_verifiable_pos - short_sw, total_verifiable_pos
print "sw_per_mb:", total_sw / total_mb, total_sw, total_mb
print "long_sw_per_mb:", long_sw / total_mb, long_sw, total_mb
print "short_sw_per_mb:", short_sw / total_mb, short_sw, total_mb
print "percent_phased:", total_phased_pos / float(total_het_pos), total_phased_pos, total_het_pos
