#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict

class BEDPE (object):
    def __init__(self, bedpeLine):
        bedpeFields = bedpeLine.split('\t')
        self.ref1    = bedpeFields[0]
        self.left1   = int(bedpeFields[1])
        self.right1  = int(bedpeFields[2])
        self.ref2    = bedpeFields[3]
        self.left2   = int(bedpeFields[4])
        self.right2  = int(bedpeFields[5])
        self.name    = bedpeFields[6]
        self.score   = bedpeFields[7]
        self.strand1 = bedpeFields[8]
        self.strand2 = bedpeFields[9]
        if (self.left1 > self.right1 or self.left2 > self.right2):
            raise Exception('BEDPE with left > right')
        

class BED (object):
    def unparse(self):
        return '\t'.join([self.ref, str(self.left), str(self.right), self.name, self.score, self.strand])
        
    def __init__(self, ref, left, right, name, score, strand):
        self.ref    = ref 
        self.left   = left
        self.right  = right
        self.name   = name
        self.score  = score
        self.strand = strand
        if (self.left > self.right):
            raise Exception('BED with left > right')
        
def bedpeToBed(bedpe, stats, maxsize):
    if (bedpe.ref1 == '.' and bedpe.ref2 == '.'):
        stats['nomatch'] += 1
    elif (bedpe.ref1 == '.' or bedpe.ref2 == '.'):
        stats['onematch'] += 1
    elif (bedpe.ref1 != bedpe.ref2):
        stats['diffchr'] += 1
    elif (bedpe.strand1 == bedpe.strand2):
        stats['samestrand'] += 1
    elif (bedpe.strand1 == '+' and bedpe.strand2 == '-'):        
        if (bedpe.left1 + maxsize <= bedpe.right2):
            stats['toolong'] += 1
        elif (bedpe.left1 <= bedpe.right2):
            stats['good'] += 1
            return BED(bedpe.ref1, bedpe.left1, bedpe.right2, bedpe.name, bedpe.score, bedpe.strand1)
        else:            
            stats['inverted'] += 1
    elif (bedpe.strand1 == '-' and bedpe.strand2 == '+'):
        if (bedpe.left2 + maxsize <= bedpe.right1):
            stats['toolong'] += 1
        elif (bedpe.left2 <= bedpe.right1):
            stats['good'] += 1
            return BED(bedpe.ref1, bedpe.left2, bedpe.right1, bedpe.name, bedpe.score, bedpe.strand1)
        else:
            stats['inverted'] += 1
    else:
        raise Exception('BEDPE confusing ' + bedpe.ref1 + ', ' + bedpe.ref2 + ', ' + bedpe.strand1 + ', ' + bedpe.strand2)
    return None
    
argparser = argparse.ArgumentParser(description='Convert BEDPE file to BED file')
argparser.add_argument('bedpe', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
argparser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
argparser.add_argument('-m', '--maxsize', type=int, default=1000000)
statgroup = argparser.add_mutually_exclusive_group()
statgroup.add_argument('-s', '--stats', type=argparse.FileType('w'), default=sys.stderr)
statgroup.add_argument('--no-stats', action="store_true")
args = argparser.parse_args()

stats = defaultdict(int)

for line in args.bedpe:
    try:
        bedpe = BEDPE(line.strip())
        bed = bedpeToBed(bedpe, stats, args.maxsize)
        if bed is not None:
            args.output.write(bed.unparse() + '\n')
    except Exception as e:
        print '** ERROR ', e, 'on line', line

if not args.no_stats:
    total = 0
    for fate in stats:
        total += stats[fate]
    for fate in stats:
        args.stats.write('{:<12} {:>9} {:.1%}\n'.format(fate, stats[fate], float(stats[fate]) / total))
    args.stats.write('{:<12} {:>9}\n'.format('TOTAL', total))
