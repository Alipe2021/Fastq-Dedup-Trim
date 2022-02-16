#!/usr/bin/env python
'''
Author: Liu Peng (sxliulian2012@hotmail.com)
Date: 2022-02-15 15:41:18
LastEditTime: 2022-02-16 16:19:58
LastEditors: Please set LastEditors
Description: Convert, Deduplication and trimming reads.
FilePath: Scripts/fastq-dedup-trim.py
Version: 0.2
'''
from collections import defaultdict
from functools import reduce
from pysam import AlignmentFile
import argparse
import hashlib
import gzip
import os

TrieNode = lambda: defaultdict(TrieNode)
class Trie:
    """ Build trie """
    def __init__(self):
        self.trie = TrieNode()
    def insert(self, word):
        reduce(dict.__getitem__, word, self.trie)['end'] = True
    def search(self, word):
        return reduce(lambda d,k: d[k] if k in d else TrieNode(), word, self.trie).get('end', False)

def phred_based_trimming(seq, qual, keep_length:int):
    """ Timming fastq based phred quality. \n
        Return trimmed sequence and quality.
    """
    assert len(seq) >= keep_length
    if len(seq) == keep_length:
        return seq, qual
    else:
        ThisQual, MaxQual, n = 0, 0, 0
        seq_qual = [ord(q) for q in qual]
        for i in range(len(seq)-keep_length):
            ThisQual = sum([q for q in seq_qual[i:i+keep_length]])
            if ThisQual > MaxQual:
                MaxQual = ThisQual
                n = i
        return seq[n:n+keep_length], qual[n:n+keep_length]

def read_pair_generator(bam, region_string=None):
    """ 
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    from collections import defaultdict
    reads_dict = defaultdict(lambda: [None, None])

    for read in bam.fetch(until_eof=True, region=region_string):  # for unmapped reads
        if not read.is_paired:
            continue

        qname = read.qname
        if qname not in reads_dict:
            if read.is_read1:
                reads_dict[qname][0] = read
            else:
                reads_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, reads_dict[qname][1]
            else:
                yield reads_dict[qname][0], read
            del reads_dict[qname]
# 
def write_read_to_fastq(fastq, read):
    """
    Write read to open FASTQ file.
    """
    info = {'index': int(not read.is_read1) + 1,
            'name':  read.qname}

    if read.is_reverse:
        info.update({'quality':  read.qual[::-1],
                     'sequence': reverse_complement2(read.seq)})
    else:
        info.update({'quality':  read.qual,
                     'sequence': read.seq})
    fastq.write('@{name}/{index}\n{sequence}\n+\n{quality}\n'.format(**info))

def write_fastq_to_dict(dict, read):
    """
    Write read to open FASTQ file.
    """
    out_dict = dict()
    info = {'index': int(not read.is_read1) + 1,
            'name':  read.qname}

    if read.is_reverse:
        info.update({'quality':  read.qual[::-1],
                     'sequence': reverse_complement2(read.seq)})
    else:
        info.update({'quality':  read.qual,
                     'sequence': read.seq})
    out_dict.update({read.qname : '@{name}/{index}\n{sequence}\n+\n{quality}'.format(**info)})

def reverse_complement1(sequence):
    """
    Return reverse complement of DNA sequence.
    """
    # DNA base complements
    COMPLEMENT = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c'
    }
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])

def reverse_complement2(sequence):
    """
    Return reverse complement of DNA sequence.
    """
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string[::-1]

def string_to_md5(string):
    """ Check md5 of string. """
    m = hashlib.md5()
    m.update(string.encode('utf-8'))
    return m.hexdigest()

def main():
    parser = argparse.ArgumentParser(description='Deduplication and Trimming fastq from bam file.')
    parser.add_argument('bam', help='Input bam file.')
    parser.add_argument('-p', '--prefix', help='Prefix of output file.')
    # parser.add_argument('--deduplication', help='Deduplication.', action="store_true")
    parser.add_argument('--keep-length', help='Sequence length for keeping.', type=int, default=145)
    # parser.add_argument('--out-format', help='Format to output.', type=str, default="bam")
    args = parser.parse_args()

    outdir = os.path.dirname(args.prefix)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    id_md5_set = set()
    fastq_r1 = gzip.open(args.prefix + "_DedupTrim_R1.fq.gz", "wt")
    fastq_r2 = gzip.open(args.prefix + "_DedupTrim_R2.fq.gz", "wt")

    with AlignmentFile(args.bam, 'rb') as in_bam:
        for read1, read2 in read_pair_generator(in_bam):
            if len(read1.seq) < args.keep_length or len(read2.seq) < args.keep_length:
                continue
            else:
                read1.seq, read1.qual = phred_based_trimming(read1.seq, read1.qual, keep_length=args.keep_length)
                read2.seq, read2.qual = phred_based_trimming(read1.seq, read1.qual, keep_length=args.keep_length)
                
                seq_md5 = string_to_md5("_".join([read1.seq, read2.seq]))
                if seq_md5 in id_md5_set:
                    continue
                else:
                    id_md5_set.add(seq_md5)
                    write_read_to_fastq(fastq_r1, read1)
                    write_read_to_fastq(fastq_r2, read2)

    fastq_r1.close()
    fastq_r2.close()

if __name__ == "__main__":
    main()