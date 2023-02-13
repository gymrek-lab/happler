#!/usr/bin/env python

import sys
import cyvcf2
import argparse
import trtools.utils.tr_harmonizer as trh


def parse_args():
    parser = argparse.ArgumentParser(
        description=
        "Create an STR genotype matrix with the sum of allele differences from the REF."
    )
    parser.add_argument(
        '-o', '--out', default=sys.stdout,
        help='path to TSV file where variants are rows and samples are cols'
    )
    parser.add_argument(
        '-i', '--id', default=None,
        help='the ID of a specific STR; ignores all others'
    )
    parser.add_argument(
        'vcf', nargs='?', default=sys.stdin,
        help='a VCF containing STRs'
    )
    args = parser.parse_args()
    return args

def main(args):
    invcf = cyvcf2.VCF(args.vcf)
    harmonizer = trh.TRRecordHarmonizer(invcf)
    for trrecord in harmonizer:
        print(trrecord.GetLengthGenotypes(), file=args.out)

if __name__ == '__main__':
    main(parse_args())
