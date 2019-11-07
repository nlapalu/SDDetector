#!/usr/bin/env python

"""
   segmental_duplication_filter
"""

import argparse
import logging
import os
import sys
import re

from SDDetector.version import __version__
from SDDetector.Entities.Chain import Chain

class GFF3Reader(object):

    def __init__(self, fh, loglevel='ERROR'):

        self.fh = fh
        self.loglevel = loglevel
        logging.basicConfig(level=self.loglevel)
        self.nbFeatures = 0
        self.sReferences = set()

        try:
            self.filehandle = open(self.fh, 'r')
        except Exception as e:
            logging.error(e)
            sys.exit(1)

    def __del__(self):
        """..."""

        pass

    @staticmethod
    def _stringToDict(string):
        """convert field 9 from string to dict"""

        dAttributes = {}
        for att in string.split(";"):
            if att and att.split():
                tag,val = att.split("=")
                dAttributes[tag] = val.split(",")
        return dAttributes

    def read(self):
        """Iterator on feature"""

        currentLine = None
        for idx, line in enumerate(self.filehandle):
            currentLine = line.strip()
            if not currentLine:
                pass
            elif re.match('^#', currentLine):
                pass
            else:
                m = re.search(r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([+-.])\s+(\S+)\s+(\S.*)$", currentLine)
                if m == None:
                    raise Exception("Error GFF format line:{}".format(idx))

                id = GFF3Reader._getFeatureTagValue('ID',m.group(9))
                dAttributes = GFF3Reader._stringToDict(m.group(9))
                score = None
                try:
                    score = float(m.group(6))
                except:
                    pass
                strand = GFF3Reader._getStrand(m.group(7))
                frame = None
                if m.group(8).isdigit():
                    frame = int(m.group(8))
                self.sReferences.add(m.group(1))
                f = Feature(id,m.group(1),m.group(2),m.group(3),int(m.group(4)),int(m.group(5)),score,strand, frame, dAttributes)
                self.nbFeatures += 1

                yield f

    @staticmethod
    def convertRowToFeature(row, idx=None):

        m = re.search(r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([+-.])\s+(\S+)\s+(\S.*)$", row)
        if m == None:
            raise Exception("Error GFF format line:{}".format(idx))

        id = GFF3Reader._getFeatureTagValue('ID',m.group(9))
        dAttributes = GFF3Reader._stringToDict(m.group(9))
        score = None
        try:
            score = float(m.group(6))
        except:
            pass
        strand = GFF3Reader._getStrand(m.group(7))
        frame = None
        if m.group(8).isdigit():
            frame = int(m.group(8))
        f = Feature(id,m.group(1),m.group(2),m.group(3),int(m.group(4))+1,int(m.group(5)),score,strand, frame, dAttributes)
        return f



    @staticmethod
    def _getFeatureTagValue(tag, line):
        """Return the fist value of the tag property"""
        m = re.search(r".*{mytag}=([^;]*);{{0,1}}.*".format(mytag = tag),line)
        if m:
            return m.group(1).split(',')[0]
        else:
            raise Exception('Cannot find tag {} in string \'{}\''.format(tag, line))


    @staticmethod
    def _getStrand(strand):
        """Return strand as integer(1,-1) instead of +,- """

        if strand == '+':
            return 1
        elif strand == '-':
            return -1
        elif strand == '.':
            return None
        else:
            raise Exception('Cannot defined strand for feature')

    def getsReferences(self):
        """getter set of references"""

        return self.sReferences

        with open(self.fh, 'r') as input:
            for line in input:
                if not re.match('^#', line):
                    line = line.rstrip('\n')
                    values = line.split('\t')

    @staticmethod
    def _getFeatureTagValues(tag, line):
        """Return the list of values of the tag property"""
        m = re.search(r".*{mytag}=([^;]*);{{0,1}}.*".format(mytag = tag),line)
        if m:
            return m.group(1).split(',')
        else:
            raise Exception('Cannot find tag {} in string \'{}\''.format(tag, line))


class Feature(object):

    def __init__(self, id, seqid, source, type, start, end, score, strand, frame, attributes):
        """Feature constructor"""

        self.id = id
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = attributes

    def __eq__(self, other):
        """Equality on all args"""

        return ((self.id,self.seqid,self.source,self.type,self.start,self.end,self.score,self.strand,self.frame,self.attributes) == (other.id, other.seqid, other.source,other.type, other.start, other.end, other.score,other.strand, other.frame,other.attributes))

    def __repr__(self):
        """Feature representation"""

        return 'Feature: {}-{}-{}-{}-{}-{}-{}-{}-{}-{}'.format(self.type,self.id,self.seqid,self.source,self.start,self.end,self.score, self.strand,self.frame,self.attributes)




class Filter(object):

    def __init__(self, SDFile, AnnotFile, bed, outputFile, feature, size_annot, logLevel='ERROR'):

        self.SDFile = SDFile
        self.AnnotFile = AnnotFile
        self.bed = bed
        self.outputFile = outputFile
        self.feature = feature
        self.size_annot = size_annot
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel)


    def annot_coords(self):
        """extract coordinates"""

        logging.info("reading {}".format(self.AnnotFile))

        self.coords = []
        with open(self.AnnotFile, 'r') as f:
            for line in f:
                if not re.match('^#', line):
                    val = line.rstrip().split("\t")
                    if self.bed:
                        if int(val[2])-int(val[1]) > self.size_annot:
                            self.coords.append((val[0],int(val[1])+1,int(val[2])))
                    else:
                        if val[2] == self.feature and int(val[4])-int(val[3])+1 > self.size_annot:
                            self.coords.append((val[0],int(val[3]),int(val[4])))

        logging.info("{} regions extracted from {}".format(len(self.coords),self.AnnotFile))
        return self.coords


    def overlap_annot_sdd(self, fraction=0.1):
        """overlap between TE and specified feature"""

        logging.info("computing overlaps")

        lchainstoremove = set()
        gffr = GFF3Reader(self.SDFile)
        for feat in gffr.read():
            if feat.type == 'match_part':
                to_remove = False
                for co in self.coords:
                    nb_overlapping_bases = 0
                    if feat.seqid != co[0]:
                        continue
                    if feat.start <= co[1] <= feat.end <= co[2]:
                        nb_overlapping_bases += feat.end-co[1]+1
                    if co[1] <= feat.start <= co[2] <= feat.end:
                        nb_overlapping_bases += co[2]-feat.start+1
                    if feat.start <= co[1] <= co[2] <= feat.end:
                        nb_overlapping_bases += co[2]-co[1]+1
                    if co[1] <= feat.start and co[2] >= feat.end:
                        nb_overlapping_bases += feat.end-feat.start+1

                    if float(nb_overlapping_bases) / (feat.end-feat.start+1) > fraction:
                        to_remove = True
                        lchainstoremove.add(feat.attributes['Parent'][0][:-2])
                        continue
                if to_remove:
                    continue
        return lchainstoremove


    def export(self, lchains):

        if self.outputFile:
            o = open(self.outputFile,'w')

        with open(self.SDFile,'r') as f:
            for line in f:
                f = GFF3Reader.convertRowToFeature(line.rstrip())
                idx = ""
                if f.type == 'match':
                    idx = f.id[:-2]
                if f.type == 'match_part':
                    idx = f.attributes['Parent'][0][:-2]
                if idx not in lchains:
                    if self.outputFile:
                        o.write(line)
                    else:
                        print line.rstrip()


if __name__ == "__main__":

    program = 'SDDetector'
    version = __version__
    description = "SDDetector: detects segmental duplications in genome"

    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version', action='version', version='{} {}'.format(program,version))
    parser.add_argument("SDFile", help="gff3 file of segmental duplications", type=str)
    parser.add_argument("AnnotFile", help="gff3 file of segmental duplications", type=str)
    parser.add_argument("--bed", help="Bed file expected instead of GFF/GTF file", action="store_true", default=False)
    parser.add_argument("-o","--out", help="Output File in gff3 format", type=str, default=None)
    parser.add_argument("-f","--feature", help="Feature type expected for filtering default=match_part", type=str, default="match_part")
    parser.add_argument("-v", "--verbosity", type=int, choices=[1,2,3],
                        help="increase output verbosity 1=error, 2=info, 3=debug")
    parser.add_argument("-s", "--size", help="Minimum required size of feature for filtering, default=300bp", type=int, default=300)
    args = parser.parse_args()

    logLevel = 'ERROR'
    if args.verbosity == 1:
        logLevel = 'ERROR'
    if args.verbosity == 2:
        logLevel = 'INFO'
    if args.verbosity == 3:
        logLevel = 'DEBUG'
    logging.getLogger().setLevel(logLevel)


    ifilter = Filter(args.SDFile, args.AnnotFile, args.bed, args.out, args.feature, args.size, logLevel)

    ifilter.annot_coords()

    lchains = ifilter.overlap_annot_sdd()

    ifilter.export(lchains)
