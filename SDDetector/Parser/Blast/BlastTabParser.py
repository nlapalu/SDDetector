#!/usr/bin/env python

import logging

from SDDetector.Entities.Alignment import Alignment


class BlastTabParser(object):

    def __init__(self, inputBlastTabFile=""):
        """Constructor"""
        self.inputBlastTabFile = inputBlastTabFile

    def getAllAlignments(self):
        """Return list of all Alignments"""

        lAlignments = []
        with open(self.inputBlastTabFile,  'r') as input :
            index = 0
            nb_hsp = 0
            for line in input:
                index += 1
                nb_hsp += 1
                qseqid,sseqid,qstart,qend,sstart,send,length,nident = line.split('\t')
                qframe = 1
                sframe = None
                if int(sstart) < int(send):
                    sframe = 1
                else:
                    sframe = -1

                if sframe == 1: 
                    lAlignments.append(Alignment(qseqid, sseqid, int(qstart), int(qend), int(sstart), int(send), int(length),
                                                 int(nident), int(qframe), int(sframe), id=index))
                if sframe == -1:
                    lAlignments.append(Alignment(qseqid, sseqid, int(qstart), int(qend), int(send), int(sstart), int(length),
                                                 int(nident), int(qframe), int(sframe), id=index))

            logging.debug('{} HSP parsed'.format(nb_hsp))

        input.closed
        return lAlignments

if __name__ == "__main__":
    blastTabParser = BlastTabParser()
