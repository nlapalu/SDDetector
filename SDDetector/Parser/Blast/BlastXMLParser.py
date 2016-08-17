#!/usr/bin/env python

import logging

from Bio.Blast import NCBIXML
from SDDetector.Entities.Alignment import Alignment


class BlastXMLParser(object):

    def __init__(self, inputBlastXMLFile=""):
        """Constructor"""
        self.inputBlastXMLFile = inputBlastXMLFile

    def getAllAlignments(self):
        """Return list of all Alignments"""

        lAlignments = []
        with open(self.inputBlastXMLFile,  'r') as input :

            blast_records = NCBIXML.parse(input)
            index = 0

            for blast_record in blast_records:
                logging.debug('QUERY: {}'.format(blast_record.query))

                for alignment in blast_record.alignments:
                    logging.debug('SUBJECT: {}'.format(alignment.hit_id))
                    nb_hsp = 0
                    for hsp in alignment.hsps:
                        nb_hsp += 1
                        index += 1
                        if hsp.frame[1] == 1:
                            lAlignments.append(Alignment(blast_record.query, alignment.hit_id,
                                               hsp.query_start, hsp.query_end, hsp.sbjct_start,
                                               hsp.sbjct_end, hsp.align_length, hsp.identities, hsp.frame[0], hsp.frame[1], id=index))
                        elif hsp.frame[1] == -1:
                            lAlignments.append(Alignment(blast_record.query, alignment.hit_id,
                                               hsp.query_start, hsp.query_end, hsp.sbjct_end,
                                               hsp.sbjct_start, hsp.align_length, hsp.identities, hsp.frame[0], hsp.frame[1], id=index))
                        else:
                            logging.error('Blast Parsing: Unknown strand')
                            raise Exception("Unknown strand")
                    logging.debug('{} HSP parsed'.format(nb_hsp))
        input.closed
        return lAlignments


    def getAlignmentsFromTupleOfRegions(self,lRegions):
        """Return list of duplications with alignment"""

        lindex = []
        for reg1,reg2 in lRegions:
            if reg1.strand == -1:
                lindex.append((reg1.seq,reg1.end,reg1.start,reg2.seq,reg2.start,reg2.end))
            else:
                lindex.append((reg1.seq,reg1.start,reg1.end,reg2.seq,reg2.start,reg2.end))
#        print "INDEX {}".format(lindex)

        lRegionAlgmts = [ () for i in range(len(lRegions)) ]

        with open(self.inputBlastXMLFile,  'r') as input :

            blast_records = NCBIXML.parse(input)
            for blast_record in blast_records:
                logging.debug('QUERY: {}'.format(blast_record.query))
                for alignment in blast_record.alignments:
                    logging.debug('SUBJECT: {}'.format(alignment.hit_id))
                    for hsp in alignment.hsps:
                        logging.debug('{}-{}-{}-{}-{}-{}'.format(alignment.hit_id,hsp.sbjct_start,hsp.sbjct_end,blast_record.query,hsp.query_start,hsp.query_end))
                        if (alignment.hit_id,hsp.sbjct_start,hsp.sbjct_end,blast_record.query,hsp.query_start,hsp.query_end) in lindex:
                            index = lindex.index((alignment.hit_id,hsp.sbjct_start,hsp.sbjct_end,blast_record.query,hsp.query_start,hsp.query_end))
                            lRegionAlgmts[index] = (hsp.sbjct,hsp.query)

        input.close()
        return lRegionAlgmts


if __name__ == "__main__":
    blastXMLParser = BlastXMLParser()
