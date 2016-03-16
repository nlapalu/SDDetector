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

if __name__ == "__main__":
    blastXMLParser = BlastXMLParser()
