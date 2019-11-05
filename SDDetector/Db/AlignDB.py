#!/usr/bin/env python

import os

from SDDetector.Db.SqliteDB import SqliteDB
from SDDetector.Entities.Alignment import Alignment


class AlignDB(SqliteDB):

    def __init__(self, dbfile='', logLevel='ERROR', copy=False):
        """AlignDB Constructor"""

        SqliteDB.__init__(self, dbfile=dbfile, logLevel=logLevel, copy=copy)
        if copy == False:
            self._createDBSchema()

    def _createDBSchema(self):
        """Create Database Schema"""

        self._createAlignmentTable()

    def _createAlignmentTable(self):
        """Create table Alignment"""

        self.conn.execute('''create table alignment (id int primary key not null,
                                                     query text not null,
                                                     sbjct text not null,
                                                     qstart int not null,
                                                     qend int not null,
                                                     sstart int not null,
                                                     send int not null,
                                                     length int not null,
                                                     identities int not null,
                                                     qstrand int not null,
                                                     sstrand int not null);''')
        self.conn.commit()

    def insertAlignment(self, algmt):
        """Insert an alignment in db"""

        self.conn.execute('''insert into alignment (id, query, sbjct,
                             qstart, qend, sstart, send, length, identities,
                             qstrand,sstrand) values (%s)''' % (','.join('\'{}\''.format(i) for i in (algmt.id,algmt.query,algmt.sbjct, algmt.qstart, algmt.qend, algmt.sstart, algmt.send, algmt.length, algmt.identities, algmt.qstrand, algmt.sstrand))))


    def insertlAlignments(self, lAlignments, length=0):
        """Insert a list of alignments"""

        if length > 0:
            lAlignments = [ algmt for algmt in lAlignments if algmt.length > length ]

        self.conn.executemany('''insert into alignment (id, query, sbjct,
                                 qstart, qend, sstart,send, length, identities,
                                 qstrand, sstrand) values (?,?,?,?,?,?,?,?,?,?,?)''', [(algmt.id, algmt.query,algmt.sbjct, algmt.qstart, algmt.qend, algmt.sstart, algmt.send, algmt.length, algmt.identities, algmt.qstrand, algmt.sstrand) for algmt in lAlignments ])
        self.commit()


    def deleteAlignment(self, id):
        """Delete an alignment"""

        cursor = self.conn.execute('''delete from alignment where id = {}''' \
                                   .format(id))


    def deletelAlignments(self, lIds):
        """Delete a list of alignments"""

        cursor = self.conn.execute('''delete from alignment where id in ({})''' \
                                   .format(','.join(str(id) for id in lIds)))


    def deleteAllAlignments(self):
        """Delete all alignments"""

        cursor = self.conn.execute('''delete from alignment''')
        self.commit()


    def deleteDB(self):
        """Delete the DB"""

        self.conn.close()
        os.remove(self.dbfile)


    def exportDbToGff3(self, fileName):
        """Export entire db in gff3 match features"""

        with open(fileName,'w') as f:
            cursor = self.conn.execute('''select id, query, sbjct, qstart,
                                          qend, sstart, send, length,identities,
                                          qstrand, sstrand from alignment''')

            for row in cursor :
                f.write(Alignment(row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], id=row[0]).convertToGff3())

        f.close()


    def selectAlignmentMaxId(self):
        """Select alignment max id"""

        cursor = self.conn.execute('''select max(id) from alignment''')
        return cursor.fetchone()[0]


    def selectAllIds(self):
        """Select all ids alignments"""

        cursor = self.conn.execute('''select id from alignment''')
        return [ row[0] for row in cursor ]


    def selectAllAlignments(self):
        """Select all alignments"""

        cursor = self.conn.execute('''select id, query, sbjct, qstart,
                                   qend, sstart, send, length, identities,
                                   qstrand, sstrand from alignment''')
        return [ Alignment(row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], id=row[0]) for row in cursor ]


    def selectAlignmentById(self, id):
        """Select alignment by id"""

        cursor = self.conn.execute('''select id, query, sbjct, qstart,
                                   qend, sstart, send, length, identities,
                                   qstrand, sstrand from alignment where id = {}''' \
                                   .format(id))
        row = cursor.fetchone()
        if row:
            return Alignment(row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], id=row[0])
        else:
            return None

    def selectAllSbjcts(self):
        """Select all sequence name used as subject"""

        cursor = self.conn.execute('''select distinct(sbjct) from alignment order by sbjct ASC''' )

        return [ row[0] for row in cursor ]

    def selectAllQueries(self):
        """Select all sequence name used as query"""

        cursor = self.conn.execute('''select distinct(query) from alignment order by query ASC''' )

        return [ row[0] for row in cursor ]

    def selectAlignmentsWithDefinedSbjctAndQueryOrderBySbjctCoord(self, sbjct, query):
        """Select ordered sub-list of alignments with defined sbjct and query"""

        cursor = self.conn.execute('''select id, query, sbjct, qstart,
                                   qend, sstart, send, length, identities,
                                   qstrand, sstrand from alignment where \
                                   sbjct = \'{}\' and query = \'{}\' \
                                   order by sstart ASC''' \
                                   .format(sbjct,query))

        return [ Alignment(row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], id=row[0]) for row in cursor ]


    def selectAlignmentsWithDefinedIdOrderBySbjctCoord(self, lIds):
        """Select ordered sub-list of alignments with defined id"""

        cursor = self.conn.execute('''select id, query, sbjct, qstart,
                                   qend, sstart, send, length, identities,
                                   qstrand, sstrand from alignment where \
                                   id in ({}) order by sstart ASC''' \
                                   .format(','.join([ '\'{}\''.format(id) for id in lIds ])))

        return [ Alignment(row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], id=row[0]) for row in cursor ]


    def selectSelfSelfMatchAlgmts(self):
        """Select all self-matches alignments"""

        cursor = self.conn.execute('''select id from alignment where query = sbjct \
                                   and qstart = sstart and qend = send ''')

        return [ row[0] for row in cursor ]


    def selectAlgmtsBelowIdentityThreshold(self, threshold=0.9):
        """Select alignments with identity below the defined threshold"""

        cursor = self.conn.execute('''select id, identities, qstart, qend, sstart, send
                                      from alignment ''')

        return [ row[0] for row in cursor if (float(row[1])/(row[3]-row[2]+1)) < threshold ]


    def selectProximalAlgmts(self, id, maxGap=3000):
        """
            Select alignments with a maximal distance of 'maxGap' between the provided
            alignment and alignments on identical (subject/query) further up the subject
            and collinear !!!!
        """

        cursor = self.conn.execute('''select al2.id, al2.query, al2.sbjct, al2.qstart,
                                   al2.qend, al2.sstart, al2.send, al2.length, al2.identities,
                                   al2.qstrand, al2.sstrand from alignment al1, alignment al2 where \
                                   (al2.sstart-al1.send) < {} and al1.sstart < al2.sstart \
                                   and (al2.sstart-al1.send) > -1 \
                                   and al1.sbjct = al2.sbjct and al1.query = al2.query \
                                   and al1.qstrand = al2.qstrand and al1.sstrand = al2.sstrand \
                                   and al1.id = {} order by al2.sstart''' \
                                   .format(maxGap,id))

        return [ Alignment(row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], id=row[0]) for row in cursor ]

    def selectSuboptimalAlgmts(self, maxOverlap):
        """
            A suboptimal alignment is an alignment in the sense that the subject and query are
            complety covered or spanned by another alignment.
        """
        cursor = self.conn.execute('''select distinct al2.id from alignment al1, alignment al2 \
                                   where \
                                   ((al1.sstart >= al2.sstart and al1.sstart < al2.send and (al2.send-al1.sstart) >= {} ) \
                                   or (al1.send > al2.sstart and al1.send <= al2.send and (al1.send-al2.sstart) >= {} )) \
                                   and ((al1.qstart >= al2.qstart and al1.qstart < al2.qend and (al2.qend-al1.qstart) >= {} )   \
                                   or (al1.qend > al2.qstart and al1.qend <= al2.qend and (al1.qend-al2.qstart) >= {})) \
                                   and  al1.sbjct = al2.sbjct and al1.query = al2.query \
                                   and al1.length > al2.length'''.format(maxOverlap, maxOverlap,maxOverlap, maxOverlap))

        return [ row[0] for row in cursor ]


    def selectQueryOnlySuboptimalAlgmts(self, lIds):

        cursor = self.conn.execute('''select distinct al2.id from alignment al1, alignment al2 \
                                   where \
                                   ((al1.qstart >= al2.qstart and al1.qstart < al2.qend) \
                                   or (al1.qstart <= al2.qstart and al1.qend > al2.qstart) \
                                   or (al1.qend > al2.qstart and al1.qend <= al2.qend)) \
                                   and al1.query = al2.query \
                                   and al1.length >= al2.length \
                                   and al1.id != al2.id \
                                   and al1.id in ({})\
                                   and al2.id in ({})'''.format( \
                                   ','.join(str(id) for id in lIds),\
                                   ','.join(str(id) for id in lIds)))

        return [ row[0] for row in cursor ]


    def selectSubjectOnlySuboptimalAlgmts(self, lIds):

        cursor = self.conn.execute('''select distinct al2.id from alignment al1, alignment al2 \
                                   where \
                                   ((al1.sstart >= al2.sstart and al1.sstart < al2.send) \
                                   or (al1.sstart <= al2.sstart and al1.send > al2.sstart) \
                                   or (al1.send > al2.sstart and al1.send <= al2.send)) \
                                   and al1.query = al2.query \
                                   and al1.length >= al2.length \
                                   and al1.id != al2.id \
                                   and al1.id in ({})\
                                   and al2.id in ({})'''.format( \
                                   ','.join(str(id) for id in lIds),\
                                   ','.join(str(id) for id in lIds)))

        return [ row[0] for row in cursor ]


    def  commit(self):
        """Commit transactions"""

        self.conn.commit()
