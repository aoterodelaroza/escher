#!/usr/bin/env python

import sqlite3

from atom import t


db = sqlite3.connect( 'escher.db')

cur = db.cursor()

def create_db():

    # the id autoincrement is performed giving a null value
    # the assigned id value is lastrowid + 1
    cur.execute("create table atoms (Z integer primary key autoincrement, name, weight, rcov, Rjmolcolor, Gjmolcolor, Bjmolcolor)")

    #t = []
    #for key in colors.keys():
    #    t.append((key, 0.0, rcov[key], colors[key][0], colors[key][1], colors[key][2]))
    #print t

    cur.executemany("insert into atoms values (null, ?, ?, ?, ?, ?, ?)", t)

create_db()

t = (12, 4)
cur.execute("update atoms set rcov=? where Z=?", t)

def setat(setfield, setvalue,  **kwargs):
    '''
    Set rcov=3.1 in the row that corresponds to Hydrogen
    setat('rcov', 3.1, Z=1)
    setat('weight', 40., name='Na')
    '''
    sql = "update atoms set {}=? where {}=?".format(setfield, kwargs.keys()[0])
    cur.execute(sql, (setvalue, kwargs.values()[0],))


setat('rcov', 3.1, Z=7)
setat('weight', 40., name='Na')


def delat(**kwargs):

    sql = "delete from atoms where {}=?".format(kwargs.keys()[0])
    cur.execute(sql, (kwargs.values()[0],))


delat(Z=1)
#cur.execute("delete from atoms where Z=1")
cur.execute("alter table atoms add column 'rvdw' real")

cur.execute("select * from atoms")

for atom in cur.fetchall():

    print(atom)



def getat(get='rcov', **kwargs):
    '''
    getat('rcov', Z=1)
    getat('rcov', name='Na')
    '''
    sql = "select {} from atoms where {}=?".format(get, kwargs.keys()[0])
    cur.execute(sql, (kwargs.values()[0],))

    for atom in cur.fetchone():
        print atom

getat('rcov', Z=7)
getat('weight', name='Na')

#for i in range(2,100):
#    print i
#    getat('rcov', Z=i)


cur.close()
db.commit()
db.close()
