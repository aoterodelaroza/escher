#!/usr/bin/env python

import sqlite3

from atom import t


#db = sqlite3.connect( 'escher.db')
#cur = db.cursor()

def create_db():

    db = sqlite3.connect( 'escher.db')
    cur = db.cursor()

    # the id autoincrement is performed giving a null value
    # the assigned id value is lastrowid + 1
    cur.execute("create table atoms (Z integer primary key autoincrement, name, weight, rcov, Rjmolcolor, Gjmolcolor, Bjmolcolor)")

    #t = []
    #for key in colors.keys():
    #    t.append((key, 0.0, rcov[key], colors[key][0], colors[key][1], colors[key][2]))
    #print t

    cur.executemany("insert into atoms values (null, ?, ?, ?, ?, ?, ?)", t)

    cur.close()
    db.commit()
    db.close()

#create_db()

#t = (12, 4)
#cur.execute("update atoms set rcov=? where Z=?", t)

def setat(setfield, setvalue,  **kwargs):
    '''
    Set rcov=3.1 in the row that corresponds to Hydrogen
    setat('rcov', 3.1, Z=1)
    setat('weight', 40., name='Na')
    '''
    db = sqlite3.connect( 'escher.db')
    cur = db.cursor()

    sql = "update atoms set {}=? where {}=?".format(setfield, kwargs.keys()[0])
    cur.execute(sql, (setvalue, kwargs.values()[0],))

    cur.close()
    db.commit()
    db.close()


#setat('rcov', 3.1, Z=7)
#setat('weight', 40., name='Na')


def delat(**kwargs):

    db = sqlite3.connect( 'escher.db')
    cur = db.cursor()

    sql = "delete from atoms where {}=?".format(kwargs.keys()[0])
    cur.execute(sql, (kwargs.values()[0],))

    cur.close()
    db.commit()
    db.close()

def newcolat(colname, coltype):

    db = sqlite3.connect( 'escher.db')
    cur = db.cursor()

    sql = "alter table atoms add column '{}' {}".format(colname, coltype)
    cur.execute(sql)

    cur.close()
    db.commit()
    db.close()

def printat():

    db = sqlite3.connect( 'escher.db')
    cur = db.cursor()

    cur.execute('select * from atoms')
    for atom in cur.fetchall():

        print(atom)

    cur.close()
    db.commit()
    db.close()

#delat(Z=1)
#cur.execute("delete from atoms where Z=1")
#cur.execute("alter table atoms add column 'rvdw' real")

#cur.execute("select * from atoms")

#for atom in cur.fetchall():

#    print(atom)



def getat(get='rcov', **kwargs):
    '''
    getat('rcov', Z=1)
    getat('rcov', name='Na')
    '''
    db = sqlite3.connect( 'escher.db')
    cur = db.cursor()

    sql = "select {} from atoms where {}=?".format(get, kwargs.keys()[0])
    cur.execute(sql, (kwargs.values()[0],))

    atom, = cur.fetchone()

    cur.close()
    db.commit()
    db.close()

    return atom

#getat('rcov', Z=7)
#getat('weight', name='Na')

#for i in range(2,100):
#    print i
#    getat('rcov', Z=i)


#cur.close()
#db.commit()
#db.close()
