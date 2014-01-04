
import os
from escherpy.escherDB import *
from nose.tools import assert_equals


def test_create_db():
    create_db()
    assert_equals(os.path.isfile('escher.db'), True)



def test_getat():

    setat('rcov', 3.1, Z=7)
    setat('weight', 40., name='Na')

    delat(Z=1)
    newcolat('rvdw', 'real')

    printat()
    assert_equals( getat('rcov', Z=7) , 3.1)
    assert_equals( getat('weight', name='Na') , 40.0)
    os.remove('escher.db')
