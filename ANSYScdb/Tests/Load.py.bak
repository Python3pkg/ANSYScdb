# -*- coding: utf-8 -*-
"""
Tests/Examples to load cdb files.
"""
from __future__ import print_function

from os.path import dirname, join, realpath
from ANSYScdb import CDB_Reader


def HexBeam():
    """ Load and plot hexahedral beam """
    
    # get location of this file
    pth = dirname(realpath(__file__))
    filename = join(pth, 'HexBeam.cdb')
    
    # Create cdb object
    cdb = CDB_Reader.Read(filename)

    print('Raw cdb loaded with the following keys:')
    for key in cdb.raw:
        print(key, end='')
        print(':\tWith shape: ', end='')
        try:
            print(cdb.raw[key].shape)
        except:
            print('')

    return cdb
    
    
def TetBeam():
    """ Load and plot tetrahedral beam """
    
    # get location of this file
    pth = dirname(realpath(__file__))
    filename = join(pth, 'TetBeam.cdb')
    
    # Create cdb object
    cdb = CDB_Reader.Read(filename)

    print('Raw cdb loaded with the following keys:')
    for key in cdb.raw:
        print(key, end='')
        print(':\tWith shape: ', end='')
        try:
            print(cdb.raw[key].shape)
        except:
            print('')

    return cdb
