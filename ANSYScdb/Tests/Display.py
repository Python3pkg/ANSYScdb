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

    # Parse into unstuctured grid
    cdb.ParseVTK()
    
    print(cdb.uGrid)

    # Plot unstuctured grid
    cdb.Plot()
    
    # Return cdb object for use to inspect
    return cdb
    
    
def TetBeam():
    """ Load and plot tetrahedral beam """
    
    # get location of this file
    pth = dirname(realpath(__file__))
    filename = join(pth, 'TetBeam.cdb')
    
    # Create cdb object
    cdb = CDB_Reader.Read(filename)

    # Parse into unstuctured grid
    cdb.ParseVTK()
    print(cdb.uGrid)

    # Plot unstuctured grid
    cdb.Plot()
    
    # Return cdb object for use to inspect
    return cdb
