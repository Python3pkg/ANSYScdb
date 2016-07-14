"""
Module to read ANSYS ASCII block formatted CDB files

USAGE

# load module
from ANSYScdb import CDB

# load ANSYS cdb file
cdb = CDB.Reader('example.cdb')

# Parse the raw data into a VTK unstructured grid
cdb.ParseVTK()

# Plot the result
cdb.Plot()

"""
import warnings

# Attempt to load VTK dependent modules
try:
    import Utilities
    import Plotting
    vtk_failed = False

except:
    vtk_failed = True

import numpy as np

# Cython modules
import CDBparser
import CDBreader
import PythonReader
import PythonParser


class Read(object):
    """ FEM object """
    
    def __init__(self, filename, use_cython=True):
        """
        Initialize cdb object by reading raw cdb from file
        
        INPUTS:
        filename (string):
            filename of block formatted cdb file
            
        use_cython (bool optional):
            boolean flag to use cython reader defaults to True
    
        """
        
        # Defaults to cython reader if user selects it
        if use_cython:
            self.raw = CDBreader.Read(filename)

        # Python reader for debug purposes
        else:
            self.raw = PythonReader.Read(filename)



    def ParseVTK(self, use_cython=True):
        """
        Parses raw data from cdb file to VTK format
        
        """
        if vtk_failed:
            raise Exception('Unable to load VTK module.  Cannot parse raw cdb data!')
            return
           
        if self.CheckRaw():
            raise Exception('Missing key data.  Cannot parse into unstructured grid.')            
           
        # Convert to vtk style arrays
        cells, offset, cell_type, numref = CDBparser.Parse(self.raw)

        # Create unstructured grid
        self.uGrid = Utilities.MakeuGrid(offset, cells, cell_type,
                                         self.raw['nodes'][:, :3])

        # Store original ANSYS cell and node numbering
        Utilities.AddPointScalars(self.uGrid, self.raw['nnum'], 'ANSYSnodenum')

        # Add node components to unstructured grid
        ibool = np.empty(self.uGrid.GetNumberOfPoints(), dtype=np.int8)
        for comp in self.raw['node_comps']:
            ibool[:] = 0 # reset component array

            # Convert to new node numbering
            nodenum = numref[self.raw['node_comps'][comp]]
            
            ibool[nodenum] = 1
            Utilities.AddPointScalars(self.uGrid, ibool, comp)
            
        return self.uGrid
        
        
    def ParseFEM(self, use_cython=True):
        """ Parses raw data from cdb file to VTK format """
        if vtk_failed:
            raise Exception('Unable to load VTK module.  Cannot parse raw cdb data!')
            return
            
        if self.CheckRaw():
            raise Exception('Missing key data.  Cannot parse into unstructured grid.')            
            
        # Convert to vtk style arrays
        if use_cython:
            self.data = CDBparser.ParseForFEM(self.raw)
        else:
            self.data = PythonParser.ParseForFEM(self.raw)
            

        # Create unstructured grid
        self.uGrid = Utilities.MakeuGrid(self.data['offset'], self.data['cells'], 
                                         self.data['cell_type'],
                                         self.data['nodes'][:, :3])

        # Store original ANSYS cell and node numbering
        Utilities.AddPointScalars(self.uGrid, self.data['orignode'], 'ANSYSnodenum')

        # Extract ANSYS element numbering and store
        ansyselem = self.raw['enum'].compress(self.data['elemused'])
        Utilities.AddCellScalars(self.uGrid, ansyselem, 'ANSYSelemnum')

        # Add node components to unstructured grid
        ibool = np.empty(self.uGrid.GetNumberOfPoints(), dtype=np.int8)
        for comp in self.data['node_comps']:
            ibool[:] = 0 # reset component array
            ibool[self.data['node_comps'][comp]] = 1          
            Utilities.AddPointScalars(self.uGrid, ibool, comp)
            
        return self.data, self.uGrid, self.data['cellarr'], self.data['ncellpts']
        
        
    def AddThickness(self):
        """
        Adds 'thickness' point scalars to uGrid
        
        Assumes that thickness is stored as SURF154 elements
        
        """
        # Extract and add to uGrid.  Compress so that unused nodes are ignored
        t = ExtractThickness(self.raw)[self.data['nodenum']]
        Utilities.AddPointScalars(self.uGrid, t, 'thickness', False)
        
        
    def Plot(self):
        """ Plot unstructured grid """
        if vtk_failed:
            raise Exception('VTK not loaded!')
        
        if hasattr(self, 'uGrid'):
            Plotting.Plot(self.uGrid)
        else:
            raise Exception('Unstructred grid not generated.  Run ParseVTK first.')

    def CheckRaw(self):
        """ Check if raw data can be converted into an unstructured grid """
        try:
            self.raw['elem'][0, 0]
            self.raw['enum'][0]
        except:
            warnings.warn('Missing EBLOCK!')
            return 1
            
def ExtractThickness(raw):
    """
    Extract thickness from raw element data:

    Assumes that thickness is stored as a real constant for SURF154 elements

    The "thickness" of nodes belonging to multiple SURF154 elements will is 
    averaged   
    
    """
    
    ekey = raw['ekey']
    etype = raw['etype']
    nnum = raw['nnum']
    e_rcon = raw['e_rcon'] # real constants
    rdat = raw['rdat']
    rnum = raw['rnum']
    
    # Assemble node thickness array (nodes belonging to SURF154)
    ety = np.in1d(ekey[:, 1], 154)
    maskC = np.in1d(etype, ekey[ety, 0])
    
    # Create thickness array
    maxnode = nnum.max() + 1
    t = np.zeros(maxnode)
    a = np.zeros(maxnode, np.int32) # number of entries in thickness array
    
    if np.any(maskC):

        # Reduce element matrix to maskC elements and first four nodes
        elem = raw['elem'][maskC, :4] 
        e_rcon = e_rcon[maskC]

        # Work by getting the nodes matching each real constant
        for i in range(len(rnum)):
            # Get all surf154 elements matching the real constant
            idx = elem[e_rcon == rnum[i]] # get node indices
            
            # Add thickness            
            t[idx] += rdat[i][6]
            a[idx] += 1
        
        # normalize thickness by number of entires
        a[a == 0] = 1 # avoid divide by zero
        t /= a
    
    # keep linear nodes
    t = t[nnum]

    return t
