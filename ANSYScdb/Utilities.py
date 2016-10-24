"""
VTK based supporting functions

"""
import vtk
from vtk.util import numpy_support as VN
import numpy as np

# Determine if using vtk > 5
new_vtk = vtk.vtkVersion().GetVTKMajorVersion() > 5

def ExtractExteriorTri(uGrid):
    """
    Return surface mesh from an unstructured grid
    
    INPUTS
    uGrid (vtkUnstructuredGrid)
        
    OUTPUTS
    surf (vtkPolyData)
        vtkPolyData surface containing array 'vtkOriginalPointIds' relating the
        points of extsurf and uGrid
    
    """    

    # Extract surface mesh
    surf = vtk.vtkDataSetSurfaceFilter()
    SetVTKInput(surf, uGrid)        
    surf.PassThroughPointIdsOn()
    surf.PassThroughCellIdsOn()
    surf.Update()
    return surf.GetOutput()
    

def MakeuGrid(offset, cells, cell_type, nodes):
    """ Create VTK unstructured grid """
    
    # Check inputs (necessary since offset and cells can int32)    
    if offset.dtype != 'int64':
        offset = offset.astype(np.int64)
    
    if cells.dtype != 'int64':
        cells = cells.astype(np.int64)
    
    # Get number of cells
    ncells = len(cell_type)
    
    # Convert to vtk arrays
    cell_type = VN.numpy_to_vtk(cell_type, deep=True)
    offset = VN.numpy_to_vtkIdTypeArray(offset, deep=True)
    
    vtkcells = vtk.vtkCellArray()
    vtkcells.SetCells(ncells, VN.numpy_to_vtkIdTypeArray(cells, deep=True))
    
    # Convert points to vtkfloat object
    vtkArray = VN.numpy_to_vtk(np.ascontiguousarray(nodes), deep=True)
    points = vtk.vtkPoints()
    points.SetData(vtkArray)
    
    # Create unstructured grid
    uGrid = vtk.vtkUnstructuredGrid()
    uGrid.SetPoints(points)
    uGrid.SetCells(cell_type, offset, vtkcells)
    
    return uGrid


def AddCellScalars(grid, scalars, name, setactive=False):
    """ Adds cell scalars to uGrid """
    vtkarr = VN.numpy_to_vtk(np.ascontiguousarray(scalars), deep=True)
    vtkarr.SetName(name)
    grid.GetCellData().AddArray(vtkarr)
    if setactive:
        grid.GetCellData().SetActiveScalars(name)
    
    
def AddPointScalars(mesh, scalars, name, setactive=False):
    """
    Adds point scalars to a VTK object or structured/unstructured grid """
    vtkarr = VN.numpy_to_vtk(np.ascontiguousarray(scalars), deep=True)
    vtkarr.SetName(name)
    mesh.GetPointData().AddArray(vtkarr)
    if setactive:
        mesh.GetPointData().SetActiveScalars(name)

def ReturnCells(grid):
    """ Returns a numpy array of cells from a vtk grid object """
    return VN.vtk_to_numpy(grid.GetCells().GetData())


def GetPointScalars(vobj, name):
    """ Returns point scalars of a vtk object """
    return VN.vtk_to_numpy(vobj.GetPointData().GetArray(name))


def CopyGrid(grid):
    """ Copies a vtk structured, unstructured, or polydata object """
    if isinstance(grid, vtk.vtkUnstructuredGrid):
        gridcopy = vtk.vtkUnstructuredGrid()
    elif isinstance(grid, vtk.vtkStructuredGrid):
        gridcopy = vtk.vtkStructuredGrid()
    elif isinstance(grid, vtk.vtkPolyData):
        gridcopy = vtk.vtkPolyData()
        
    gridcopy.DeepCopy(grid)
    return gridcopy
    

def PerformLandmarkTrans(sland, tland):
    """ Performs a landmark transformation between two sets of points """
    
    slandvtk = MakevtkPoints(sland)
    tlandvtk = MakevtkPoints(tland)

    landtrans = vtk.vtkLandmarkTransform()
    landtrans.SetSourceLandmarks(slandvtk)
    landtrans.SetTargetLandmarks(tlandvtk)
    landtrans.SetModeToRigidBody()
    landtrans.Update()

    return landtrans


def MakevtkPoints(numpypts):
    """ Convert numpy points to vtkPoints """
    vtkArray = VN.numpy_to_vtk(np.ascontiguousarray(numpypts), deep=True)
    vpts = vtk.vtkPoints()
    vpts.SetData(vtkArray)
    return vpts
    

def SetVTKInput(obj, inp):
    """ Accounts for version discrepancy between VTK versions in input method """
    if new_vtk:
        obj.SetInputData(inp)
    else:
        obj.SetInput(inp)
        
        
def GetMeshAreaVol(mesh):
    """ Returns volume and area of a triangular mesh """
    mprop = vtk.vtkMassProperties()
    SetVTKInput(mprop, mesh) 
    mprop.Update() 
    return mprop.GetSurfaceArea(), mprop.GetVolume()


def GetPoints(mesh, datatype=[]):
    """ returns points from a mesh as numpy array """
    points = VN.vtk_to_numpy(mesh.GetPoints().GetData())
    
    if datatype:
        if points.dtype != datatype:
            return points.astype(datatype)
    
    return points


def GetFaces(mesh):
    """ returns points from a polydata polydata object and return as a numpy int array """
    return VN.vtk_to_numpy(mesh.GetPolys().GetData()).reshape(-1, 4)[:, 1:]


def TriFilter(mesh):
    """ Returns an all triangle mesh of a polydata """
    trifilter = vtk.vtkTriangleFilter()
    SetVTKInput(trifilter, mesh)
    trifilter.PassVertsOff()
    trifilter.PassLinesOff()
    trifilter.Update()
    
    # Return triangular exterior surface mesh and nontriangular
    return trifilter.GetOutput()    
        
    
def SetVTKCoords(VTKobject, pts):
    """ Sets points on a VTK object """
    print('SetVTKCoords degenerated.  Use SetPoints')
    vtkfloat = VN.numpy_to_vtk(np.ascontiguousarray(pts), deep=True)
    vtkpointArray = vtk.vtkPoints()
    vtkpointArray.SetData(vtkfloat)
    VTKobject.SetPoints(vtkpointArray)
    
    
def SetPoints(VTKobject, points):
    """ Sets points on a VTK object.  Same as SetVTKCoords, but better named """
    vtkfloat = VN.numpy_to_vtk(np.ascontiguousarray(points), deep=True)
    vtkpointArray = vtk.vtkPoints()
    vtkpointArray.SetData(vtkfloat)
    VTKobject.SetPoints(vtkpointArray)
        
    
def GetPointNormals(mesh):
    """ Compute point normals froma a vtkPolyData object """
    raise Exception('Depreciated: Use Meshutil GetPointNormals')
    vtknormals = vtk.vtkPolyDataNormals()
    SetVTKInput(vtknormals, mesh)
    vtknormals.SplittingOff()
    vtknormals.ComputeCellNormalsOff()
    vtknormals.ComputePointNormalsOn()
    try:
        vtknormals.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)
    except:
        pass
    vtknormals.Update()
    return VN.vtk_to_numpy(vtknormals.GetOutput().GetPointData().GetNormals())
    
    
def GetCurvature(mesh, curvtype='Mean'):
    """
    Returns the pointwise curvature of a mesh
    
    Availble options include:
        Mean
        Gaussian
        Maximum
    """
    
    # Create curve filter and compute curvature
    curvefilter = vtk.vtkCurvatures()
    SetVTKInput(curvefilter, mesh)
    if curvtype == 'Mean':
        curvefilter.SetCurvatureTypeToMean()
    elif curvtype == 'Gaussian': 
        curvefilter.SetCurvatureTypeToGaussian()
    elif curvtype == 'Maximum':
        curvefilter.SetCurvatureTypeToMaximum()
    else:
        curvefilter.SetCurvatureTypeToMinimum()
    curvefilter.Update()
    
    # Store curvature
    curves = curvefilter.GetOutput()
    curvature = VN.vtk_to_numpy(curves.GetPointData().GetScalars())
#    del curvefilter
    
    return curvature, curves
    

def VertFacetoPoly(v, f, removeduplicates=True):
    """ Creates a vtk polydata object given points and triangular faces """
    
    # Check inputs
    if f.dtype != 'int64':
        f = f.astype(np.int64)
        
    if not v.flags['C_CONTIGUOUS']:
        v = np.ascontiguousarray(v)
        
    # Remove duplicate verticies and renumber v and f
    if removeduplicates:
        b = v.view(np.dtype((np.void, v.dtype.itemsize*v.shape[1])))
        _, idx, idx2 = np.unique(b, return_index=True, return_inverse=True)
        if idx.size != v.shape[0]:
            v = v.take(idx, 0)
            f = idx2[f]
        else:
            idx = np.arange(v.shape[0])
    
    # Convert points to vtkfloat object
    vtkArray = VN.numpy_to_vtk(v, deep=True)
    points = vtk.vtkPoints()
    points.SetData(vtkArray) 
    
    # Convert faces to vtk cells object
    cells = np.empty((f.shape[0], 4), dtype=np.int64)
    cells[:, -3:] = f
    cells[:, 0] = 3
    cells = cells.ravel()
    if not cells.flags['C_CONTIGUOUS']:
        cells = np.ascontiguousarray(cells)
    vtkcells = vtk.vtkCellArray()
    vtkcells.SetCells(cells.shape[0], VN.numpy_to_vtkIdTypeArray(cells, deep=True))

    # Create polydata object
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(points)
    mesh.SetPolys(vtkcells)
    
     # Return original indices if requested
    if removeduplicates:
        return mesh, idx
    else:
        return mesh


def CreateVectorPolyData(orig, vec):
    """ Creates a vtkPolyData object composed of vectors """
        
    
    # Create vtk points and cells objects
    vpts = vtk.vtkPoints()
    vpts.SetData(VN.numpy_to_vtk(np.ascontiguousarray(orig), deep=True))
    
    npts = orig.shape[0]
    cells = np.hstack((np.ones((npts, 1), 'int'),
                       np.arange(npts).reshape((-1, 1))))
                       
    cells = np.ascontiguousarray(cells)
    vcells = vtk.vtkCellArray()
    vcells.SetCells(npts, VN.numpy_to_vtkIdTypeArray(cells, deep=True))
    
    # Create vtkPolyData object
    pdata = vtk.vtkPolyData()
    pdata.SetPoints(vpts);
    pdata.SetVerts(vcells)
    
    # Add vectors to polydata
    name='vectors'
    vtkfloat = VN.numpy_to_vtk(np.ascontiguousarray(vec), deep=True)
    vtkfloat.SetName(name)
    pdata.GetPointData().AddArray(vtkfloat)
    pdata.GetPointData().SetActiveVectors(name)
    
    # Add magnitude of vectors to polydata
    name='mag'
    scalars = (vec*vec).sum(1)**0.5
    vtkfloat = VN.numpy_to_vtk(np.ascontiguousarray(scalars), deep=True)
    vtkfloat.SetName(name)
    pdata.GetPointData().AddArray(vtkfloat)
    pdata.GetPointData().SetActiveScalars(name)
    
    return pdata    
    
    
def ApplyTransformation(mesh, trans):
    """ Apply vtk transformation to vtk mesh """

    # Convert 4x4 matrix to a transformation object if applicable
    if trans.IsA('vtkMatrix4x4'):
        transform = vtk.vtkTransform()
        transform.SetMatrix(trans)
        trans = transform
    
    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)
    SetVTKInput(transformFilter, mesh)
    transformFilter.SetTransform(trans)    
    transformFilter.Update()
    return transformFilter.GetOutput()
    
    
def ApplyTransformationSolid(grid, trans):
    """ Apply vtk transformation to vtk mesh """

    # Convert 4x4 matrix to a transformation object if applicable
    if trans.IsA('vtkMatrix4x4'):
        transform = vtk.vtkTransform()
        transform.SetMatrix(trans)
        trans = transform
    
    transformFilter = vtk.vtkTransformFilter()
    SetVTKInput(transformFilter, grid)
    transformFilter.SetTransform(trans)    
    transformFilter.Update()
    return transformFilter.GetOutput()
    
    
def AlignMesh(fixed, moving, nlan=[], niter=[], mmdist=[]):
    """ Aligns the moving mesh to the fixed mesh """
    
    # Perform ICP
    vtkICP = vtk.vtkIterativeClosestPointTransform()
    vtkICP.SetSource(moving)
    vtkICP.SetTarget(fixed)
    vtkICP.GetLandmarkTransform().SetModeToRigidBody()
    
    # Set settings if applicable
    if nlan:
        vtkICP.SetMaximumNumberOfLandmarks(nlan)
    if niter:
        vtkICP.SetMaximumNumberOfIterations(niter)
    if mmdist:
        vtkICP.SetMaximumMeanDistance(mmdist)
    vtkICP.StartByMatchingCentroidsOn()
    vtkICP.Modified()
    vtkICP.Update()
    
    return vtkICP

    # Update target
#    ApplyTransformation(moving, vtkICP)
    
    
###############################################################################
# File I/O
###############################################################################
    
def WriteMesh(filename, mesh, ftype=[], binary=True):
    """ Writes a VTK mesh to one of the following formats:
        ply, stl, vtk
    """
    
    if not ftype:
        ftype = filename[-3:]
    
    # Get filetype
    if ftype == 'ply':
        writer = vtk.vtkPLYWriter()
    elif ftype == 'stl':
        writer = vtk.vtkSTLWriter()
    elif ftype == 'vtk':
        writer = vtk.vtkPolyDataWriter()
    else:
        raise Exception('Unknown file type')
        
    # Write
    writer.SetFileName(filename)
    SetVTKInput(writer, mesh)
    if binary:
        writer.SetFileTypeToBinary()
    else:
        writer.SetFileTypeToASCII()
    writer.Write()
    
    
def LoadMesh(filename):
    """ Reads mesh from file """

    # Get extension
    fext = filename[-3:].lower()

    # Select reader
    if fext == 'ply':
        reader = vtk.vtkPLYReader()
        
    elif fext == 'stl':
        reader = vtk.vtkSTLReader()
        
    else:
        raise Exception('Unknown file type')
    
    # Load file
    reader.SetFileName(filename) 
    reader.Update()
    return reader.GetOutput()
    