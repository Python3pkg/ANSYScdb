import vtk
import numpy as np
import colorsys

import Utilities
from Utilities import new_vtk
from vtk.util import numpy_support as VN

    
class PlotClass(object):
    """ Simple interface to VTK's underlying ploting """
    
    def __init__(self):

        # Add FEM Actor to renderer window
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(0.3, 0.3, 0.3)
        
        self.renWin = vtk.vtkRenderWindow()
        self.renWin.AddRenderer(self.ren)
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.renWin)
        
        # Allow user to interact
        istyle = vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(istyle)


    def AddMesh(self, meshin, color=[1, 1, 1], style='', scalars=[], name='',
                rng=[], stitle='', showedges=True, psize=5, opacity=1,
                linethick=[]):
        """ Adds an actor to the renderwindow """
                
        # Create mapper
        mapper = vtk.vtkDataSetMapper()
                
        # Add scalars if they exist
        isscalars = False
        nscalars = len(scalars)
        if nscalars == meshin.GetNumberOfPoints():
            mesh = Utilities.CopyGrid(meshin)
            Utilities.AddPointScalars(mesh, scalars, name)
            isscalars = True
            mapper.SetScalarModeToUsePointData()
            

        elif nscalars == meshin.GetNumberOfCells():
            mesh = Utilities.CopyGrid(meshin)
            Utilities.AddCellScalars(mesh, scalars, name)
            isscalars = True
            mapper.SetScalarModeToUseCellData()
            
        else:
            mesh = meshin
                    
        # Set scalar range
        if isscalars:
            if not rng:
                rng = [np.min(scalars), np.max(scalars)]
            mapper.SetScalarRange(rng[0], rng[1])
        
        # Set Scalar
        Utilities.SetVTKInput(mapper, mesh)
        
        # Create Actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        if style == 'wireframe':
            actor.GetProperty().SetRepresentationToWireframe()
        elif style == 'points':
            actor.GetProperty().SetRepresentationToPoints()
            actor.GetProperty().SetPointSize(psize)
        else:
            actor.GetProperty().SetRepresentationToSurface()
            
        if showedges:
            actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetOpacity(opacity)
        actor.GetProperty().LightingOff()
        
        if style == 'wireframe' and linethick:
            actor.GetProperty().SetLineWidth(linethick) 

        
        # Add to renderer
        self.ren.AddActor(actor)
        
        # Add scalar bar
        if stitle:
            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(mapper.GetLookupTable())
            scalarBar.SetTitle(stitle)
            scalarBar.SetNumberOfLabels(5)    
            self.ren.AddActor(scalarBar)


    def AddLines(self, lines, color=[1, 1, 1], width=5):
        """ Adds an actor to the renderwindow """
                
        # Create mapper and add lines
        mapper = vtk.vtkDataSetMapper()
        Utilities.SetVTKInput(mapper, lines)
        
        # Create Actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetLineWidth(width); 
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetColor(color)
        actor.GetProperty().LightingOff()
        
        # Add to renderer
        self.ren.AddActor(actor)
        

    def AddPoints(self, points, color=[1, 1, 1], psize=5):
        
        # Convert to points actor if points is a numpy array
        if type(points) == np.ndarray:
            npoints = points.shape[0]
            
            # Make VTK cells array
            cells = np.hstack((np.ones((npoints, 1)), 
                               np.arange(npoints).reshape(-1, 1)))
            cells = np.ascontiguousarray(cells, dtype=np.int64)
            vtkcells = vtk.vtkCellArray()
            vtkcells.SetCells(npoints, VN.numpy_to_vtkIdTypeArray(cells, deep=True))
            
            # Convert points to vtk object
            vtkPoints = Utilities.MakevtkPoints(points)
            
            # Create polydata
            pdata = vtk.vtkPolyData()
            pdata.SetPoints(vtkPoints)
            pdata.SetVerts(vtkcells)
            
        # Create mapper and add lines
        mapper = vtk.vtkDataSetMapper()
        Utilities.SetVTKInput(mapper, pdata)
        
        # Create Actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetPointSize(psize); 
        actor.GetProperty().SetColor(color)
        actor.GetProperty().LightingOff()
        
        self.ren.AddActor(actor)
                
        
    def GetCameraPosition(self):
        """ Returns camera position of active render window """
        camera = self.ren.GetActiveCamera()
        pos = camera.GetPosition()
        fpt = camera.GetFocalPoint()
        vup = camera.GetViewUp()
        return [pos, fpt, vup]
        

    def SetCameraPosition(self, cameraloc):
        """ Set camera position of active render window """
        camera = self.ren.GetActiveCamera()
        camera.SetPosition(cameraloc[0])
        camera.SetFocalPoint(cameraloc[1]) 
        camera.SetViewUp(cameraloc[2])        
        

    def SetBackground(self, bcolor):
        """ Sets background color """
        self.ren.SetBackground(bcolor)
        
        
    def AddLegend(self, entries, bcolor=[0.5, 0.5, 0.5], border=False):
        """
        Adds a legend to render window.  Entries must be a list containing
        one string and color entry for each item
        """
        
        legend = vtk.vtkLegendBoxActor()
        legend.SetNumberOfEntries(len(entries))
        
        c = 0
        nulldata = vtk.vtkPolyData()
        for entry in entries:
            legend.SetEntry(c, nulldata, entry[0], entry[1])
            c += 1
        
        legend.UseBackgroundOn()
        legend.SetBackgroundColor(bcolor)
        if border:
            legend.BorderOn()
        else:
            legend.BorderOff()
        
        # Add to renderer
        self.ren.AddActor(legend)
        
        
    def Plot(self, title=''):
        """ Renders """
        if title:
            self.renWin.SetWindowName(title)
            
        # Render
        self.iren.Initialize()
        self.renWin.Render()
        self.iren.Start()
        
        
    def AddActor(self, actor):
        """ Adds actor to render window """
        self.ren.AddActor(actor)
        
        
    def AddAxes(self):
        """ Add axes widget """
        axes = vtk.vtkAxesActor()
        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOrientationMarker(axes)
        widget.SetInteractor(self.iren)
        widget.SetViewport(0.0, 0.0, 0.4, 0.4)
        widget.SetEnabled(1)
        widget.InteractiveOn()
        


def CreateArrowsActor(pdata):
    """ Creates an actor composed of arrows """
    
    # Create arrow object
    arrow = vtk.vtkArrowSource()
    arrow.Update()
    glyph3D = vtk.vtkGlyph3D()
    if new_vtk:
        glyph3D.SetSourceData(arrow.GetOutput())
        glyph3D.SetInputData(pdata)
    else:
        glyph3D.SetSource(arrow.GetOutput())
        glyph3D.SetInput(pdata)
    glyph3D.SetVectorModeToUseVector()
    glyph3D.Update()
    
    # Create mapper    
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection(glyph3D.GetOutputPort())
    
    # Create actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().LightingOff()

    return actor
    
        
def PlotCurvature(mesh, curvtype):
    """
    Plots curvature
    Availble options for curvtype:
        'Mean'
        'Gaussian'
        'Maximum  '  
    
    """
    
    # Get curvature values and plot
    c = Utilities.GetCurvature(mesh, curvtype)[0]
    pobj = PlotClass()
    pobj.AddMesh(mesh, scalars=c)
    pobj.Plot(); del pobj

    
def PlotGrids(grids, wFEM=False):
    """
    Creates a plot of several grids as wireframes.  When wFEM is true, the first
    grid is a white solid
    """
    
    # Make grid colors
    N = len(grids)
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    colors = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    
    pobj = PlotClass()
    for i in range(len(grids)):
        if not i and wFEM: # Special plotting for first grid
            pobj.AddMesh(grids[i])
        else:
            pobj.AddMesh(grids[i], color=colors[i], style='wireframe')
    
    # Render plot and delete when finished
    pobj.SetBackground([0.8, 0.8, 0.8])
    pobj.Plot(); del pobj


def PlotPoly(mesh, representation='surface', color=[1, 1, 1]):
    """ Plots vtk unstructured grid or poly object """
    pobj = PlotClass()
    pobj.AddMesh(mesh, color, style=representation)
    pobj.Plot()
    del pobj


def Plot(mesh, representation='surface', color=[1, 1, 1]):
    """ calls PlotPoly """
    PlotPoly(mesh, representation, [1, 1, 1])
    
    
def PlotEdges(mesh, angle, width=10):
    """ Plots edges of a mesh """
    
    # Extract edge points from a mesh
    edges = Utilities.GetEdgePoints(mesh, angle, False)
        
    # Render
    pobj = PlotClass()
    pobj.AddLines(edges, [0, 1, 1], width)
    pobj.AddMesh(mesh)
    pobj.Plot(); del pobj
    
    
def PlotBoundaries(mesh):
    """ Plots boundaries of a mesh """
    featureEdges = vtk.vtkFeatureEdges()
    Utilities.SetVTKInput(featureEdges, mesh)
    
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOff()
    featureEdges.ManifoldEdgesOff()
    
    edgeMapper = vtk.vtkPolyDataMapper();
    edgeMapper.SetInputConnection(featureEdges.GetOutputPort());
    
    edgeActor = vtk.vtkActor();
    edgeActor.GetProperty().SetLineWidth(5);
    edgeActor.SetMapper(edgeMapper)

    mapper = vtk.vtkDataSetMapper()
    Utilities.SetVTKInput(mapper, mesh)

    # Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().LightingOff()    
        
    # Render
    pobj = PlotClass()
    pobj.AddActor(actor)
    pobj.Plot(); del pobj
    
    