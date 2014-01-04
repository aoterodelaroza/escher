# ESCHERpy -- A computational chemistry workflow tool
# Copyright (C) 2013 Daniel Menendez
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from math import sqrt, acos, degrees
import vtk
import numpy as np
from logging import getLogger
from templates.blender import blenderinput
log = getLogger('escherlog')

vtkversion = vtk.vtkVersion().GetVTKSourceVersion().split()[-1]
vtkversion = int(vtkversion.split(".")[0])
vtk6 = (vtkversion == 6)



class Representation(object):

    '''
    function rep = representation()

    representation - create an empty rep structure.

    Output:
    {rep}: the empty representation.
    '''

    def __init__(self,repi=None,camangle=np.matrix([80,75,45]),zoom=3):

        self.name = ""

        self.nball = 0
        self.nstick = 0

        self.ntriangle = 0
        self.nvertex = 0
        self.triangle = []
        self.vertex = []

        self.cam = type('Cam', (), {})
        self.nlight = 0
        self.light = []
        self.bgcolor = np.zeros([1,3])

    
    @staticmethod
    def ball(radius, position, color): #, name):
        '''
        Make a sphere of a certain radius, position
        and color. This will represent an atom.
        '''
        # create the object to be drawn
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius( radius )
        sphere.SetThetaResolution(40)
        sphere.SetPhiResolution(40)
        sphere.SetCenter(position)
        # Convert the sphere into polygons
        sphereMapper = vtk.vtkPolyDataMapper()
        sphereMapper.SetInputConnection(sphere.GetOutputPort())

        #Create an actor for the sphere 
        sphereActor = vtk.vtkActor()
        sphereActor.SetMapper(sphereMapper)
        (sphereActor.GetProperty()).SetColor(color)

        ## Add atomic name label
        #atext = vtk.vtkVectorText()
        #atext.SetText(name)
        #textMapper = vtk.vtkPolyDataMapper()
        #textMapper.SetInputConnection(atext.GetOutputPort())
        #textActor = vtk.vtkFollower()
        #textActor.SetMapper(textMapper)
        #textActor.SetScale(0.8, 0.8, 0.8)
        ##postxt = position
        ##postxt[1] = position[1] - radius
        #textActor.AddPosition(position)
        #textActor.GetProperty().SetColor(0,0,0) # (R,G,B)
        #return sphereActor, textActor

        return sphereActor

    @staticmethod
    def stick(v1, v2):
        '''
        It represents a bond within a cylinder
        '''

        # vector that represents the bond
        v = [0.,0.,0.]
        center = [0.,0.,0.]
        for i in range(3):
            v[i] = v2[i] - v1[i]
            center[i] = v1[i] + v[i]/2.

        height = sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        if (v[0]**2 + v[1]**2 + v[2]**2) > 1e-8:
            phi = acos(v[1]/sqrt(v[0]**2 + v[1]**2 + v[2]**2))
        else:
            phi = 0.

        # create cylinder
        source = vtk.vtkCylinderSource()
        source.SetRadius(0.08)
        source.SetHeight(height)
        source.SetResolution(60)
        source.SetCenter(0,0,0)
         
        # create a transform that rotates the bond cylinder
        transform = vtk.vtkTransform()
        # translate axes to the center of bond
        transform.Translate(center)
        # Rotatees phi around axis perpendicular to bond and axis y
        transform.RotateWXYZ(degrees(phi),[v[2],0,-v[0]])
        transformFilter=vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(source.GetOutputPort())
        transformFilter.Update()
         
        # mapper for original cylinder
        cylinMapper1 = vtk.vtkPolyDataMapper()
        cylinMapper1.SetInputConnection(source.GetOutputPort())
         
        # another mapper for the rotated cylinder
        cylinMapper2 = vtk.vtkPolyDataMapper()
        cylinMapper2.SetInputConnection(transformFilter.GetOutputPort())
         
        # actor for original cylinder
        #actor1 = vtk.vtkActor()
        #actor1.SetMapper(cylinMapper1)
         
        # actor for rotated cylinder
        actor2 = vtk.vtkActor()
        actor2.SetMapper(cylinMapper2)
         
        #actor1.GetProperty().SetColor(1,0,0) # (R,G,B)
        actor2.GetProperty().SetColor(1,1,1) # (R,G,B)

        #return actor1, actor2
        return actor2

    def window(self):
        '''
        Creates a VTK rendering window
        '''
        log.debug('Setting VTK rendering window')

        # create a rendering window and renderer
        ren = vtk.vtkRenderer()
        # camera (optional)
        global camera
        camera = vtk.vtkCamera()
        camera.SetPosition(-30,0,0)
        camera.SetFocalPoint(0,0,0)
        camera.ParallelProjectionOn()
        ren.SetActiveCamera(camera)
        log.info('Initial camera orientation: {}'.format(camera.GetOrientation()))
        log.info('Initial camera orientation: {}'.format(camera.GetOrientationWXYZ()))
        log.info('Initial camera-focus distance: {}'.format(camera.GetDistance()))
        log.info('Initial camera position: {}'.format(camera.GetPosition()))
        log.info('Initial camera focal point: {}'.format(camera.GetFocalPoint()))
#        ren.ResetCamera()
#        ren.GetActiveCamera().Azimuth(30)
#        ren.GetActiveCamera().Elevation(20)
#        ren.GetActiveCamera().Dolly(2.8)
#        ren.ResetCameraClippingRange()
        ren.SetBackground(1,1,1)

        global renWin
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        renWin.SetSize(1700,1200)
        renWin.StereoCapableWindowOn()
        #self.renWin = renWin

        iren = vtk.vtkRenderWindowInteractor()
        iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        iren.RemoveObservers('MiddleButtonPressEvent')
        iren.AddObserver("MiddleButtonPressEvent",self.exporter)
        iren.SetRenderWindow(renWin)
        self.ren = ren
        self.iren = iren

    def start(self):
        '''
        Starts rendering the VTK window
        '''
        log.debug('VTK render window starting')
        # enable user interface interactor
        self.iren.Initialize()
        self.iren.Start()


        ## image
        #wif = vtk.vtkWindowToImageFilter()
        #wif.SetInputData(renWin)
        #png = vtk.vtkPNGWriter()
        #png.SetInputConnection(wif.GetOutputPort())
        #png.SetFileName('img.png')
        #png.Write()
    
    def imagedata(self, orig, delta, dims):
        '''
        Creates a VTKImageData object from
        given voxel specifications

        :param list orig: corner origin of the voxel
                          [x0, y0, z0]
        :param list delta: spacing between points of the voxel
                          [dx, dy, dz]
        :param list dims: number of points in each direction of the voxel
                          [nx, ny, nz]


        '''

        log.debug('Scalar grid volume data')
        # construct the volume
        grid = vtk.vtkImageData()
        grid.SetOrigin(orig) 
        grid.SetSpacing(delta)
        grid.SetDimensions(dims) # number of points in each direction

        self.grid = grid

    def array(self, dims, scalars, name):

        log.debug('vtkDoubleArray creation')
        array = vtk.vtkDoubleArray()
        array.SetNumberOfComponents(1) # this is 3 for a vector
        array.SetNumberOfTuples(self.grid.GetNumberOfPoints())

        nd = 0
        for k in range(dims[2]):
            for j in range(dims[1]):
                for i in range(dims[0]):
                    array.SetValue(nd, scalars[i,j,k])
                    nd = nd + 1

        array.SetName(name)

        return array

    def outline(self):
        '''
        Outline box that indicates the
        volumetric limits
        '''

        bounds = vtk.vtkOutlineFilter()
        if vtk6:
            bounds.SetInputData(self.grid)
        else:
            bounds.SetInput(self.grid)

        boundsMapper = vtk.vtkPolyDataMapper()
        boundsMapper.SetInputConnection(bounds.GetOutputPort())

        boundsActor = vtk.vtkActor()
        boundsActor.SetMapper(boundsMapper)
        boundsActor.GetProperty().SetColor(0, 0, 0)
        self.ren.AddActor(boundsActor)

    def lookuptable(self, rangecolor, colortable):
        '''
        A color map creator
        '''

        #lutNCI = vtk.vtkLookupTable()
        #lutNCI.SetNumberOfColors(3)
        #lutNCI.SetHueRange(0.667, 0.0)
        #lutNCI.SetSaturationRange(1.0, 1.0)
        #lutNCI.SetValueRange(0.8, 0.8)
        #lutNCI.Build()

        # color palette
        # -3.0 -> blue
        #  0.0 -> green
        #  3.0 -> red
        colorNCI = vtk.vtkColorTransferFunction()
        colorNCI.AddRGBPoint(rangecolor[0],*colortable[0])
        colorNCI.AddRGBPoint(0.,*colortable[1])
        colorNCI.AddRGBPoint(rangecolor[1],*colortable[2])
        #colorNCI.AddRGBPoint(0.0,0.0,1.0,0.0)
        #colorNCI.AddRGBPoint(rangecolor[1],1.0,0.0,0.0)
        self.colorNCI = colorNCI

    def scalarbar(self):
        '''
        Renders a scalar bar showing the colorbar used
        with :py:meth:`lookuptable`
        '''

        scalarBar = vtk.vtkScalarBarActor()
        scalarBar.SetLookupTable(self.mapper.GetLookupTable())
        scalarBar.SetTitle("sign(lambda_2)rho")
        scalarBar.SetOrientationToHorizontal()
        scalarBar.GetLabelTextProperty().SetColor(0,0,1)
        scalarBar.GetTitleTextProperty().SetColor(0,0,1)
        # position it in window
        coord = scalarBar.GetPositionCoordinate()
        coord.SetCoordinateSystemToNormalizedViewport()
        coord.SetValue(0.1,0.05)
        scalarBar.SetWidth(.8)
        scalarBar.SetHeight(.085)

        self.ren.AddActor(scalarBar)


    def marchingcubes(self, isovalue):
        '''
        Extracts an isosurface from volumetric data
        '''

        log.debug('Marching Cubes algorithm contour filter')
        # ContourFilter or MarchingCubes
        #isoExtractor = vtk.vtkContourFilter()
        isoExtractor = vtk.vtkMarchingCubes()
        if vtk6:
            isoExtractor.SetInputData(self.grid)
        else:
            isoExtractor.SetInput(self.grid)
        # isosurface id and value
        isoExtractor.SetValue(0, isovalue)
        #isoExtractor.GenerateValues(4, isovalue-0.05, isovalue+0.05)
        isoExtractor.ComputeGradientsOn()

        isoNormals = vtk.vtkPolyDataNormals()
        isoNormals.SetInputConnection(isoExtractor.GetOutputPort())
        isoNormals.SetFeatureAngle(35.0)
        isoNormals.ConsistencyOn()
        isoNormals.SplittingOn()

        isoStripper = vtk.vtkStripper()
        isoStripper.SetInputConnection(isoNormals.GetOutputPort())

        isoLocator = vtk.vtkCellLocator()
        isoLocator.SetDataSet(isoExtractor.GetOutput())
        isoLocator.LazyEvaluationOn()

        isoMapper = vtk.vtkPolyDataMapper()
        isoMapper.SetInputConnection(isoStripper.GetOutputPort())
        isoMapper.ScalarVisibilityOn()
        self.mapper = isoMapper
        #isoMapper.InterpolateScalarsBeforeMappingOn()
        #isoMapper.SetScalarModeToUsePointData()
        #isoMapper.ColorByArrayComponent('dens', 0) # doesn't work
        #isoMapper.GetArrayName('dens') # doesn't work

        isoProperty = vtk.vtkProperty()
        isoProperty.SetColor(0.2,1,0.2) # ~green
        #isoProperty.SetInterpolationToPhong()

        iso = vtk.vtkActor()
        iso.SetMapper(isoMapper)
        iso.SetProperty(isoProperty)
        #iso.GetProperty().SetRepresentationToWireframe()
        #iso.GetProperty().SetOpacity(.6)

        #self.ren.AddViewProp(iso)
        self.ren.AddActor(iso)

    def surface(self,nv, vxyz, color):
        """
        Construct a surface from a set of vertices
        """

        points = vtk.vtkPoints()
        for i in range(0, nv):
            points.InsertPoint(i, vxyz[i])
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)

        delny = vtk.vtkDelaunay3D()
        if vtk6:
            delny.SetInputData(polydata)
        else:
            delny.SetInput(polydata)
        delny.SetTolerance(0.01)
        delny.SetAlpha(0.0)
        delny.BoundingTriangulationOff()

        surfmap = vtk.vtkDataSetMapper()
        surfmap.SetInputConnection(delny.GetOutputPort())

        triangulation = vtk.vtkActor()
        triangulation.SetMapper(surfmap)
        triangulation.GetProperty().SetColor(color)

        self.ren.AddActor(triangulation)



    def _gridvolume(self, orig, delta, dims, scalars):

        log.debug('Scalar grid volume data')
        # construct the volume
        grid = vtk.vtkImageData()
        grid.SetOrigin(orig) 
        grid.SetSpacing(delta)
        grid.SetDimensions(dims) # number of points in each direction
        array = vtk.vtkDoubleArray()
        array.SetNumberOfComponents(1) # this is 3 for a vector
        array.SetNumberOfTuples(grid.GetNumberOfPoints())

        nd = 0
        for k in range(dims[2]):
            for j in range(dims[1]):
                for i in range(dims[0]):
                    array.SetValue(nd, scalars[i,j,k])
                    nd = nd + 1

        array.SetName("scalars")
        grid.GetPointData().SetScalars(array)
        #grid.GetPointData().AddArray(array)
        self.grid = grid

        bounds = vtk.vtkOutlineFilter()
        bounds.SetInputData(grid)

        boundsMapper = vtk.vtkPolyDataMapper()
        boundsMapper.SetInputConnection(bounds.GetOutputPort())

        boundsActor = vtk.vtkActor()
        boundsActor.SetMapper(boundsMapper)
        boundsActor.GetProperty().SetColor(0, 0, 0)
        self.ren.AddActor(boundsActor)

#        array2 = vtk.vtkDoubleArray()
#        array2.SetNumberOfComponents(1) # this is 3 for a vector
#        array2.SetNumberOfTuples(self.grid.GetNumberOfPoints())
#
#        nd = 0
#        for k in range(dims2[2]):
#            for j in range(dims2[1]):
#                for i in range(dims2[0]):
#                    array2.SetValue(nd, scalars2[i,j,k])
#                    nd = nd + 1
#
#        array2.SetName("dens")
#        self.grid.GetPointData().SetScalars(array2)

    def _volumescalar(self, dims, scalars):

        array = vtk.vtkDoubleArray()
        array.SetNumberOfComponents(1) # this is 3 for a vector
        array.SetNumberOfTuples(self.grid.GetNumberOfPoints())

        nd = 0
        for k in range(dims[2]):
            for j in range(dims[1]):
                for i in range(dims[0]):
                    array.SetValue(nd, scalars[i,j,k])
                    nd = nd + 1

        array.SetName("dens")
        self.array = array
        self.grid.GetPointData().AddArray(array)

    @staticmethod
    def exporter(obj, ev):
        '''
        Exports representation state to several
        common 3D file formats while pressing the
        mouse middle button on the VTK window

        '''
        log.debug("Middle Button pressed")

        log.info('Current camera orientation: {}'.format(camera.GetOrientation()))
        log.info('Current camera orientation: {}'.format(camera.GetOrientationWXYZ()))
        log.info('Current camera-focus distance: {}'.format(camera.GetDistance()))
        log.info('Current camera position: {}'.format(camera.GetPosition()))
        log.info('Current camera focal point: {}'.format(camera.GetFocalPoint()))

        log.debug('Exporting current structure to POVRay, OBJ, VRML 2.0, OOGL.')

        # POVRay output 
        pov = vtk.vtkPOVExporter()
        pov.SetRenderWindow(renWin)
        pov.SetFileName('out.pov')
        pov.Write()

        # OBJ Wavefront output 
        # main file: <prefix>.obj 
        # textures:  <prefix>.mtl
        # It can be imported with g3dviewer
        obj = vtk.vtkOBJExporter()
        obj.SetInput(renWin)
        obj.SetFilePrefix('out')
        obj.Write()

        # VRML 2.0 output
        # It can be imported with Blender
        vrml = vtk.vtkVRMLExporter()
        vrml.SetRenderWindow(renWin)
        vrml.SetFileName('out.wrl')
        vrml.Write()

        # Geomview OOGL ouput
        oogl = vtk.vtkOOGLExporter()
        oogl.SetRenderWindow(renWin)
        oogl.SetFileName('out.oogl')
        oogl.Write()

        log.debug('Generating Blender script.')

        # Blender script to view VRML file
        # Runnable with blender -P blender.py
        blenderf = open('vrml.bpy', 'w')
        blenderf.write(blenderinput)

        return

