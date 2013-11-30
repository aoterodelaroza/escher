# Copyright (C) 2011 Victor Lua~na and Alberto Otero-de-la-Roza
#
# This octave routine is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version. See <http://www.gnu.org/licenses/>.
#
# The routine distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.

from math import sqrt, acos, degrees
import vtk
import numpy as np
from logging import getLogger
log = getLogger('escherlog')

class Representation():

    '''
    function rep = representation()

    representation - create an empty rep structure.

    Output:
    {rep}: the empty representation.
    '''

    def __init__(self,repi=None,camangle=np.matrix([80,75,45]),zoom=3,LOG=0):

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
    def ball(radius, position, color):
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
         
        # create a transform that rotates the cone
        transform = vtk.vtkTransform()
        # translate axes to the center of bond
        transform.Translate(center)
        # Rotatees phi around axis perpendicular to bond and axis y
        transform.RotateWXYZ(degrees(phi),[v[2],0,-v[0]])
        transformFilter=vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(source.GetOutputPort())
        transformFilter.Update()
         
        # mapper for original cone
        coneMapper1 = vtk.vtkPolyDataMapper()
        coneMapper1.SetInputConnection(source.GetOutputPort())
         
        # another mapper for the rotated cone
        coneMapper2 = vtk.vtkPolyDataMapper()
        coneMapper2.SetInputConnection(transformFilter.GetOutputPort())
         
        # actor for original cone
        actor1 = vtk.vtkActor()
        actor1.SetMapper(coneMapper1)
         
        # actor for rotated cone
        actor2 = vtk.vtkActor()
        actor2.SetMapper(coneMapper2)
         
        actor1.GetProperty().SetColor(1,0,0) # (R,G,B)
        actor2.GetProperty().SetColor(1,1,1) # (R,G,B)

        #return actor1, actor2
        return actor2

    def window(self):
        log.debug('Setting VTK rendering window')

        # create a rendering window and renderer
        ren = vtk.vtkRenderer()
        # camera (optional)
        camera = vtk.vtkCamera()
        camera.SetPosition(15,0,0)
        camera.SetFocalPoint(0,0,0)
        camera.ParallelProjectionOn()
        ren.SetActiveCamera(camera)
#        ren.ResetCamera()
#        ren.GetActiveCamera().Azimuth(30)
#        ren.GetActiveCamera().Elevation(20)
#        ren.GetActiveCamera().Dolly(2.8)
#        ren.ResetCameraClippingRange()
        ren.SetBackground(1,1,1)

        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        renWin.SetSize(1700,1200)
        renWin.StereoCapableWindowOn()
        self.renWin = renWin

        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        self.ren = ren
        self.iren = iren

    def start(self):
        log.debug('VTK render window starting')
        # enable user interface interactor
        self.iren.Initialize()
        self.iren.Start()

        pov = vtk.vtkPOVExporter()
        pov.SetRenderWindow(self.renWin)
        pov.SetFileName('out.pov')
        pov.Write()

        ## image
        #wif = vtk.vtkWindowToImageFilter()
        #wif.SetInputData(renWin)
        #png = vtk.vtkPNGWriter()
        #png.SetInputConnection(wif.GetOutputPort())
        #png.SetFileName('img.png')
        #png.Write()
    
    def imagedata(self, orig, delta, dims):

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

        bounds = vtk.vtkOutlineFilter()
        bounds.SetInputData(self.grid)

        boundsMapper = vtk.vtkPolyDataMapper()
        boundsMapper.SetInputConnection(bounds.GetOutputPort())

        boundsActor = vtk.vtkActor()
        boundsActor.SetMapper(boundsMapper)
        boundsActor.GetProperty().SetColor(0, 0, 0)
        self.ren.AddActor(boundsActor)

    def lookuptable(self):

        #lutNCI = vtk.vtkLookupTable()
        #lutNCI.SetNumberOfColors(3)
        #lutNCI.SetHueRange(0.667, 0.0)
        #lutNCI.SetSaturationRange(1.0, 1.0)
        #lutNCI.SetValueRange(0.8, 0.8)
        #lutNCI.Build()

        # color palette
        # -3.0 -> blue
        #  0.1 -> green
        #  3.0 -> red
        colorNCI = vtk.vtkColorTransferFunction()
        colorNCI.AddRGBPoint(-3.0,0.0,0.0,1.0)
        colorNCI.AddRGBPoint(0.1,0.0,1.0,0.0)
        colorNCI.AddRGBPoint(3.0,1.0,0.0,0.0)
        self.colorNCI = colorNCI

    def scalarbar(self):

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

        log.debug('Marching Cubes algorithm contour filter')
        # ContourFilter or MarchingCubes
        #isoExtractor = vtk.vtkContourFilter()
        isoExtractor = vtk.vtkMarchingCubes()
        isoExtractor.SetInputData(self.grid)
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



