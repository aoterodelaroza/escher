#!/usr/bin/env python

##
# This is an example of how to use the vtkVolumePicker.
##

import math
import vtk
from vtk.util.misc import vtkGetDataRoot
VTK_DATA_ROOT = vtkGetDataRoot()

from molecule import Molecule


a = Molecule()
a.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')
#a.readcube('/home/daniel/pkg/vis/VTKData/Data/m4_TotalDensity.cube')
#a.readcube('/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.5.cube')

#---------------------------------------------------------
# renderer and interactor
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(640, 480)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

#---------------------------------------------------------
## read the volume
#reader = vtk.vtkImageReader2()
#reader.SetDataExtent(0,63,0,63,0,92)
#reader.SetFileNameSliceOffset(1)
#reader.SetDataScalarTypeToUnsignedShort()
#reader.SetDataByteOrderToLittleEndian()
#reader.SetFilePrefix(str(VTK_DATA_ROOT) + "/Data/headsq/quarter")
#reader.SetDataSpacing(3.2,3.2,1.5)


# construct the volume
dx = 2.0
grid = vtk.vtkStructuredPoints()
#grid.SetOrigin(0, 0, 0) # default values
#grid.SetDimensions(5, 8, 10) # number of points in each direction
#grid.SetOrigin(-6.512752, -6.512752, -7.752607) # default values
#grid.SetSpacing(0.33, 0.33, 0.33)
#grid.SetDimensions(40, 40, 48) # number of points in each direction
grid.SetOrigin(a.x0) # default values
grid.SetSpacing(a.dx[0,0],a.dx[1,1],a.dx[2,2])
grid.SetDimensions(a.n) # number of points in each direction
print grid.GetNumberOfPoints()
# print grid.GetNumberOfCells()
array = vtk.vtkDoubleArray()
array.SetNumberOfComponents(1) # this is 3 for a vector
array.SetNumberOfTuples(grid.GetNumberOfPoints())
#for i in range(grid.GetNumberOfPoints()):
#    array.SetValue(i, a.alldata[i])
nd = 0
for k in range(a.n[2]):
    for j in range(a.n[1]):
        for i in range(a.n[0]):
            array.SetValue(nd, a.f[i,j,k])
            nd = nd + 1


array.SetName("density")
grid.GetPointData().SetScalars(array)
# print grid.GetPointData().GetNumberOfArrays()

reader = vtk.vtkGaussianCubeReader()
reader.SetFileName('../../escher_data/mol/ethylene_iso/c2h4.cube')

reader.Update()


#---------------------------------------------------------
# Do the surface rendering
boneExtractor = vtk.vtkMarchingCubes()
#boneExtractor.SetInputConnection(reader.GetOutputPort())
boneExtractor.SetInputData(grid)
boneExtractor.SetValue(0,0.1)

boneNormals = vtk.vtkPolyDataNormals()
boneNormals.SetInputConnection(boneExtractor.GetOutputPort())
boneNormals.SetFeatureAngle(35.0)

boneStripper = vtk.vtkStripper()
boneStripper.SetInputConnection(boneNormals.GetOutputPort())

boneLocator = vtk.vtkCellLocator()
boneLocator.SetDataSet(boneExtractor.GetOutput())
boneLocator.LazyEvaluationOn()

boneMapper = vtk.vtkPolyDataMapper()
boneMapper.SetInputConnection(boneStripper.GetOutputPort())
boneMapper.ScalarVisibilityOff()

boneProperty = vtk.vtkProperty()
boneProperty.SetColor(1.0,0.0,0.0)

bone = vtk.vtkActor()
bone.SetMapper(boneMapper)
bone.SetProperty(boneProperty)
bone.GetProperty().SetOpacity(.4)


#---------------------------------------------------------
ren.AddViewProp(bone)


#renWin.Render()
iren.Render()
iren.Start()

pov = vtk.vtkPOVExporter()
pov.SetRenderWindow(renWin)
pov.SetFileName('iso.pov')
pov.Write()
