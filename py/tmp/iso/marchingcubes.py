#!/usr/bin/env python

##
# This is an example of how to use the vtkVolumePicker.
##

import math
import vtk
from vtk.util.misc import vtkGetDataRoot
VTK_DATA_ROOT = vtkGetDataRoot()

#---------------------------------------------------------
# renderer and interactor
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

#---------------------------------------------------------
# read the volume
#reader = vtk.vtkImageReader2()
#reader.SetDataExtent(0,63,0,63,0,92)
#reader.SetFileNameSliceOffset(1)
#reader.SetDataScalarTypeToUnsignedShort()
#reader.SetDataByteOrderToLittleEndian()
#reader.SetFilePrefix(str(VTK_DATA_ROOT) + "/Data/headsq/quarter")
#reader.SetDataSpacing(3.2,3.2,1.5)

# image reader
filename = "hydrogen.vtk"
reader = vtk.vtkStructuredPointsReader()
reader.SetFileName( filename )
# must call Update() before we fetch the dimensions
reader.Update()


#---------------------------------------------------------
# Do the surface rendering
boneExtractor = vtk.vtkContourFilter()
boneExtractor.SetInputConnection(reader.GetOutputPort())
boneExtractor.SetValue(0,0.3)

boneNormals = vtk.vtkPolyDataNormals()
boneNormals.SetInputConnection(boneExtractor.GetOutputPort())
boneNormals.SetFeatureAngle(60.0)

boneStripper = vtk.vtkStripper()
boneStripper.SetInputConnection(boneNormals.GetOutputPort())

boneLocator = vtk.vtkCellLocator()
boneLocator.SetDataSet(boneExtractor.GetOutput())
boneLocator.LazyEvaluationOn()

boneMapper = vtk.vtkPolyDataMapper()
boneMapper.SetInputConnection(boneStripper.GetOutputPort())
boneMapper.ScalarVisibilityOff()

boneProperty = vtk.vtkProperty()
boneProperty.SetColor(1.0,1.0,0.9)

bone = vtk.vtkActor()
bone.SetMapper(boneMapper)
bone.SetProperty(boneProperty)


#---------------------------------------------------------
ren.AddViewProp(bone)


#renWin.Render()
iren.Render()
iren.Start()
