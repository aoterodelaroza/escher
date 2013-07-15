#!/usr/bin/env python
import os

from vtk import *
from vtk.util.colors import *

ren = vtkRenderer()
renWin = vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# create pipeline
#
v16 = vtkVolume16Reader()
v16.SetDataDimensions(256,256)
v16.GetOutput().SetOrigin(0.0,0.0,0.0)
v16.SetFilePrefix(
  os.environ['HOME'] + "/python/examples/vtk/images/r")
v16.SetFilePattern( '%s%d.ima')
v16.SetDataByteOrderToBigEndian()
v16.SetImageRange(1001,1060)
v16.SetDataSpacing(1.0,1.0,3.5)
v16.Update()

#vtkImageMarchingCubes iso
iso = vtkMarchingCubes()
iso.SetInput(v16.GetOutput())
iso.SetValue(0,30)
#120 vessles near cerebellum
#100 cortex
#20 face
#iso SetStartMethod {puts "Start Marching"}



isoMapper = vtkPolyDataMapper()
isoMapper.SetInput(iso.GetOutput())
isoMapper.ScalarVisibilityOff()

isoActor = vtkActor()
isoActor.SetMapper(isoMapper)
isoActor.GetProperty().SetColor(antique_white)

outline = vtkOutlineFilter()
outline.SetInput(v16.GetOutput())
outlineMapper = vtkPolyDataMapper()
outlineMapper.SetInput(outline.GetOutput())
outlineActor = vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.VisibilityOff()

# Add the actors to the renderer, set the background and size
#
ren.AddActor(outlineActor)
ren.AddActor(isoActor)
ren.SetBackground(0.2,0.3,0.4)
renWin.SetSize(450,450)
## ren.GetActiveCamera().Elevation(235)
## ren.GetActiveCamera().SetViewUp(0,.5,-1)
## ren.GetActiveCamera().Azimuth(90)


iren.Initialize()

iren.Start()
