#!/usr/bin/env python

# This example demonstrates the use of the vtkTransformPolyDataFilter
# to reposition a 3D text string.

import vtk
from vtk.util.colors import *

def marching(values, coords):
    # Define a Single Cube

    Scalars = vtk.vtkFloatArray()
    for i in values:
        Scalars.InsertNextValue(i)

    Points = vtk.vtkPoints()
    for i in coords:
        Points.InsertNextPoint(i)

    Ids = vtk.vtkIdList()
    for i in range(8):
        Ids.InsertNextId(i)

    Grid = vtk.vtkUnstructuredGrid()
    Grid.Allocate(20, 20)
    Grid.InsertNextCell(12, Ids)
    Grid.SetPoints(Points)
    Grid.GetPointData().SetScalars(Scalars)

    # Find the triangles that lie along the 0.5 contour in this cube.
    Marching = vtk.vtkContourFilter()
    Marching.SetInputData(Grid)
    Marching.SetValue(0, 0.5)
    Marching.Update()

    # Extract the edges of the triangles just found.
    triangleEdges = vtk.vtkExtractEdges()
    triangleEdges.SetInputConnection(Marching.GetOutputPort())
    # Draw the edges as tubes instead of lines.  Also create the associated
    # mapper and actor to display the tubes.
    triangleEdgeTubes = vtk.vtkTubeFilter()
    triangleEdgeTubes.SetInputConnection(triangleEdges.GetOutputPort())
    triangleEdgeTubes.SetRadius(.005)
    triangleEdgeTubes.SetNumberOfSides(6)
    triangleEdgeTubes.UseDefaultNormalOn()
    triangleEdgeTubes.SetDefaultNormal(.577, .577, .577)
    triangleEdgeMapper = vtk.vtkPolyDataMapper()
    triangleEdgeMapper.SetInputConnection(triangleEdgeTubes.GetOutputPort())
    triangleEdgeMapper.ScalarVisibilityOff()
    triangleEdgeActor = vtk.vtkActor()
    triangleEdgeActor.SetMapper(triangleEdgeMapper)
    triangleEdgeActor.GetProperty().SetDiffuseColor(lamp_black)
    triangleEdgeActor.GetProperty().SetSpecular(.4)
    triangleEdgeActor.GetProperty().SetSpecularPower(10)

    # Shrink the triangles we found earlier.  Create the associated mapper
    # and actor.  Set the opacity of the shrunken triangles.
    aShrinker = vtk.vtkShrinkPolyData()
    aShrinker.SetShrinkFactor(1)
    aShrinker.SetInputConnection(Marching.GetOutputPort())
    aMapper = vtk.vtkPolyDataMapper()
    aMapper.ScalarVisibilityOff()
    aMapper.SetInputConnection(aShrinker.GetOutputPort())
    Triangles = vtk.vtkActor()
    Triangles.SetMapper(aMapper)
    Triangles.GetProperty().SetDiffuseColor(banana)
    Triangles.GetProperty().SetOpacity(.4)

    # Draw a cube the same size and at the same position as the one
    # created previously.  Extract the edges because we only want to see
    # the outline of the cube.  Pass the edges through a vtkTubeFilter so
    # they are displayed as tubes rather than lines.
    CubeModel = vtk.vtkCubeSource()
    CubeModel.SetCenter(.5, .5, .5)
    Edges = vtk.vtkExtractEdges()
    Edges.SetInputConnection(CubeModel.GetOutputPort())
    Tubes = vtk.vtkTubeFilter()
    Tubes.SetInputConnection(Edges.GetOutputPort())
    Tubes.SetRadius(.01)
    Tubes.SetNumberOfSides(6)
    Tubes.UseDefaultNormalOn()
    Tubes.SetDefaultNormal(.577, .577, .577)
    # Create the mapper and actor to display the cube edges.
    TubeMapper = vtk.vtkPolyDataMapper()
    TubeMapper.SetInputConnection(Tubes.GetOutputPort())
    CubeEdges = vtk.vtkActor()
    CubeEdges.SetMapper(TubeMapper)
    CubeEdges.GetProperty().SetDiffuseColor(khaki)
    CubeEdges.GetProperty().SetSpecular(.4)
    CubeEdges.GetProperty().SetSpecularPower(10)

    # Create a sphere to use as a glyph source for vtkGlyph3D.
    Sphere = vtk.vtkSphereSource()
    Sphere.SetRadius(0.04)
    Sphere.SetPhiResolution(20)
    Sphere.SetThetaResolution(20)
    # Remove the part of the cube with data values below 0.5.
    ThresholdIn = vtk.vtkThresholdPoints()
    ThresholdIn.SetInputData(Grid)
    ThresholdIn.ThresholdByUpper(.5)
    # Display spheres at the vertices remaining in the cube data set after
    # it was passed through vtkThresholdPoints.
    Vertices = vtk.vtkGlyph3D()
    Vertices.SetInputConnection(ThresholdIn.GetOutputPort())
    Vertices.SetSourceConnection(Sphere.GetOutputPort())
    # Create a mapper and actor to display the glyphs.
    SphereMapper = vtk.vtkPolyDataMapper()
    SphereMapper.SetInputConnection(Vertices.GetOutputPort())
    SphereMapper.ScalarVisibilityOff()
    CubeVertices = vtk.vtkActor()
    CubeVertices.SetMapper(SphereMapper)
    CubeVertices.GetProperty().SetDiffuseColor(tomato)
    CubeVertices.GetProperty().SetDiffuseColor(tomato)

    # Add the actors to the renderer
    ren.AddActor(triangleEdgeActor)
    ren.AddActor(CubeEdges)
    ren.AddActor(CubeVertices)
    ren.AddActor(Triangles)


# Create the Renderer, RenderWindow, and RenderWindowInteractor
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(640, 480)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

values = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0]
coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
marching(values, coords)
values = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0]
coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, -1], [1, 0, -1], [1, 1, -1], [0, 1, -1]]
marching(values, coords)
values = [0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0]
coords = [[0, 0, -1], [1, 0, -1], [1, 1, -1], [0, 1, -1], [0, 0, -2], [1, 0, -2], [1, 1, -2], [0, 1, -2]]
marching(values, coords)

# Set the background color.
ren.SetBackground(slate_grey)

# Position the camera.
ren.ResetCamera()
ren.GetActiveCamera().Dolly(1.2)
ren.GetActiveCamera().Azimuth(30)
ren.GetActiveCamera().Elevation(20)
ren.ResetCameraClippingRange()

iren.Initialize()
renWin.Render()
iren.Start()
