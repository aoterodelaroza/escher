# functions for reading positions, radii, and connections from file
# written by Filip Malmberg and modified by Erik Vidholm

import os
import string
import vtk

# read coordinates
def readPoints(filename):
    points = vtk.vtkPoints()
    file = open(filename)

    line = file.readline()
    while line:
        data = string.split(line)
      
        if data and data[0] != '#':
            x, y, z = float(data[0]), float(data[1]), float(data[2])
            points.InsertNextPoint(x, y, z)
        line = file.readline()

    return points;

# read scalars (radii)
def readScalars(filename):
    scalars = vtk.vtkFloatArray()
    file = open(filename)

    line = file.readline()
    while line:
        data = string.split(line)
        if data and data[0] != '#':
            x= float(data[0])
            scalars.InsertNextValue(x)
        line = file.readline()

    return scalars;

# read connections
def readConnections(filename):
    connections=vtk.vtkCellArray()
    file = open(filename)

    line = file.readline()
    while line:
        data = string.split(line)
        if data and data[0] != '#':
            a, b = int(data[0]), int(data[1])
            connections.InsertNextCell(2)
            connections.InsertCellPoint(a)
            connections.InsertCellPoint(b)
        line = file.readline()
        
    return connections

    

