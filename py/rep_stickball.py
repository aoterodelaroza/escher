#!/usr/bin/env python

from vtk import *
from math import cos, sin, pi, sqrt, acos, degrees
from molecule import Molecule
###### Functions for making objects ################


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
    sphere.SetCenter(position[0], position[1], position[2])
    # Convert the sphere into polygons
    sphereMapper = vtk.vtkPolyDataMapper()
    sphereMapper.SetInputConnection(sphere.GetOutputPort())

    #Create an actor for the sphere 
    sphereActor = vtk.vtkActor()
    sphereActor.SetMapper(sphereMapper)
    (sphereActor.GetProperty()).SetColor(color[0], color[1], color[2])
    return sphereActor

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

    # create cone
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


###### The Main Program Starts Here ################

# camera (optional)
camera = vtk.vtkCamera ();
camera.SetPosition(10, 0,10);
camera.SetFocalPoint(0, 0, 0);

# create a rendering window and renderer
ren = vtkRenderer()
# camera (optional)
ren.SetActiveCamera(camera);
renWin = vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(1700,1200)
renWin.StereoCapableWindowOn()

iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)





a = Molecule()
a.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')
#a.readcube('/home/daniel/pkg/vis/VTKData/Data/m4_TotalDensity.cube')


for i in range(len(a.atxyz)):
    if a.atname[i] == 'H':
        color = [1,1,1]
    if a.atname[i] == 'C':
        color = [0.5,0.5,0.5]
    if a.atname[i] == 'N':
        color = [0,0,1,255]
    if a.atname[i] == 'O':
        color = [1,0,0,255]
    ren.AddActor(ball(0.6, a.atxyz[i], color))

for i in range(len(a.atxyz)):
    for j in range(i):
        dist = sqrt((a.atxyz[i][0] - a.atxyz[j][0])**2 + (a.atxyz[i][1] - a.atxyz[j][1])**2 + (a.atxyz[i][2] - a.atxyz[j][2])**2)
        if (dist) < 3. :
            ren.AddActor(stick(a.atxyz[i], a.atxyz[j]))







# enable user interface interactor
iren.Initialize()


iren.Start()

pov = vtk.vtkPOVExporter()
pov.SetRenderWindow(renWin)
pov.SetFileName('out.pov')
pov.Write()

## image
#wif = vtk.vtkWindowToImageFilter()
#wif.SetInputData(renWin)
#png = vtk.vtkPNGWriter()
#png.SetInputConnection(wif.GetOutputPort())
#png.SetFileName('img.png')
#png.Write()

