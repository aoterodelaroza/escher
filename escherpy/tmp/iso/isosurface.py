# File:        isosurface.py
# Description: Iso-surface extraction from a volume
#              Run this example from a command prompt by typing: 
#              "python isosurface.py"

import vtk

# image reader
filename = "hydrogen.vtk"
reader = vtk.vtkStructuredPointsReader()
reader.SetFileName( filename )
 # must call Update() before we fetch the dimensions
reader.Update()

# just for illustration:
# get the extent of the data and print it
W,H,D = reader.GetOutput().GetDimensions()
# string formatting is similar to the sprintf style in C
print "Reading '%s', width=%i, height=%i, depth=%i" %(filename, W, H, D)

# create an outline of the dataset
outline = vtk.vtkOutlineFilter()
outline.SetInputConnection(reader.GetOutputPort())
outlineMapper = vtk.vtkPolyDataMapper()
outlineMapper.SetInputConnection(outline.GetOutputPort())
outlineActor = vtk.vtkActor()
outlineActor.SetMapper( outlineMapper )

# the actors property defines color, shading, line width,...
outlineActor.GetProperty().SetColor(0.0,0.0,1.0)
outlineActor.GetProperty().SetLineWidth(2.0)

#
#
#
#
# add your code here...
iso = vtk.vtkContourFilter()
iso.SetInputConnection(reader.GetOutputPort())
iso.SetValue(0,0.01)

normals = vtk.vtkPolyDataNormals()
normals.SetInputConnection(iso.GetOutputPort())
normals.SetFeatureAngle(45)
#
#
#
#

# renderer and render window 
ren = vtk.vtkRenderer()
ren.SetBackground(.8, .8, .8)
renWin = vtk.vtkRenderWindow()
renWin.SetSize( 400, 400 )
renWin.AddRenderer( ren )

# render window interactor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow( renWin )

# add the actors
ren.AddActor( outlineActor )
renWin.Render()

## create window to image filter to get the window to an image
w2if = vtk.vtkWindowToImageFilter()
w2if.SetInput(renWin)

## create png writer
wr = vtk.vtkPNGWriter()
wr.SetInputConnection(w2if.GetOutputPort())

# Python function for the keyboard interface
# count is a screenshot counter
count = 0
def Keypress(obj, event):
    global count, iv
    key = obj.GetKeySym()
    if key == "s":
        renWin.Render()     
        w2if.Modified() # tell the w2if that it should update
        fnm = "screenshot%02d.png" %(count)
        wr.SetFileName(fnm)
        wr.Write()
        print "Saved '%s'" %(fnm)
        count = count+1
    # add your keyboard interface here
    # elif key == ...

# add keyboard interface, initialize, and start the interactor
iren.AddObserver("KeyPressEvent", Keypress)
iren.Initialize()
iren.Start()

