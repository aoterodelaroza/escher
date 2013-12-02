# File: molecules.py

import vtk

# import the functions in ReadPoints.py
from ReadPoints import *

#Read the data into a vtkPolyData using the functions in ReadPoints.py
data = vtk.vtkPolyData()
data.SetPoints( readPoints("coordinates.txt") )
data.GetPointData().SetScalars( readScalars("radii.txt") )
data.SetLines( readConnections("connections.txt") )














# renderer and render window 
ren = vtk.vtkRenderer()
ren.SetBackground(.2, .2, .2)
renWin = vtk.vtkRenderWindow()
renWin.SetSize(400, 400)
renWin.AddRenderer( ren )

# render window interactor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow( renWin )

# add the actors to the renderer
# ren.AddActor( ... )

# render
renWin.Render()

## create window to image filter to get the window to an image
#w2if = vtk.vtkWindowToImageFilter()
#w2if.SetInput(renWin)

# create png writer
#wr = vtk.vtkPNGWriter()
#wr.SetInput(w2if.GetOutput())

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
