#!/usr/bin/python

import vtk

# three 'ones' in their common face
st = vtk.vtkStructuredPoints()
st.SetDimensions(2,3,2)
st.AllocateScalars()
sc = (st.GetPointData()).GetScalars()
sc.SetValue(  0, 0.0 )
sc.SetValue(  1, 0.0 )
sc.SetValue(  2, 1.0 )
sc.SetValue(  3, 1.0 )
sc.SetValue(  4, 0.0 )
sc.SetValue(  5, 0.0 )
sc.SetValue(  6, 0.0 )
sc.SetValue(  7, 0.0 )
sc.SetValue(  8, 0.0 )
sc.SetValue(  9, 1.0 )
sc.SetValue( 10, 0.0 )
sc.SetValue( 11, 0.0 )

mc = vtk.vtkMarchingCubes()
mc.SetNumberOfContours(1)
mc.SetValue( 0, 1.0 )
mc.SetInput(st)
mc.Update()

pd = mc.GetOutput()
print "Number of Cells: %d\n" % pd.GetNumberOfCells() # return 2L
print pd.GetCell(0) # return a triangle, coincident with the face
print pd.GetCell(1) # return another triangle, in the same position
