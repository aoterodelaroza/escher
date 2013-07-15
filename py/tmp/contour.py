
import vtk

 

# Prepare to read the file

readerVolume = vtk.vtkImageReader()
readerVolume.SetDataScalarType( vtk.VTK_UNSIGNED_SHORT )
readerVolume.SetFileDimensionality( 3 )
readerVolume.SetDataExtent ( 0,255, 0,255, 0,576)
readerVolume.SetDataSpacing( 1,1,1 )
readerVolume.SetNumberOfScalarComponents( 1 )
readerVolume.SetDataByteOrderToBigEndian()
readerVolume.SetFileName("/home/daniel/pkg/vis/VTKData/Data/HeadMRVolume.raw")

 
# Extract the region of interest
voiHead = vtk.vtkExtractVOI()
voiHead.SetInputConnection(readerVolume.GetOutputPort() )
voiHead.SetVOI( 0,255, 60,255, 0,100 )
voiHead.SetSampleRate( 1,1,1 )
 

# Generate an isosurface
# UNCOMMENT THE FOLLOWING LINE FOR CONTOUR FILTER
# contourBoneHead = vtk.vtkContourFilter()
contourBoneHead = vtk.vtkMarchingCubes()
contourBoneHead.SetInputConnection(voiHead.GetOutputPort() )
contourBoneHead.ComputeNormalsOn()
contourBoneHead.SetValue( 0, 120 )  # Bone isovalue
 

# Take the isosurface data and create geometry
geoBoneMapper = vtk.vtkPolyDataMapper()
geoBoneMapper.SetInputConnection(contourBoneHead.GetOutputPort() )
geoBoneMapper.ScalarVisibilityOff()
 

# Take the isosurface data and create geometry
actorBone = vtk.vtkLODActor()
actorBone.SetNumberOfCloudPoints( 1000000 )
actorBone.SetMapper( geoBoneMapper )
actorBone.GetProperty().SetColor( 1, 1, 1 )
 

# Create renderer
ren = vtk.vtkRenderer()
ren.SetBackground( 0.329412, 0.34902, 0.427451 ) #Paraview blue
ren.AddActor(actorBone)
 

# Create a window for the renderer of size 250x250
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(250, 250)
 

# Set an user interface interactor for the render window
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
 

# Start the initialization and rendering
iren.Initialize()
renWin.Render()
iren.Start()
