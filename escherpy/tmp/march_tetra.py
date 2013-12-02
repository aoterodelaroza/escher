class Vector: # struct XYZ
 def __init__(self,x,y,z):
  self.x=x
  self.y=y
  self.z=z
 
 def __str__(self):
  return str(self.x)+" "+str(self.y)+" "+str(self.z)
  
class Gridcell: # struct GRIDCELL
 def __init__(self,p,n,val):
  self.p   = p   # p=[8]
  self.n   = n   # n=[8]
  self.val = val # val=[8]

class Triangle: # struct TRIANGLE
 def __init__(self,p1,p2,p3):
  self.p = [p1, p2, p3] # vertices

# return triangle as an ascii STL facet  
 def __str__(self):
  return """facet normal 0 0 0
outer loop
vertex %s
vertex %s
vertex %s
endloop
endfacet"""%(self.p[0],self.p[1],self.p[2])

# return a 3d list of values
def readdata(f=lambda x,y,z:x*x+y*y+z*z,size=5.0,steps=11):
 m=int(steps/2)
 ki = []
 for i in range(steps):
  kj = []
  for j in range(steps):
   kd=[]
   for k in range(steps):
    kd.append(f(size*(i-m)/m,size*(j-m)/m,size*(k-m)/m))
   kj.append(kd)
  ki.append(kj)
 return ki

from math import cos,exp,atan2

def lobes(x,y,z):
 try:
  theta = atan2(x,y)         # sin t = o 
 except:
  theta = 0
 try:
  phi = atan2(z,y)
 except:
  phi = 0
 r = x*x+y*y+z*z
 ct=cos(theta)
 cp=cos(phi)
 return ct*ct*cp*cp*exp(-r/10)
 
def main():

 data = readdata(lobes,5,41)
 isolevel = 0.1

 #print(data)
 
 triangles=[]
 for i in range(len(data)-1):
  for j in range(len(data[i])-1):
   for k in range(len(data[i][j])-1):
    p=[None]*8
    val=[None]*8
    #print(i,j,k)
    p[0]=Vector(i,j,k)
    val[0] = data[i][j][k]
    p[1]=Vector(i+1,j,k)
    val[1] = data[i+1][j][k]
    p[2]=Vector(i+1,j+1,k)
    val[2] = data[i+1][j+1][k]
    p[3]=Vector(i,j+1,k)
    val[3] = data[i][j+1][k]
    p[4]=Vector(i,j,k+1)
    val[4] = data[i][j][k+1]
    p[5]=Vector(i+1,j,k+1)
    val[5] = data[i+1][j][k+1]
    p[6]=Vector(i+1,j+1,k+1)
    val[6] = data[i+1][j+1][k+1]
    p[7]=Vector(i,j+1,k+1)
    val[7] = data[i][j+1][k+1]
   
    grid=Gridcell(p,[],val)
    triangles.extend(PolygoniseTri(grid,isolevel,0,2,3,7))
    triangles.extend(PolygoniseTri(grid,isolevel,0,2,6,7))
    triangles.extend(PolygoniseTri(grid,isolevel,0,4,6,7))
    triangles.extend(PolygoniseTri(grid,isolevel,0,6,1,2))
    triangles.extend(PolygoniseTri(grid,isolevel,0,6,1,4))
    triangles.extend(PolygoniseTri(grid,isolevel,5,6,1,4))
    
 export_triangles(triangles)

def export_triangles(triangles): # stl format
 print("solid points")
 for tri in triangles:
  print(tri)
 print("endsolid points")
 
def t000F(g, iso, v0, v1, v2, v3):
 return []

def t0E01(g, iso, v0, v1, v2, v3):
 return [Triangle(
 VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]),
 VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]),
 VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]))
 ]

def t0D02(g, iso, v0, v1, v2, v3):
 return [Triangle(
 VertexInterp(iso,g.p[v1],g.p[v0],g.val[v1],g.val[v0]),
 VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]),
 VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]))
 ]

def t0C03(g, iso, v0, v1, v2, v3):
 tri=Triangle(
 VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]),
 VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]),
 VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]))
 return [tri,Triangle(
 tri.p[2],
 VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]),
 tri.p[1])
 ]

def t0B04(g, iso, v0, v1, v2, v3):
 return [Triangle(
 VertexInterp(iso,g.p[v2],g.p[v0],g.val[v2],g.val[v0]),
 VertexInterp(iso,g.p[v2],g.p[v1],g.val[v2],g.val[v1]),
 VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]))
 ]

def t0A05(g, iso, v0, v1, v2, v3):
 tri = Triangle(
 VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]),
 VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]),
 VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]))
 return [tri,Triangle(
 tri.p[0],
 VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]),
 tri.p[1])
 ]

def t0906(g, iso, v0, v1, v2, v3):
 tri=Triangle(
 VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]),
 VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]),
 VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]))
 return [tri,
 Triangle(
 tri.p[0],
 VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]),
 tri.p[2])
 ]

def t0708(g, iso, v0, v1, v2, v3):
 return [Triangle(
 VertexInterp(iso,g.p[v3],g.p[v0],g.val[v3],g.val[v0]),
 VertexInterp(iso,g.p[v3],g.p[v2],g.val[v3],g.val[v2]),
 VertexInterp(iso,g.p[v3],g.p[v1],g.val[v3],g.val[v1]))
 ]

trianglefs = {7:t0708,8:t0708,9:t0906,6:t0906,10:t0A05,5:t0A05,11:t0B04,4:t0B04,12:t0C03,3:t0C03,13:t0D02,2:t0D02,14:t0E01,1:t0E01,0:t000F,15:t000F}

def PolygoniseTri(g, iso, v0, v1, v2, v3):
 triangles = []

 #   Determine which of the 16 cases we have given which vertices
 #   are above or below the isosurface

 triindex = 0;
 if g.val[v0] < iso: triindex |= 1
 if g.val[v1] < iso: triindex |= 2
 if g.val[v2] < iso: triindex |= 4
 if g.val[v3] < iso: triindex |= 8

 return trianglefs[triindex](g, iso, v0, v1, v2, v3)

def VertexInterp(isolevel,p1,p2,valp1,valp2):
 if abs(isolevel-valp1) < 0.00001 :
  return(p1);
 if abs(isolevel-valp2) < 0.00001 :
  return(p2);
 if abs(valp1-valp2) < 0.00001 :
  return(p1);
 mu = (isolevel - valp1) / (valp2 - valp1)
 return Vector(p1.x + mu * (p2.x - p1.x), p1.y + mu * (p2.y - p1.y), p1.z + mu * (p2.z - p1.z))

if __name__ == "__main__":
 main()
 
