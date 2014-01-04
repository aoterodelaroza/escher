
blenderinput = \
'''
import os
import bpy

bpy.data.objects['Cube'].select = True
bpy.ops.object.delete()
bpy.data.objects['Camera'].select = True
bpy.ops.object.delete()

bpy.ops.import_scene.x3d(filepath=os.getcwd() + "/out.wrl")

cam = bpy.data.cameras.get('Viewpoint')
cam.type = 'ORTHO'
cam.ortho_scale = 20.
'''
