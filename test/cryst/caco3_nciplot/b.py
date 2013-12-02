import bpy
from mathutils import Vector

def look_at(obj_camera, point):
    loc_camera = obj_camera.matrix_world.to_translation()

    direction = point - loc_camera
    # point the cameras '-Z' and use its 'Y' as up
    rot_quat = direction.to_track_quat('-Z', 'Y')

    # assume we're using euler rotation
    obj_camera.rotation_euler = rot_quat.to_euler()

# Test
obj_camera = bpy.data.objects["Camera"]
obj_other = bpy.data.objects["Cube"]

obj_camera.location = (5.0, 2.0, 3.0)
look_at(obj_camera, Vector((0.,0.,0.)))
