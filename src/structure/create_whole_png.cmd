#feedback disable,all,actions
#feedback disable,all,results
#feedback disable,all,blather
#feedback disable,all,debugging
#feedback disable,all,details
#feedback disable,all,warnings
#feedback disable,all,errors

delete all
reset

load STRUCT

viewport 350,350

hide

# here, a section to load PSSM sec struct instead of lousy Pymol one

bg_color white
util.cbc
color red, (resn DUM and name O)
color blue, (resn DUM and name N)

set antialias=1
set line_smooth=1.00000
set depth_cue=1
set specular=1.00000
set surface_quality=1.00000
set stick_quality=25.00000
set sphere_quality=2.00000
set cartoon_sampling=14.00000
set ribbon_sampling=10.00000
set ray_trace_fog=1

rotate [1,0,0], angle=-90
show cartoon
set sphere_transparency=0.75
show spheres, (resname DUM)


png f0001.png

quit
