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

turn y, 0.03125
png f0001.png, ray=1
turn y, 0.0625
png f0002.png, ray=1
turn y, 0.125
png f0003.png, ray=1
turn y, 0.25
png f0004.png, ray=1
turn y, 0.25
png f0005.png, ray=1
turn y, 0.25
png f0006.png, ray=1
turn y, 0.25
png f0007.png, ray=1
turn y, 0.25
png f0008.png, ray=1
turn y, 0.25
png f0009.png, ray=1
turn y, 0.25
png f0010.png, ray=1
turn y, 0.25
png f0011.png, ray=1
turn y, 0.25
png f0012.png, ray=1
turn y, 0.25
png f0013.png, ray=1
turn y, 0.25
png f0014.png, ray=1
turn y, 0.25
png f0015.png, ray=1
turn y, 0.25
png f0016.png, ray=1
turn y, 0.25
png f0017.png, ray=1
turn y, 0.25
png f0018.png, ray=1
turn y, 0.25
png f0019.png, ray=1
turn y, 0.25
png f0020.png, ray=1
turn y, 0.25
png f0021.png, ray=1
turn y, 0.25
png f0022.png, ray=1
turn y, 0.25
png f0023.png, ray=1
turn y, 0.25
png f0024.png, ray=1
turn y, 0.25
png f0025.png, ray=1
turn y, 0.25
png f0026.png, ray=1
turn y, 0.25
png f0027.png, ray=1
turn y, 0.25
png f0028.png, ray=1
turn y, 0.25
png f0029.png, ray=1
turn y, 0.25
png f0030.png, ray=1
turn y, 0.25
png f0031.png, ray=1
turn y, 0.25
png f0032.png, ray=1
turn y, 0.25
png f0033.png, ray=1
turn y, 0.25
png f0034.png, ray=1
turn y, 0.25
png f0035.png, ray=1
turn y, 0.25
png f0036.png, ray=1
turn y, 0.25
png f0037.png, ray=1
turn y, 0.25
png f0038.png, ray=1
turn y, 0.25
png f0039.png, ray=1
turn y, 0.25
png f0040.png, ray=1
turn y, 0.25
png f0041.png, ray=1
turn y, 0.25
png f0042.png, ray=1
turn y, 0.25
png f0043.png, ray=1
turn y, 0.25
png f0044.png, ray=1
turn y, 0.25
png f0045.png, ray=1
turn y, 0.25
png f0046.png, ray=1
turn y, 0.25
png f0047.png, ray=1
turn y, 0.25
png f0048.png, ray=1
turn y, 0.125
png f0049.png, ray=1
turn y, 0.0625
png f0050.png, ray=1
turn y, 0.03125
png f0051.png, ray=1

quit
