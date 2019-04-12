import os


with open("isosurface.1.node") as myfile:
    lines = [list(filter(len,line.strip().split(' '))) for line in myfile]

numnodes = int(lines[0][0] )

nodes=[]
for iii in  range(numnodes):
   nodes.append( (float(lines[iii+1][1]),float(lines[iii+1][2]),float(lines[iii+1][3])) )

#image data
dim = [110, 165, 5] 
dim = [20, 20, 5] 
bb = [[49.217,10.5465,0], [126.558,126.558,25]]
bb = [ [1.e-3*bb[0][0] , 1.e-3*bb[0][1] , 1.e-3*bb[0][2]], [1.e-3*bb[1][0] , 1.e-3*bb[1][1] , 1.e-3*bb[1][2]] ]
vox = [0.7031, 0.7031, 5]

for iii in range(dim[0]):
  for jjj in range(dim[1]):
     for kkk in range(dim[2]):
         newnode = (bb[0][0] + (bb[1][0]-bb[0][0])/dim[0] * iii,
                    bb[0][1] + (bb[1][1]-bb[0][1])/dim[1] * jjj,
                    bb[0][2] + (bb[1][2]-bb[0][2])/dim[2] * kkk)
         nodes.append(newnode)

f = open("mytest.node", "w")
f.write("%d 3 0 0\n" % len(nodes))
for iii, node in enumerate(nodes):
  f.write("%d %f %f %f\n" % (iii, node[0],node[1],node[2]))
f.close()

os.system("tetgen -k mytest.node")
