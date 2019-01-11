# cmd line
# $ /opt/apps/cubit/v13.2/bin/claro -nographics -batch python meshVessel.py
# 
# from cubit GUI cmd script editor
#   exec(file("file.py"))
cubit.cmd('set developer on')
cubit.cmd('   reset')

exitradius = 3.

# bounding box for ROI (mm)'
cubit.cmd('create brick  width 70 height 80 depth 70')
cubit.cmd('volume 1 name "healthy" ')
idhealthy = cubit.get_id_from_name('healthy')
print "id" ,idhealthy

# taper vessel
cubit.cmd('create frustum height 40 radius 5 top %f' % exitradius)
cubit.cmd('volume 2 name "vessel" ')
idvessel = cubit.get_id_from_name('vessel')
print "id" ,idvessel
cubit.cmd('volume %d  move z -20' % idvessel)

# exit driven by pressure at the surface
cubit.cmd('create sphere radius  %f' % exitradius)
cubit.cmd('volume 3 name "vesselexit" ')

cubit.cmd('webcut volume 1 tool volume 2')
cubit.cmd('webcut volume 1 tool volume 3')
cubit.cmd('delete volume 2 3')

# merge geometry
cubit.cmd('imprint volume 1 4 5')
cubit.cmd('merge volume   1 4 5')


# mesh
cubit.cmd('volume 1 4 5 scheme tetmesh' )
cubit.cmd('volume      1     size 4.' )
cubit.cmd('volume        4 5 size 2.' )
cubit.cmd('mesh volume 1 4 5')

# BC
cubit.cmd('sideset 2 surface 17')
#cubit.cmd('nodeset 3 node in surface 11')
cubit.cmd('nodeset 3 node in surface 16')
#cubit.cmd('nodeset 2 name "marker"') TODO: name note used in DMPlex

cubit.cmd('skin volume 1 make group 7')
cubit.cmd('group 7 name "appface"')
cubit.cmd('nodeset 4 node in tri in group 7')
cubit.cmd('nodeset 4 node in surface 17  remove')

cubit.cmd('export mesh "meshVessel.e" overwrite' )

## cubit.cmd('webcut volume %d with plane xplane' % idhealthy)
## cubit.cmd('webcut volume %d with plane yplane' % idhealthy)
## #   
## # a single probe containing
## # both charged surfaces within a single needle of diameter
## # 1.65mm, with distal and proximal conducting surfaces of
## # 0.70 and 0.75cm long, respectively, separated by a 0.80 cm
## # long insulation area. These
## # applicator consists of two electrodes 7mm in length placed at 8cm apart.
## #   
## #    total length = 22mm
## ActiveDistanceFromTip = 3.0 # mm
## ActiveLength = 11.0 # mm
## cubit.cmd('   webcut volume %d with plane zplane offset %f  ' % (idhealthy,-ActiveLength-ActiveDistanceFromTip) )
## cubit.cmd('   webcut volume %d with plane zplane offset %f  ' % (idhealthy,-ActiveLength) )
## cubit.cmd('   webcut volume %d with plane zplane offset %f  ' % (idhealthy, ActiveLength) )
## #   
## # applicator is a 17 gauge needle = 1.15mm diameter
## # 
## outerRadius  = 1.65 / 2.0 #mm
## print "fiber radius = %f " % ( outerRadius )
## # create outer shell       
## cubit.cmd('create cylinder radius %f height 100' % outerRadius)
## cubit.cmd('webcut volume 1 4 5 6 tool volume 7 ')
## cubit.cmd('delete volume 7 8 10 11')
## 
## resolutionID = 'lores'
## resolutionID = 'midres'
## meshResolutions = {'lores':(1.0 , 5.8 ), 'midres':(0.6 , 3.0 ), 'hires':(0.25, 2.5 )}
## (fineSize,coarseSize) = meshResolutions[resolutionID]
## 
## # mesh 
## cubit.cmd('volume 1 4 5 6 9 scheme tetmesh' )
## cubit.cmd('volume 1 4 5 6 9 size %f' % coarseSize )
## curveBiasTemplate = 'curve %d scheme bias fine size %f coarse size %f start vertex %d'
## #volume 1 radial bias
## cubit.cmd( curveBiasTemplate  %(101,fineSize,coarseSize, 51))
## cubit.cmd( curveBiasTemplate  %(102,fineSize,coarseSize, 52))
## cubit.cmd( curveBiasTemplate  %( 99,fineSize,coarseSize, 53))
## cubit.cmd( curveBiasTemplate  %(100,fineSize,coarseSize, 54))
## #volume 1 axial bias
## cubit.cmd( curveBiasTemplate  %( 98,fineSize,coarseSize, 52))
## cubit.cmd( curveBiasTemplate  %( 96,fineSize,coarseSize, 51))
## cubit.cmd('mesh volume 1')
## # volume 6 radial bias
## cubit.cmd( curveBiasTemplate  %( 131,fineSize,coarseSize, 69))
## cubit.cmd( curveBiasTemplate  %( 132,fineSize,coarseSize, 70))
## # volume 6 axial equal 
## cubit.cmd('curve 144 146 size %f' % fineSize )
## cubit.cmd('mesh volume 6')
## # volume 5 radial bias
## cubit.cmd( curveBiasTemplate  %( 117,fineSize,coarseSize, 59))
## cubit.cmd( curveBiasTemplate  %( 118,fineSize,coarseSize, 60))
## cubit.cmd('mesh volume 5')
## # volume 4 radial bias
## cubit.cmd( curveBiasTemplate  %( 116,fineSize,coarseSize, 62))
## cubit.cmd( curveBiasTemplate  %( 115,fineSize,coarseSize, 61))
## # volume 4 axial bias
## ## cubit.cmd( curveBiasTemplate  %( 59,fineSize,coarseSize, 25))
## ## cubit.cmd( curveBiasTemplate  %( 60,fineSize,coarseSize, 27))
## ## cubit.cmd( curveBiasTemplate  %( 58,fineSize,coarseSize, 28))
## cubit.cmd( curveBiasTemplate  %(114,fineSize,coarseSize, 60))
## cubit.cmd( curveBiasTemplate  %(112,fineSize,coarseSize, 59))
## cubit.cmd('mesh volume 4')
## # volume 9 axial bias
## cubit.cmd( curveBiasTemplate  %(57 ,fineSize,coarseSize, 30))
## cubit.cmd('mesh volume 9')
## #   
## cubit.cmd('group "badhex"  equals quality hex  in volume  all   jacobian high 0.0')
## cubit.cmd('group "badtet"  equals quality tet  in volume  all   jacobian high 0.0')
## # export mesh in distinct pieces
## cubit.cmd('reset genesis')
## 
## # reflect to get full mesh
## cubit.cmd('volume all copy reflect x')
## cubit.cmd('volume all copy reflect y')
## cubit.cmd('imprint volume all')
## cubit.cmd('merge volume all')
## ## 
## ## 
## ## ## # create second cylinder for mesh indepenent nodeset
## ## ## # offset by epsilon to avoid divide by zero in source evaluation
## ## ## epsilon = 0.05
## ## ## cubit.cmd('create cylinder radius %f height 10' % ( outerRadius - epsilon ))
## ## ## cubit.cmd('volume      70 size 0.1' )
## ## ## cubit.cmd('mesh volume      70 ' )
## ## 
## ## ## # bc
## ## cubit.cmd('skin volume all make sideset 4')
## ## cubit.cmd('sideset 4 name "neumann" ')
## ## ## cubit.cmd('sideset 2 face in surface 48 49 61 35 36 71 43 44 55 56 remove')
## ## ## cubit.cmd('sideset 3         surface 48 49 61 35 36 71 43 44 55 56')
## ## ## cubit.cmd('sideset 3 name "cauchy" ')
## ## ## #   ##sideset 3 surface 48 49 61 35 36 71 43 44 55 56 69 79
## ## ## #   # specify nodeset for applicator
## ## ## #   # skin the applicator to remove the boundary from the nodeset
## ## ## #   skin volume 8 10 11 16 18 19 24 26 27 32 34 35 make group 3
## ## ## #   group 3 name "appface"
## ## ## #   group "appnode" add node in face in group 3
## ## cubit.cmd('nodeset 1 volume 8 10 11 13 15 16 18 20 26 28 29 30 32 33 34 36 42 44 45 46 48 49 50 52 58 60 61 62 64 65 66 68 21 37 53 69')
## ## cubit.cmd('nodeset 1 name "dirichlet"')
## # mesh dependent node set
## cubit.cmd('skin volume all make group 7')
## cubit.cmd('group 7 name "appface"')
## cubit.cmd('nodeset 4 node in tri in group 7')
## #cubit.cmd('nodeset 4 node in surface 60 100 129 158 172 143 114 80 90 149 178 120 remove')
## cubit.cmd('nodeset 4 node in surface 90 149 178 120 172 143 114 80 remove')
## ## specify active electrodes by location
## ## ground electrode
## cubit.cmd('nodeset 3 node with z_coord >  4 in surface 172 143 114 80 90 149 178 120')
## ## active electrode
## cubit.cmd('nodeset 2 node with z_coord < -4 in surface 172 143 114 80 90 149 178 120')
## ## rest of applicator is insulated
## cubit.cmd('nodeset 4 node with z_coord > -4 and z_coord <  4 in surface 172 143 114 80 90 149 178 120')
## #cubit.cmd('nodeset 2 node in surface 90 149 178 120')
## 
## #cubit.cmd('nodeset 2 name "marker"') TODO: name note used in DMPlex
## ## # mesh independent node set
## ## cubit.cmd('nodeset 3 surface 402 ')
## ## cubit.cmd('nodeset 3 name "IndependentSources"')
## ## ## #   #nodeset 1 node in group 4 remove
## ## ## #   nodeset 2 node in group 4
## ## ## #   nodeset 2 surface 60 80 90
## ## ## #   nodeset 2 name "appboundarynode"
## # volume
## cubit.cmd('block 1 volume all ')
## cubit.cmd('block 1 name "tissue"  ')
## ## cubit.cmd('block 2 name "catheder"  ')
## ## cubit.cmd('block 3 volume 11 29 45 61')
## ## cubit.cmd('block 3 name "laserTip"  ')
## ## ##cubit.cmd('block 4 volume 70 ')
## ## ##cubit.cmd('block 4 name "laserSources"  ')
## cubit.cmd('volume all scale 0.001')
## cubit.cmd('export mesh "meshIRE%s.e" overwrite' % resolutionID)
## ## ##cubit.cmd('export abaqus "meshIRE%s.inp" overwrite' % resolutionID)
