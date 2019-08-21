#!/usr/bin/env python 

## \file 3DChannelPeriodic.py
#  \brief Python script for box meshing with periodicity
#  \author T. A. Oliver
#  \version 3.2.8 "eagle"
#
# SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2015 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from optparse import OptionParser

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename", default="channel.su2",
                  help="write mesh to FILE", metavar="FILE")
parser.add_option("-n", "--nNode", dest="nNode", default=5,
                  help="use this NNODE in x direction", metavar="NNODE")
parser.add_option("-m", "--mNode", dest="mNode", default=5,
                  help="use this MNODE in y direction", metavar="MNODE")
parser.add_option("-l", "--lNode", dest="lNode", default=5,
                  help="use this LNODE in z direction", metavar="LNODE")
parser.add_option("-x", "--xLength", dest="xLength", default=1.0,
                  help="use this XLENGTH", metavar="XLENGTH")
parser.add_option("-y", "--yLength", dest="yLength", default=1.0,
                  help="use this YLENGTH", metavar="YLENGTH")
parser.add_option("-z", "--zLength", dest="zLength", default=1.0,
                  help="use this ZLENGTH", metavar="ZLENGTH")
parser.add_option("-d", "--delta", dest="delta", default=0.0,
                  help="delta in tanh spacing in z dir", metavar="DELTA")
parser.add_option("-s", "--deltaZ", dest="deltaZ", default=-1.0,
                  help="dz at the wall", metavar="DELTAZ")
(options, args)=parser.parse_args()

KindElem = 12
KindBound = 9
nNode = int(options.nNode)
mNode = int(options.mNode)
lNode = int(options.lNode)
xLength = float(options.xLength)
yLength = float(options.yLength)
zLength = float(options.zLength)
delta = float(options.delta)
dZ = float(options.deltaZ)

# Add periodic halo
Lx = xLength
Ly = yLength

dx = Lx/(nNode-1)
dy = Ly/(mNode-1)

xLength = Lx+dx
yLength = Ly+dy

nNode = nNode + 1
mNode = mNode + 1

# vertical spacing
zz = np.linspace(0, zLength, lNode)
if (delta > 0.0):
    ztmp = np.linspace(0,1,lNode);
    f2 = 0.5*(1.0 + np.tanh((ztmp-0.5)*delta)/np.tanh(0.5*delta))
    zz = zLength*f2

def res(delta, z1, dxi):
    f2 = 0.5*(1.0 + np.tanh((z1-0.5)*delta)/np.tanh(0.5*delta))
    return f2 - dxi;

def res_delta(delta, z1, dxi):
    td2 = np.tanh(0.5*delta)
    tzd = np.tanh((z1-0.5)*delta)
    f2_delta = 0.5*( (z1-0.5)*(1 - tzd*tzd)/td2 - 0.5*tzd*(1 - td2*td2)/(td2*td2) )
    return f2_delta;

if (dZ > 0.0):

    # invert for delta that gives dz
    z1 = 1.0/(lNode-1)
    dxi = dZ/zLength
    print "dxi = " , dxi

    delta = 1.0
    r = res(delta, z1, dxi)
    r_d = res_delta(delta, z1, dxi)

    cnt = 0
    
    while (np.abs(r) > 1e-10 and cnt<100):
        delta = delta - r/r_d
        r = res(delta, z1, dxi)
        r_d = res_delta(delta, z1, dxi)
        cnt = cnt+1

    print "Found delta = ", delta
        
    # use delta to compute spacing
    ztmp = np.linspace(0,1,lNode);
    f2 = 0.5*(1.0 + np.tanh((ztmp-0.5)*delta)/np.tanh(0.5*delta))
    zz = zLength*f2
    print "zz[1] = " , zz[1]
    

# generate mapping from natural, obvious point numbering to
# one where receiver points are at the end
renumber = np.zeros(nNode*mNode*lNode, dtype=np.int)
cnt = 0

# interior
for kNode in range(lNode):
    for jNode in range(1,mNode-1):
        # iNode = 0 and iNode = nNode-1 are receivers
        for iNode in range(1,nNode-1):
            ind = kNode*mNode*nNode + jNode*nNode + iNode # original numbering
            renumber[ind] = cnt
            cnt = cnt + 1

nPointDomain = cnt
print "nPointDomain = ", nPointDomain
            
# periodic fringe on left
for kNode in range(lNode):
    for jNode in range(mNode):
        iNode = 0
        ind = kNode*mNode*nNode + jNode*nNode + iNode # original numbering
        renumber[ind] = cnt
        cnt = cnt + 1

# periodic fringe on right
for kNode in range(lNode):
    for jNode in range(mNode):
        iNode = nNode-1
        ind = kNode*mNode*nNode + jNode*nNode + iNode # original numbering
        renumber[ind] = cnt
        cnt = cnt + 1

# periodic fringe on front
for kNode in range(lNode):
    jNode = 0
    for iNode in range(1,nNode-1):
        ind = kNode*mNode*nNode + jNode*nNode + iNode # original numbering
        renumber[ind] = cnt
        cnt = cnt + 1

# periodic fringe on right
for kNode in range(lNode):
    jNode = mNode-1
    for iNode in range(1,nNode-1):
        ind = kNode*mNode*nNode + jNode*nNode + iNode # original numbering
        renumber[ind] = cnt
        cnt = cnt + 1

nPoint = cnt
print "nPoint = ", nPoint
print 'renumber[0] = ', renumber[0]
    
Mesh_File = open(options.filename,"w")

Mesh_File.write( "%\n" )
Mesh_File.write( "% Problem dimension\n" )
Mesh_File.write( "%\n" )
Mesh_File.write( "NDIME=3\n" )
Mesh_File.write( "%\n" )
Mesh_File.write( "% Inner elements\n" )
Mesh_File.write( "%\n" )
Mesh_File.write( "NELEM=%s\n" % ((lNode-1)*(nNode-1)*(mNode-1)))

iElem = 0
for kNode in range(lNode-1):
    for jNode in range(mNode-1):
        for iNode in range(nNode-1):
            Point0 = kNode*mNode*nNode + jNode*nNode + iNode
            Point1 = kNode*mNode*nNode + jNode*nNode + iNode + 1
            Point2 = kNode*mNode*nNode + (jNode+1)*nNode + (iNode+1)
            Point3 = kNode*mNode*nNode + (jNode+1)*nNode + iNode
            Point4 = (kNode+1)*mNode*nNode + jNode*nNode + iNode
            Point5 = (kNode+1)*mNode*nNode + jNode*nNode + iNode + 1
            Point6 = (kNode+1)*mNode*nNode + (jNode + 1)*nNode + (iNode + 1)
            Point7 = (kNode+1)*mNode*nNode + (jNode + 1)*nNode + iNode
            Mesh_File.write( "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n" %
                             (KindElem,
                              renumber[Point0], renumber[Point1],
                              renumber[Point2], renumber[Point3],
                              renumber[Point4], renumber[Point5],
                              renumber[Point6], renumber[Point7], iElem) )
            iElem = iElem + 1

Mesh_File.write( "%\n" )
Mesh_File.write( "NPOIN= %s \t %s\n" % (nPoint, nPointDomain))
iPoint = 0
XYZI = np.zeros((nPoint, 3))
for kNode in range(lNode):
    for jNode in range(mNode):
        for iNode in range(nNode):
            #Mesh_File.write( "%15.14f \t %15.14f \t %15.14f \t %s\n" % (xLength*float(iNode)/float(nNode-1), yLength*float(jNode)/float(mNode-1),zz[kNode], iPoint) )
            XYZI[iPoint,0] = xLength*float(iNode)/float(nNode-1)
            XYZI[iPoint,1] = yLength*float(jNode)/float(mNode-1)
            XYZI[iPoint,2] = zz[kNode]
            
            iPoint = iPoint + 1

ind = np.argsort(renumber)
            
for iPoint in range(nPoint):
    Mesh_File.write( "%15.14f \t %15.14f \t %15.14f \t %s\n" %
                     (XYZI[ind[iPoint],0], XYZI[ind[iPoint],1], XYZI[ind[iPoint],2], iPoint) )
            
Mesh_File.write( "%\n" )
Mesh_File.write( "% Boundary elements\n" )
Mesh_File.write( "%\n" )
Mesh_File.write( "NMARK=8\n" )

Mesh_File.write( "MARKER_TAG= left\n" )
elem = (nNode-1)*(mNode-1);
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem)
for jNode in range(mNode-1):
    for iNode in range(nNode-1):
        p0 = jNode*nNode + iNode
        p1 = jNode*nNode + (iNode+1)
        p2 = (jNode + 1)*nNode + (iNode + 1)
        p3 = (jNode + 1)*nNode + iNode
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" %
                         (KindBound,
                          renumber[p0], renumber[p1], renumber[p2], renumber[p3]))

Mesh_File.write( "MARKER_TAG= right\n" )
elem = (nNode-1)*(mNode-1)
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
for jNode in range(mNode-1):
    for iNode in range(nNode-1):
        p0 = nNode*mNode*(lNode - 1) + jNode*nNode + iNode
        p1 = nNode*mNode*(lNode - 1) + jNode*nNode + iNode + 1
        p2 = nNode*mNode*(lNode - 1) + (jNode + 1)*nNode + (iNode + 1)
        p3 = nNode*mNode*(lNode - 1) + (jNode + 1)*nNode + iNode
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" %
                         (KindBound,
                          renumber[p0], renumber[p1], renumber[p2], renumber[p3]))

Mesh_File.write( "MARKER_TAG= lower\n" )
elem = (nNode-1)*(lNode-1)
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
for iNode in range(nNode-1):
    for kNode in range(lNode-1):
        p0 = iNode + kNode*nNode*mNode
        p1 = iNode + (kNode+1)*nNode*mNode
        p2 = iNode + 1 + (kNode+1)*nNode*mNode
        p3 = iNode + 1 + kNode*nNode*mNode
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" %
                         (KindBound,
                          renumber[p0], renumber[p1], renumber[p2], renumber[p3]))

Mesh_File.write( "MARKER_TAG= upper\n" )
elem = (nNode-1)*(lNode-1)
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
for iNode in range(nNode-1):
    for kNode in range(lNode-1):
        p0 = (nNode*mNode - 1) - iNode + kNode*nNode*mNode
        p1 = (nNode*mNode - 1) - iNode + (kNode+1)*nNode*mNode
        p2 = (nNode*mNode - 1) - (iNode + 1) + (kNode+1)*nNode*mNode
        p3 = (nNode*mNode - 1) - (iNode + 1) + kNode*nNode*mNode
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" %
                         (KindBound,
                          renumber[p0], renumber[p1], renumber[p2], renumber[p3]))

Mesh_File.write( "MARKER_TAG= outlet\n" )
elem = (mNode-1)*(lNode-1)
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
for jNode in range(mNode-1):
    for kNode in range(lNode-1):
        p0 = jNode*nNode + (nNode - 1) + kNode*nNode*mNode
        p1 = (jNode + 1)*nNode + (nNode - 1) + kNode*nNode*mNode
        p2 = (jNode + 1)*nNode + (nNode - 1)+ (kNode+1)*nNode*mNode
        p3 = jNode*nNode + (nNode - 1)+ (kNode+1)*nNode*mNode
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" %
                         (KindBound,
                          renumber[p0], renumber[p1], renumber[p2], renumber[p3]))
        
Mesh_File.write( "MARKER_TAG= inlet\n" )
elem = (mNode-1)*(lNode-1)
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
for jNode in range(mNode-2, -1, -1):
    for kNode in range(lNode-1):
        p0 = (jNode + 1)*nNode + kNode*nNode*mNode
        p1 = jNode*nNode + kNode*nNode*mNode
        p2 = jNode*nNode+ (kNode+1)*nNode*mNode
        p3 = (jNode + 1)*nNode+ (kNode+1)*nNode*mNode
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" %
                         (KindBound,
                          renumber[p0], renumber[p1], renumber[p2], renumber[p3]))
        
Mesh_File.write( "MARKER_TAG= SEND_RECEIVE\n" )
print lNode*2*(nNode+mNode-2)
Mesh_File.write( "MARKER_ELEMS= %d\n" % (lNode*2*(nNode+mNode-2)) )
#Mesh_File.write( "MARKER_ELEMS= %d\n" % (lNode*2*(nNode+mNode-2)-4) )
Mesh_File.write( "SEND_TO= 1\n" ) # these nodes SEND!
for kNode in range(0,lNode):
    offset = kNode*mNode*nNode
    Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+(mNode-1)*nNode-2], 1) ) # corner circle
    Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+2*nNode-2        ], 2) ) # corner square
    Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+nNode+1          ], 3) ) # corner triangle
    Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+(mNode-2)*nNode+1], 4) ) # corner star
    
    for ind in range(nNode+1, 2*nNode-1):
        Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+ind], 7) ) # line A'

    for ind in range((mNode-2)*nNode+1, (mNode-1)*nNode-1):
        Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+ind], 8) ) # line B'
        
    for ind in range(nNode+1, (mNode-2)*nNode+2, nNode):
        Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+ind], 5) ) # line C'

    for ind in range(2*nNode-2, (mNode-1)*nNode-1, nNode):
        Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+ind], 6) ) # line D'


Mesh_File.write( "MARKER_TAG= SEND_RECEIVE\n" )
Mesh_File.write( "MARKER_ELEMS=%s\n" % (lNode*2*(nNode+mNode-2)) )
#Mesh_File.write( "MARKER_ELEMS= %d\n" % (lNode*2*(nNode+mNode-2)-4) )
Mesh_File.write( "SEND_TO=-1\n" ) # these nodes RECEIVE!
for kNode in range(0,lNode):
    offset = kNode*mNode*nNode
    Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+0                ], 3) ) # corner o
    Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+(mNode-1)*nNode  ], 4) ) # corner s
    Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+nNode*mNode-1    ], 1) ) # corner t
    Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+nNode-1          ], 2) ) # corner star

    for ind in range((mNode-1)*nNode+1,nNode*mNode-1):
        Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+ind], 8) ) # line A

    for ind in range(1, nNode-1):
        Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+ind], 7) ) # line B

    for ind in range(2*nNode-1, (mNode-1)*nNode, nNode):
        Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+ind], 6) ) # line C

    for ind in range(nNode, (mNode-2)*nNode+1, nNode):
        Mesh_File.write( "%s \t %s \t %s\n" % (1, renumber[offset+ind], 5) ) # line D

Mesh_File.write( "NPERIODIC= 9\n" )

Mesh_File.write( "PERIODIC_INDEX= 0\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )

Mesh_File.write( "PERIODIC_INDEX= 1\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (-Lx, -Ly, 0.0) )

Mesh_File.write( "PERIODIC_INDEX= 2\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (-Lx, Ly, 0.0) )

Mesh_File.write( "PERIODIC_INDEX= 3\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (Lx, Ly, 0.0) )

Mesh_File.write( "PERIODIC_INDEX= 4\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (Lx, -Ly, 0.0) )

Mesh_File.write( "PERIODIC_INDEX= 5\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (Lx, 0.0, 0.0) )

Mesh_File.write( "PERIODIC_INDEX= 6\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (-Lx, 0.0, 0.0) )

Mesh_File.write( "PERIODIC_INDEX= 7\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, Ly, 0.0) )

Mesh_File.write( "PERIODIC_INDEX= 8\n" )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, 0.0, 0.0) )
Mesh_File.write( "%s \t %s \t %s\n" % (0.0, -Ly, 0.0) )


#Mesh_File.write( "FFD_NBOX=1\n")
#Mesh_File.write( "FFD_NLEVEL=1\n")
#Mesh_File.write( "FFD_TAG=0\n")
#Mesh_File.write( "FFD_LEVEL=0\n")
#Mesh_File.write( "FFD_DEGREE_I=6\n")
#Mesh_File.write( "FFD_DEGREE_J=6\n")
#Mesh_File.write( "FFD_DEGREE_K=1\n")
#Mesh_File.write( "FFD_PARENTS=0\n")
#Mesh_File.write( "FFD_CHILDREN=0\n")
#Mesh_File.write( "FFD_CORNER_POINTS=8\n")
#Mesh_File.write( "4.0	0	-0.1\n")
#Mesh_File.write( "6.0	0	-0.1\n")
#Mesh_File.write( "6.0	2.0	-0.1\n")
#Mesh_File.write( "4.0	2.0	-0.1\n")
#Mesh_File.write( "4.0	0	0.1\n")
#Mesh_File.write( "6.0	0	0.1\n")
#Mesh_File.write( "6.0	2.0	0.1\n")
#Mesh_File.write( "4.0	2.0	0.1\n")
#Mesh_File.write( "FFD_CONTROL_POINTS=0\n")
#Mesh_File.write( "FFD_SURFACE_POINTS=0\n")


    
Mesh_File.close()
