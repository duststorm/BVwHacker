#!/usr/bin/python
# -*- coding: utf-8 -*-

## Modified version of the BVH importer borrowed from the makehuman project
## some parts borrowed from bvhplay


"""
**Project Name:**      MakeHuman

**Product Home Page:** http://www.makehuman.org/

**Code Home Page:**    http://code.google.com/p/makehuman/

**Authors:**           Marc Flerackers, Jonas Hauquier

**Copyright(c):**      MakeHuman Team 2001-2012

**Licensing:**         GPL3 (see also http://sites.google.com/site/makehumandocs/licensing)

**Coding Standards:**  See http://sites.google.com/site/makehumandocs/developers-guide

Abstract
--------

BVH importer
"""

# Libraries from makehuman
from aljabr import vadd, makeUnit, degree2rad, makeTranslation, mmul, euler2matrix

# Standard libraries
from math import radians, cos, sin, pi
import math

# Numpy
from numpy import array, dot


#######################################
# JOINT class (formerly BONE)
# A BVH "joint" is a single vertex with potentially MULTIPLE
# edges.  It's not accurate to call these "bones" because if
# you rotate the joint, you rotate ALL attached bones.

class Joint:
  def __init__(self, name):
    self.name = name
    self.children = []
    # Ordered list of channels: each
    # list entry is one of [XYZ]position, [XYZ]rotation
    self.channels = []
    self.frames = []
    self.parent = None
    self.strans = array([0.,0.,0.]) # static translation vector

    self.order = None

    # Transformation matrix:
    self.stransmat = array([ [0.,0.,0.,0.],[0.,0.,0.,0.],    \
                               [0.,0.,0.,0.],[0.,0.,0.,0.] ])


  def info(self):
    print "Joint name:", self.name
    print " %s is connected to " % self.name,
    if(len(self.children) == 0):
        print "nothing"
    else:
      for child in self.children:
        print "%s " % child.name,
      print
    for child in self.children:
      child.info()


  def updateFrame(self, frame, scale=1):
    xpos = 0.0
    ypos = 0.0
    zpos = 0.0

    xrot = 0.0
    yrot = 0.0
    zrot = 0.0

    # We have to build up drotmat one rotation value at a time so that
    # we get the matrix multiplication order correct.
    drotmat = array([ [1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.] ])

    # Build a translation matrix for this keyframe
    dtransmat = array([ [1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.] ])

    # Apply orientation for frame of motion
    if frame >= 0 and frame < len(self.frames):
        for index, channel in enumerate(self.channels):
            if channel == 'Xposition':
                xpos = scale*self.frames[frame][index]
                dtransmat[0,3] = xpos
            elif channel == 'Yposition':
                ypos = scale*self.frames[frame][index]
                dtransmat[1,3] = ypos
            elif channel == 'Zposition':
                zpos = scale*self.frames[frame][index]
                dtransmat[2,3] = zpos
              
            if channel == 'Xrotation':
                xrot = self.frames[frame][index]
                theta = radians(xrot)
                mycos = cos(theta)
                mysin = sin(theta)
                drotmat2 = array([ [1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.], [0.,0.,0.,1.] ])
                drotmat2[1,1] = mycos
                drotmat2[1,2] = -mysin
                drotmat2[2,1] = mysin
                drotmat2[2,2] = mycos
                drotmat = dot(drotmat, drotmat2)
            elif channel == 'Yrotation':
                yrot = self.frames[frame][index]
                theta = radians(yrot)
                mycos = cos(theta)
                mysin = sin(theta)
                drotmat2 = array([ [1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.], [0.,0.,0.,1.] ])
                drotmat2[0,0] = mycos
                drotmat2[0,2] = mysin
                drotmat2[2,0] = -mysin
                drotmat2[2,2] = mycos
                drotmat = dot(drotmat, drotmat2)
            elif channel == 'Zrotation':
                zrot = self.frames[frame][index]
                theta = radians(zrot)
                mycos = cos(theta)
                mysin = sin(theta)
                drotmat2 = array([ [1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.], [0.,0.,0.,1.] ])
                drotmat2[0,0] = mycos
                drotmat2[0,1] = -mysin
                drotmat2[1,0] = mysin
                drotmat2[1,1] = mycos
                drotmat = dot(drotmat, drotmat2)

    self.drotmat = drotmat
    self.dtransmat = dtransmat

    # At this point we should have computed:
    # Static (s)
    #  stransmat  (computed previously in __readJoint subroutine)
    # Dynamic (d)
    #  dtransmat (only non-zero if we're the hips)
    #  drotmat
    # We now have enough to compute joint.trtr and also to convert
    # the position of this joint (vertex) to worldspace.

    # trtr[time]  A premultiplied series of translation and rotation matrices

    if self.parent:
        self.localtoworld = dot(self.parent.trtr, self.stransmat)
    else:
        self.localtoworld = dot(self.stransmat, self.dtransmat)

    # Add rotation of this joint to stack to use for determining children positions
    # Note that position of this joint is not affected by its rotation
    self.trtr = dot(self.localtoworld,self.drotmat)

    # Position is the translation part of the mat (fourth column)
    self.worldpos = array([ self.localtoworld[0,3],
                            self.localtoworld[1,3],
                            self.localtoworld[2,3], 
                            self.localtoworld[3,3] ])

    for child in self.children:
        child.updateFrame(frame)


  def parseChannels(self, channels):
        self.order = ""     # order in which rotation angles are stored
        for channel in channels:
            if channel == "Xrotation":
                self.order += "x"
            elif channel == "Yrotation":
                self.order += "z"
            elif channel == "Zrotation":
                self.order += "y"
        return


class Skeleton:
    
    def __init__(self, filename, scale):
        
        self.file = open(filename, 'r')

        # Read hierarchy
        self.__expectKeyword('HIERARCHY')
            
        items = self.__expectKeyword('ROOT')
        
        self.root = Joint(items[1])
        
        self.__readJoint(self.root, scale)
            
        # Read motion
        self.__expectKeyword('MOTION')
        
        items = self.__expectKeyword('Frames:')
        self.frames = int(items[1])
        items = self.__expectKeyword('Frame') # Time:
        self.frameTime = float(items[2])

        for i in range(self.frames):
            line = self.file.readline()
            items = line.split()
            data = [float(item) for item in items]
            data = self.__getChannelData(self.root, data)
                
    def getJoint(self, name):
        
        return self.__getJoint(self.root, name)
        
    def __getJoint(self, joint, name):
        
        if joint.name == name:
            return joint
            
        for child in joint.children:
            j = self.__getJoint(child, name)
            if j:
                return j
                
        return None

    def __calcPosition(self, joint, scale):
        # Get OFFSET
        items = self.__expectKeyword('OFFSET')
        joint.offset = [scale*float(x) for x in items[1:]]
        joint.strans = joint.offset[:]
        joint.stransmat[0,3] = joint.strans[0]
        joint.stransmat[1,3] = joint.strans[1]
        joint.stransmat[2,3] = joint.strans[2]
        
        if joint.parent:
            joint.position = vadd(joint.parent.position, joint.offset)
        else:
            joint.position = joint.offset[:]

        # Calculate static transformation matrix
        joint.stransmat = array([ [1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.] ])
        joint.stransmat[0,3] = joint.offset[0]
        joint.stransmat[1,3] = joint.offset[1]
        joint.stransmat[2,3] = joint.offset[2]
                
    def __readJoint(self, joint, scale=0.25):
        
        self.__expectKeyword('{')

        self.__calcPosition(joint, scale)
        
        items = self.__expectKeyword('CHANNELS')
        joint.channels = items[2:]
        joint.parseChannels(items[2:])

        if int(items[1]) != len(joint.channels):
            RuntimeError('Expected %d channels found %d' % (items[1], len(joint.channels)))
        
        # Read child joints
        while 1:
            line = self.file.readline()
            items = line.split()
            
            if items[0] == 'JOINT':
                
                child = Joint(items[1])
                joint.children.append(child)
                child.parent = joint
                self.__readJoint(child, scale)
                
            elif items[0] == 'End': # Site
                
                child = Joint('End effector')
                joint.children.append(child)
                child.channels = []
                child.parent = joint
                
                self.__expectKeyword('{')

                self.__calcPosition(child, scale)
                
                self.__expectKeyword('}')
                
            elif items[0] == '}':
                
                break
                
            else:
                
                raise RuntimeError('Expected %s found %s' % ('JOINT, End Site or }', items[0]))
                    
    def __expectKeyword(self, keyword):
        
        line = self.file.readline()
        items = line.split()
        
        if items[0] != keyword:
            raise RuntimeError('Expected %s found %s' % (keyword, items[0]))
                
        return items


    def __getChannelData(self, joint, data):
        
        channels = len(joint.channels)
        joint.frames.append(data[0:channels])
        data = data[channels:]
        
        for child in joint.children:
            data = self.__getChannelData(child, data)
        
        return data
        
    def updateFrame(self, frame, scale = 1):
        self.root.updateFrame(frame, scale)
