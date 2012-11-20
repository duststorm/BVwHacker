#!/usr/bin/python
# -*- coding: utf-8 -*-

import bvh

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys
import math


# Some api in the chain is translating the keystrokes to this octal string
# so instead of saying: ESCAPE = 27, we use the following.
ESCAPE = '\033'

# Number of the glut window.
window = 0


def InitGL(Width, Height):
    glClearColor(0.5, 0.5, 0.5, 0.0)
    glClearDepth(1.0)
    glDepthFunc(GL_LESS)
    glEnable(GL_DEPTH_TEST)
    glShadeModel(GL_SMOOTH)
    
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()

    gluPerspective(45.0, float(Width)/float(Height), 0.1, 100.0)

    glMatrixMode(GL_MODELVIEW)


def ReSizeGLScene(Width, Height):
    if Height == 0:                        # Prevent A Divide By Zero If The Window Is Too Small 
        Height = 1

    glViewport(0, 0, Width, Height)        # Reset The Current Viewport And Perspective Transformation
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45.0, float(Width)/float(Height), 0.1, 100.0)
    glMatrixMode(GL_MODELVIEW)


def DrawGLScene():
    global skeleton

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()

    # Place camera backward
    glTranslatef(0.0, 1.0, -15.0)

    drawFloorPlane(-20, 20, 10, -5, True)
    drawBVHRig(skeleton)

    glFlush()


def drawFloorPlane(pmin, pmax, lines = 10, y = 0, faces = False):
    OFFSET = 0.01
    if faces:
        glColor3f(.3,.3,.3);
        glBegin(GL_QUADS);
        glVertex3f(pmin, y-OFFSET, pmin);
        glVertex3f(pmin, y-OFFSET, pmax);
        glVertex3f(pmax, y-OFFSET, pmax);
        glVertex3f(pmax, y-OFFSET, pmin);
        glEnd();

    size = pmax-pmin

    glLineWidth(2)
    glBegin(GL_LINES);
    for i in range(lines):
        if i == 0:
            glColor3f(.6,.3,.3)
        else:
            glColor3f(.25,.25,.25)
        pos = pmin + i*(size/lines)
        glVertex3f(pos,y,pmin)
        glVertex3f(pos,y,pmax)
        if i == 0:
            glColor3f(.3,.3,.6)
        else:
            glColor3f(.25,.25,.25)
        glVertex3f(pmin,y,pos)
        glVertex3f(pmax,y,pos)
    glEnd()


def keyPressed(*args):
    global window
    global skeleton
    global frame

    # > (on numpad)
    if ord(args[0]) == 54:
        frame = frame + 1
        print "frame %s" % frame
        skeleton.updateFrame(frame)
        return
        
    # < (on numpad)
    if ord(args[0]) == 52:
        frame = frame - 1
        print "frame %s" % frame
        skeleton.updateFrame(frame)
        return
    
    # If escape is pressed, kill everything.
    if args[0] == ESCAPE:
        sys.exit()



def drawBVHRig(skeleton):
    glColor3f(0.7, 0.3, 0.2);
    drawJoint(skeleton.root)


def drawJoint(joint):
    pos = getPosition(joint)
    #if pos[0] == pos[1] == pos[2] == 0:
    #    print joint.name
    RADIUS = 0.15

    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    glColor3f(0.7, 0.3, 0.2);
    drawSphere( RADIUS, 8, 8 )
    #glutSolidSphere( RADIUS, 8, 8 )
    glColor3f(0.0, 0.0, 0.0);
    #glutWireSphere( RADIUS, 8, 8 )
    drawSphere( RADIUS, 8, 8, True)

    glColor3f(0.7, 0.3, 0.2);
    glPopMatrix()

    if joint.parent:
        head = getPosition(joint.parent)
        tail = pos
        glLineWidth(3)
        glBegin(GL_LINES)
        glVertex3f(head[0], head[1], head[2])
        glVertex3f(tail[0], tail[1], tail[2])
        glEnd();

    for child in joint.children:
        drawJoint(child)


def getPosition(joint):
    return [joint.worldpos[0], joint.worldpos[1], joint.worldpos[2]]


def drawSphere(r, lats, longs, wireFrame = False):
    i = 0
    j = 0
    for i in range(lats+1):
        lat0 = math.pi * ((-0.5 + float(i) - 1) / lats)
        z0  = math.sin(lat0)
        zr0 =  math.cos(lat0)

        lat1 = math.pi * (-0.5 + float(i) / lats)
        z1 = math.sin(lat1)
        zr1 = math.cos(lat1)

        if wireFrame:
            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE)
            glLineWidth(1)
            r = r + 0.0001

        glBegin(GL_QUAD_STRIP)
        for j in range(longs+1):
            lng = 2 * math.pi * float(j - 1) / longs
            x = math.cos(lng)
            y = math.sin(lng)

            glNormal3f(x * zr0 * r, y * zr0 * r, z0 * r)
            glVertex3f(x * zr0 * r, y * zr0 * r, z0 * r)
            glNormal3f(x * zr1 * r, y * zr1 * r, z1 * r)
            glVertex3f(x * zr1 * r, y * zr1 * r, z1 * r)
        glEnd()

        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL)


if __name__ == '__main__':
    frame = -1
    skeleton = bvh.Skeleton("cmu_mb_01_01.bvh", 0.25)
    skeleton.updateFrame(frame)


    print "Hit ESC key to quit."

    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
    glutInitWindowSize(640, 480)
    glutInitWindowPosition(0, 0)
    window = glutCreateWindow("BVwHacker")
    glutDisplayFunc(DrawGLScene)
    glutIdleFunc(DrawGLScene)
    glutReshapeFunc(ReSizeGLScene)
    glutKeyboardFunc(keyPressed)
    InitGL(640, 480)
    glutMainLoop()

