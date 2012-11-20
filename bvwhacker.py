#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
@author: Jonas Hauquier
@license:  GPL3
'''

#Uncomment if you have multiple wxWidgets versions
#import wxversion
#wxversion.select('2.8')

import math, wx
from wx import xrc
from wx.glcanvas import GLCanvas
from OpenGL.GLU import *
from OpenGL.GL import *

import bvwhacker_base 

class WxGLTest(GLCanvas):
    def __init__(self, parent):
        GLCanvas.__init__(self, parent,-1, 
                          attribList=[wx.glcanvas.WX_GL_DOUBLEBUFFER])
        wx.EVT_PAINT(self, self.OnDraw)
        wx.EVT_SIZE(self, self.OnSize)
        wx.EVT_MOTION(self, self.OnMouseMotion)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_RIGHT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_RIGHT_UP, self.OnMouseUp)
        wx.EVT_WINDOW_DESTROY(self, self.OnDestroy)
        
        self.init = True

        # initial mouse position
        self.lastx = self.x = 30
        self.lasty = self.y = 30
        self.mouseDown = False
        self.rMouseDown = False

        self.size = None


    def OnDraw(self,event):
        self.SetCurrent()

        if self.init:
            self.InitGL()
            self.init = False

            glLoadIdentity()
            # Set initial camera position
            glTranslatef(0.0, 1.0, -15.0)


        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # Mouse camera rotation
        if self.mouseDown:
            if self.size is None:
                self.size = self.GetClientSize()
            w, h = self.size
            w = max(w, 1.0)
            h = max(h, 1.0)
            xScale = 180.0 / w
            yScale = 180.0 / h
            rotX = (self.x - self.lastx) * xScale
            rotY = (self.y - self.lasty) * yScale
            #print "Rotated camera (%s, %s)" % (rotX, rotY)
            glRotatef(rotY, 1.0, 0.0, 0.0)
            glRotatef(rotX, 0.0, 1.0, 0.0)
        elif self.rMouseDown:
            if self.size is None:
                self.size = self.GetClientSize()
            _, h = self.size
            h = max(h, 1.0)
            yScale = 10.0 / h
            zoomY = (self.y - self.lasty) * yScale
            glTranslatef(0.0, 0.0, zoomY)

        # Draw skeleton
        bvwhacker_base.drawFloorPlane(-20, 20, 10, -5, True)
        bvwhacker_base.drawBVHRig(bvwhacker_base.skeleton)

        glFlush()
        self.SwapBuffers()
        #print "drawn"

        return
        


    def InitGL(self):
        '''
        Initialize GL
        '''

        glClearColor(0.5, 0.5, 0.5, 0.0)
        glClearDepth(1.0)
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        glShadeModel(GL_SMOOTH)
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()

        Width = 640
        Height = 480
        gluPerspective(45.0, float(Width)/float(Height), 0.1, 100.0)

        glMatrixMode(GL_MODELVIEW)


    def OnSize(self, event):

        try:
            width, height = event.GetSize()
        except:
            width = event.GetSize().width
            height = event.GetSize().height
        
        bvwhacker_base.ReSizeGLScene(width, height)

        self.Refresh()
        self.Update()

    def OnMouseDown(self, evt):
        #print "Mouse down"
        if evt.LeftIsDown():
            self.mouseDown = True
            self.rMouseDown = False
        elif evt.RightIsDown():
            self.rMouseDown = True
            self.mouseDown = False
        self.CaptureMouse()
        self.x, self.y = self.lastx, self.lasty = evt.GetPosition()

    def OnMouseUp(self, evt):
        #print "Mouse up"
        self.mouseDown = False
        self.rMouseDown = False
        self.ReleaseMouse()

    def OnMouseMotion(self, evt):
        if evt.Dragging() and (evt.LeftIsDown() or evt.RightIsDown()):
            self.lastx, self.lasty = self.x, self.y
            self.x, self.y = evt.GetPosition()
            self.Refresh(False)
        #print "Mouse moved (%s, %s)" % (self.x, self.y)

    def OnDestroy(self, event):
        print "Destroying Window"


def sliderUpdate(event):
    fr = slider.GetValue()
    bvwhacker_base.frame = fr
    bvwhacker_base.skeleton.updateFrame(fr, 0.1)
    status.SetLabel("frame: %s/%s" % ( bvwhacker_base.frame, 
                         len(bvwhacker_base.skeleton.root.frames) ))
    canvas.Refresh()
    canvas.Update()


if __name__ == '__main__':

    app = wx.App()

    res = xrc.XmlResource('gui/bvwhacker_gui.xrc')
    frame = res.LoadFrame(None, 'BVwHacker')
    #panel = res.LoadPanel(frame, "BVhWhackerPanel")

    #canvasPanel = xrc.XRCCTRL(frame, 'oGLCanvasPanel')
    canvas = WxGLTest(frame)
    canvas.SetSize((640,480))

    #frame.SetSize((640,480))

    bvwhacker_base.frame = -1
    bvwhacker_base.skeleton = bvwhacker_base.bvh.Skeleton("cmu_mb_01_01.bvh", 0.25)
    bvwhacker_base.skeleton.updateFrame(frame)

    slider = xrc.XRCCTRL(frame, 'frameSlider')
    print slider
    print xrc.XRCCTRL(frame, 'beginBtn')
    print xrc.XRCCTRL(frame, 'frameLbl')
    slider.SetMin(-1)
    slider.SetMax(len(bvwhacker_base.skeleton.root.frames))
    frame.Bind(wx.EVT_SLIDER, sliderUpdate)

    status = xrc.XRCCTRL(frame, 'frameLbl')
    status.SetLabel("frame: %s/%s" % ( bvwhacker_base.frame, 
                         len(bvwhacker_base.skeleton.root.frames) ))

    frame.Show()
    app.MainLoop()
