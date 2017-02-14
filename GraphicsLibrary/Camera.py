"""
This is a camera class
It creates the view of cameras
"""
from pyglet.gl import *

class Camera():
    x,y,z=0,0,512
    rx,ry,rz=0,0,0
    w,h=2*640,2*480
    far=8192
    fov=60
    keep = True
    def __init__(self,mode):
        self.mode = mode
          
    def view1(self,width,height):
        self.w,self.h=width,height
        self.w,self.h=width,height
        glViewport(0, 0, width, height)
        print "Viewport "+str(width)+"x"+str(height)
        if self.mode==2: self.isometric()
        elif self.mode==3: self.perspective()
        else: self.default()
        
    def view2(self,width,height):
        self.w,self.h=width,height
        glViewport(0, 0, width, height)
        print "Viewport "+str(width)+"x"+str(height)
        if self.mode==2: self.isometric()
        elif self.mode==3: self.perspective()
        else: self.default()
            
    def default(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(-self.w/2, self.w/2, -self.h/2, self.h/2, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        
    def isometric(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glScalef(1,1,1)
        glOrtho(-self.w/2.,self.w/2.,-self.h/2.,self.h/2.,0,self.far)
        glMatrixMode(GL_MODELVIEW)
        
    def perspective(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glScalef(1,1,1)
        gluPerspective(self.fov, float(self.w)/self.h, 0.1, self.far)
        glMatrixMode(GL_MODELVIEW)
      
    def key_down(self, symbol, modifiers):
        global clip_plane
        if symbol==key.F1:
            self.mode=1
            self.default()
            print "Projection: Pyglet default" 
        elif symbol==key.F2: 
            print "Projection: 3D Isometric" 
            self.mode=2
            self.isometric()
        elif symbol==key.F3:
            print "Projection: 3D Perspective" 
            self.mode=3
            self.perspective()
        elif self.mode==3 and symbol==key.MINUS:
            self.fov-=1
            self.perspective()
        elif self.mode==3 and symbol==key.EQUAL:
            self.fov+=1     
            self.perspective() 
        elif symbol == key.R:
            global BatchDraw
            BatchDraw = not BatchDraw
        elif symbol == key.I:
            clip_plane.addXnegative(10)
            clip_plane.addXpositive(-10)
        elif symbol == key.O:
            clip_plane.addXnegative(-10)
            clip_plane.addXpositive(10)
        elif symbol == key.J:
            clip_plane.addYnegative(10)
            clip_plane.addYpositive(-10)
        elif symbol == key.K:
            clip_plane.addYnegative(-10)
            clip_plane.addYpositive(10)
        elif symbol == key.N:
            clip_plane.addZnegative(10)
            clip_plane.addZpositive(-10)
        elif symbol == key.M:
            clip_plane.addZnegative(-10)
            clip_plane.addZpositive(10)
            
        else: print "KEY "+key.symbol_string(symbol)
        
    def motion(self,motion):
        if motion == key.MOTION_UP and self.mode==3:
            self.fov-=5    
            self.perspective() 
        if motion == key.MOTION_DOWN and self.mode==3:
            self.fov+=5     
            self.perspective() 
        if motion == key.MOTION_LEFT:
            self.rz+=10 
        if motion == key.MOTION_RIGHT:
            self.rz-=10 
        if motion == key.MOTION_UP and self.mode==2:
            self.rx+=10 
        if motion == key.MOTION_DOWN and self.mode==2:
            self.rx-=10 
        
    def drag(self, x, y, dx, dy, button, modifiers):
        if button==1:
            self.x-=dx*2
            self.y-=dy*2
        elif button==2:
            self.x-=dx*2
            self.z-=dy*2
        elif button==4:
            self.ry+=dx/4.
            self.rx-=dy/4.
        
    def apply(self):
        glLoadIdentity()
        if self.mode==1: return
        glTranslatef(-self.x,-self.y,-self.z)
        glRotatef(self.rx,1,0,0)
        glRotatef(self.ry,0,1,0)
        glRotatef(self.rz,0,0,1) 