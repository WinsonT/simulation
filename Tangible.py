from sys import path
#path.append(path[0]+"/Maps")  
#path.append(path[0]+"/NPC")  

import pygame
import os
from pygame.locals import *
from sys import exit
from time import time
from math import pi, sin, cos

from random import random, randint, randrange

pygame.init()

#TO DO
#make a class that has the following property. YES: Tangible. NO: Image, Brain, Animation
#called class Tangible(), divided into two: Polygon and Circle (no image). a subclass of Polygon and Circle (no image) will be Polygon and Circle (with image)
#
#A subclass of Polygon will be one that not only has Tangible, but also CollisionResponse, e.g., Impenetrable.
#
#An equivalent class to Polygon, but with different shapes, i.e., class Circle()

#DIAGRAM: A --> B means B import A
#                                           
# World <-------------------------->         Tangible
#   &                                        |      |
# World                                      v      v
#Template                               Circle      Polygon  <---------------------------- ProjAxis
#                                            |      |                                     in Polygon.__init__(...) to store .axes
#                                            v      v
#                                    imgCircle      imgPolygon
#
#          PASSED TO TANGIBLE                             
#           Polygon.__init__(world, ...) [right now, .world is only used to add/remove oneself from world.entity_list and world.tangible_list]
#           .render(screen, alpha)
#           .process(event)
#          ENTITIES HAVE
#           .pos
#          TANGIBLES HAVE
#           .midAABB, .width, .height (for broad phase collision detection)
#           "shape type": Polygon or Circle
#          POLYGONS HAVE
#           .vert, .normal, .axes (for narrow phase collision detection)
#

#README
#A subclass that inherits Polygon can be defined for objects that need the quality that is classified as NO for Polygon. Such
#classes automatically inherits all the quality classified as YES in Polygon.

#DESCRIPTION
#A Polygon, with many sides. YES: Tangible, Image. NO: Brain, Animation, CollisionResponse
#Tangible does not mean Impenetrable. It is so if the collision response is rigid body non-penetration constraint,
#but Tangible here only means interacting with others through touch. Collision response can be sth. else, e.g., burn the object upon 
#contact. In fact, Polygon does not have a collision response. Hence, Polygon is only a template to be used by other classes
#that own a particular collision response. Also, Polygon does not have an animation. Animation is specific to object type, hence objects that do have animation
#should be defined in its own class, e.g., class Box has an animation for when the Box is broken.

#RULES on defining vertices0:
#1. vertices0 are defined CLOCKWISE in direction
#3. object must be convex (all corners have angles LESS than 180 degrees)
#4. preferrably oriented such that at theta=0 AABB is the tightest == OBB

#PARAMETERS
#world: The World (a class) in which the Polygon is
#pos [x,y,theta]: The position of the Polygon, theta is the degree of CLOCKWISE rotation.
#vertices0 [[x1,y1],[x2,y2],...]: Coordinates of the vertices at theta == 0. Origin is also the Center of Mass (CM)
#image: file name of the image of the Polygon
#imgCM [x,y]: the coordinate of the center of mass (in pixel) in the image, using the top left as the origin.
#speed [vx,vy,w]: time derivatives of pos
class Polygon():
    def __init__ (self, world, pos, vertices0, image, imgCM, speed = [0,0,0]):
        #to implement in C++, use typename template since Entity needs World and World needs Entity
        self.world = world
        
        self.pos = pos
        self.speed = speed
        
#        self.prevPos = copy_list2(self.pos)
        
        self.massInv = 1.0
        self.inertiaInv = self.massInv/1000 #0.001

        self.vert0 = vertices0
        self.vert = copy_list2(self.vert0)
        
        '''number of vertices, >= 3'''
        self.n = 0
        for i in self.vert0:
            self.n += 1

        '''find top left and bottom right coordinates of AABB using center of mass as origin'''
        x = [self.vert0[0][0],self.vert0[0][0]]
        y = [self.vert0[0][1],self.vert0[0][1]]
        for i in self.vert0: #the first run is actually a waste
            if i[0] < x[0]:
                x[0] = i[0]
            elif i[0] > x[1]:
                x[1] = i[0]
            if i[1] < y[0]:
                y[0] = i[1]
            elif i[1] > y[1]:
                y[1] = i[1]
                    
        self.width0 = x[1] - x[0]
        self.height0 = y[1] - y[0]
        
        #The middle of the AABB
        self.midAABB0 = ((x[0]+x[1])/2.0, (y[0]+y[1])/2.0)
        self.update_wh() #this defines .width, .height, and .midAABB

    '''create a surface such that the middle of the image is the center of mass'''
        #imgCM is the center of mass position relative to the top left of image. Now we make CM to also be the center of image
        imgWidth0 = imgCM[0]*2-1
        imgHeight0 = imgCM[1]*2-1
        self.image0 = pygame.Surface([imgWidth0, imgHeight0, pygame.SRCALPHA, 32)
        self.image0 = self.image0.convert_alpha()
        image0 = pygame.image.load(image).convert_alpha()
        self.image0.blit(image0, (0,0))

        self.updateImg(self.pos[2])  #this defines .imgWidth, .imgHeight, and .lastImgTheta
        
    
    
#        self.active = 1
#        self.active_state = [(self.pos[0], self.pos[1], self.pos[2]), 0]


#        self.forces = []

        
'''defining normal vectors in local coordinates, including max min projections of vertices onto each normal'''
        self.axes = []
        self.nn = 0 #the number of normals

        skippedPoints = []
        for i in range(self.n):
            if i in skippedPoints: #in cpp skippedPoints can be an array of length n, intialized as zeros. if skippedPoints[i] == 1, continue
                continue #this statement can be used in cpp
                        
            j = self.nextVertex(i)

            a = self.vert0[i][0] - self.vert0[j][0]
            b = self.vert0[i][1] - self.vert0[j][1]
            
            c = (a*a + b*b)**0.5
            
            #Now vertices are defined CLOCKWISE, so that [b/c, -a/c] is the unit vector pointing OUTWARDS. Hence, vert0[i] and [j] are MAX projections
            maxValue = (self.vert0[i][1]*self.vert0[j][0]-self.vert0[i][0]*self.vert0[j][1])/c
            
            #to find the MIN, we traverse the vertices until a turning point is found. We shall traverse CLOCKWISE.
            k = self.nextVertex(j)
            #find the projection of k
            minValue = self.vert0[k][0]*b/c - self.vert0[k][1]*a/c
            
            #special parameter. If loop exits normally, check == 0 and this guarantees that there is only one MIN vertex
            check = 0
            
            k_final = self.prevVertex(i)
            while (k != k_final):
                k = self.nextVertex(k)
            
                #The projection of vert0[k]
                value = self.vert0[k][0]*b/c - self.vert0[k][1]*a/c
                if value < minValue:
                    minValue = value
                else:
                    check = 1
                    break
            
            #If there are two MIN vertices. SHORT CIRCUIT: value might be undefined, but in this case check is FALSE.
            if check and int(1000*(value - minValue)) == 0:
                skippedPoints.append(prevVertex(k))
                self.axes.append( ProjAxis((b/c,-a/c), maxValue, i, j, minValue, self.prevVertex(k), k) )
            else:
                self.axes.append( ProjAxis((b/c,-a/c), maxValue, i, j, minValue, k) )

            self.nn += 1

        '''normal vectors in world coordinates'''
        self.normal = []
        for i in self.axes:
            self.normal.append(i.normal0)

        '''for .vert and .normal'''
        self.updateRotation()

    def nextVertex(self, i):
        #0 <= i <= self.n-1
        if i == self.n-1:
            return 0
        else:
            return i+1

    def prevVertex(self, i):
        #0 <= i <= self.n-1
        if i == 0:
            return self.n-1
        else:
            return i-1

    #AABB size update every time .pos[2] is updated
    def update_wh(self):
        theta = self.pos[2]
        c = cos(theta)
        s = sin(theta)
        self.midAABB = [c*self.midAABB0[0] - s*self.midAABB0[1], c*self.midAABB0[1] + s*self.midAABB0[0]]
        
        if theta > pi: #rotation by pi doesnt change width/height, so theta and theta-pi are equivalent
            theta -= pi
        if theta > pi/2.0: #finally rotation by pi/2 + alpha is equivalent to pi/2-alpha. effectively, in the end 0 < theta < pi/2
            theta = pi - theta
    
        c = cos(theta)
        s = sin(theta)
        self.width = c*self.width0 + s*self.height0
        self.height = c*self.height0 + s*self.width0
    
    #image update prior rendering
    def updateImg(self, theta):
        self.image = pygame.transform.rotate(self.image0, -180*theta/pi)
        self.imgWidth, self.imgHeight = self.image.get_size()
        self.lastImgTheta = theta
                                      
    #vertices and normals update, done only prior to separating axis test IF .speed[2] != 0
    def updateRotation(self):
        c = cos(self.pos[2])
        s = sin(self.pos[2])
        for i in range(self.n):
            vert0_ix = self.vert0[i][0]
            vert0_iy = self.vert0[i][1]
            self.vert[i] = ( vert0_ix*c - vert0_iy*s , \
                             vert0_ix*s + vert0_iy*c )
        for i in range(self.nn):
            normal0_ix = (self.axes[i]).normal0[0]
            normal0_iy = (self.axes[i]).normal0[1]
            self.normal[i] = ( normal0_ix*c - normal0_iy*s , \
                               normal0_ix*s + normal0_iy*c )

    def render(self, surface, alpha = 0):
        x = self.pos[0] - (self.speed[0] * self.world.dt)*(1-alpha)
        y = self.pos[1] - (self.speed[1] * self.world.dt)*(1-alpha)
        theta = self.pos[2] - (self.speed[2] * self.world.dt)*(1-alpha)
        
        if theta != self.lastImgTheta:
            theta = resetTheta(theta)
            self.updateImg()
        
        #there is a rounding off of the coordinates to integer
        surface.blit(self.image,(x - width/2,y - height/2))
            
   
    def trans(self): #With euler's method, prevPos is not needed as it can simply be calculated using the velocity
        self.pos[0] = self.pos[0]+self.speed[0]*self.world.dt
        self.pos[1] = self.pos[1]+self.speed[1]*self.world.dt
        
        #note: image rotation is done prior to rendering, to avoid unnecessarily updating image every physics time step
        if self.speed[2] != 0:
            self.pos[2] = resetTheta(self.pos[2]+self.speed[2]*self.world.dt)
            self.update_wh()

'''For objects that can be destroyed,
create a new class where in __init__ define
self.health = 1
    def damaged(self, dmg = 0):
        self.health -= dmg
        if self.health <= 0:
            self.death()

    def death(self):
        self.world.remove_entity(self.id_num)

Alternative death(). Object can be divided into pieces (e.g., a box breaking apart) where the pieces can not collide
with anything and will simply collapse into the ground, where they will be deleted as soon
as they are off screen.
Note: If the box is hit using a sword or something, the pieces will only have a very small
initial velocity and will fall almost vertically into the ground. If the box is blown up
using a magic/bomb, then the pieces might fly around'''

                                      
    def process(self, event):
        pass
##        self.contacts = {}
##        self.forces = []
#        
#        
#        #activation and deactivation still need some work out to prevent objects from floating in mid air
#
#        if abs(self.pos[0] - self.active_state[0][0]) >= 1 or\
#           abs(self.pos[1] - self.active_state[0][1]) >= 1 or\
#           abs(self.pos[2] - self.active_state[0][2]) >= 0.1:
#            if self.id_num == 0:
#                print 'negated'
#            self.active = 1
#            self.active_state = [(self.pos[0], self.pos[1], self.pos[2]), 0]
#                
#        elif self.active:
#            
#            self.active_state[1] += self.world.dt
#                
#            if self.active_state[1] > deactive_time:
#                if self.id_num == 0:
#                    print '!!!!!!!!!!!!!!'
#                self.active = 0
#                self.speed = [0,0,0]
#                self.active_state = [(self.pos[0], self.pos[1], self.pos[2]), 0]
#
#
#        if self.active: pass
#
##            else: self.trans()
#        self.register_tile() #omit this step if the object does not collide. Thus far only class Polygon can collide                


    #the RK4 integrator
    def trans_new(self):
        #the tuples a, b, c, and d have two components. the first has a unit of velocity, the second acceleration. 
        a = self.eval_derivatives(self.world.t, 0, ((0,0,0),(0,0,0)))
        b = self.eval_derivatives(self.world.t + 0.5*self.world.dt, 0.5*self.world.dt, a)
        c = self.eval_derivatives(self.world.t + 0.5*self.world.dt, 0.5*self.world.dt, b)
        d = self.eval_derivatives(self.world.t + self.world.dt, self.world.dt, c)
        
        dxdt = ((a[0][0] + 2*(b[0][0] + c[0][0]) + d[0][0])/6.0,\
                (a[0][1] + 2*(b[0][1] + c[0][1]) + d[0][1])/6.0,\
                (a[0][2] + 2*(b[0][2] + c[0][2]) + d[0][2])/6.0)
        dvdt = ((a[1][0] + 2*(b[1][0] + c[1][0]) + d[1][0])/6.0,\
                (a[1][1] + 2*(b[1][1] + c[1][1]) + d[1][1])/6.0,\
                (a[1][2] + 2*(b[1][2] + c[1][2]) + d[1][2])/6.0)

        self.pos[0] += dxdt[0] * self.world.dt
        self.pos[1] += dxdt[1] * self.world.dt
        self.pos[2] += dxdt[2] * self.world.dt
        
        self.speed[0] += dvdt[0] * self.world.dt
        self.speed[1] += dvdt[1] * self.world.dt
        self.speed[2] += dvdt[2] * self.world.dt

        if self.pos[2] >= 2*pi:
            self.pos[2] -= 2*pi
            self.prevPos = (self.prevPos[0], self.prevPos[1], self.prevPos[2]-2*pi) #to prevent stutter due to anti-(temporal)aliasing
        elif self.pos[2] <= -2*pi:
            self.pos[2] += 2*pi
            self.prevPos = (self.prevPos[0], self.prevPos[1], self.prevPos[2]+2*pi) #to prevent stutter due to anti-(temporal)aliasing

        if self.prevPos[2] != self.pos[2]: #note: image rotation is done prior to blitting, so as to not unnecessarily update image every physics time step
            self.update_wh()
            self.normal_poly_update_needed = 1 #so that normal_poly is updated when checking collision


    def eval_derivatives(self, t, dt, derivatives):
        x = (self.pos[0] + derivatives[0][0]*dt, self.pos[1] + derivatives[0][1]*dt, self.pos[2] + derivatives[0][2]*dt)
        v = (self.speed[0] + derivatives[1][0]*dt, self.speed[1] + derivatives[1][1]*dt, self.speed[2] + derivatives[1][2]*dt)
        return (v, self.accel(x, v, t+dt))
        
             
    def accel(self, x, v, t):
        force = [0,0]
        torque = 0
        
        #think about how to implement friction, because it oftens need to know what other forces are acting on the body
        for i in self.forces:
            ft = i.calculate(x, v, t)
            force[0] += ft[0][0]
            force[1] += ft[0][1]
            torque += ft[1]
            
  


#keeping 0 <= theta <= 2*pi
def resetTheta(theta):
    if theta > 2*pi:
        theta -= 2*pi
    elif theta < 0:
        theta += 2*pi
    return theta

        
#returns the coordinates of the vertices of an n-sided regular polygon of side length == sideLength,
#with the center as the origin.
def regularPolygon(n, sideLength): 
    theta = 2*pi/n
    c = cos(theta)
    s = sin(theta)
    
    #r is the radius of the common circle on which all vertices lie
    r = sideLength/(2*sin(theta/2.0))
    
    #the vertices are listed CLOCKWISE, as required by the class Polygon
    #let the first vertex be
    vertices0 = [(sideLength/2.0, r* cos(theta/2.0))] 

    #we add the next coordinate, which is the previous rotated by theta CLOCKWISE. Remembering y-axis points downwards,
    for i in range(n-1):
        x = vertices0[i][0]*c - vertices0[i][1]*s,
        y = vertices0[i][1]*c + vertices0[i][0]*s
        vertices0.append((x,y))

    return vertices0

'''
to tuple vertices0, use tuple(). i dont because for indexing, which is what vertices and normal in Polygon are only for,
there's no noticeable difference in performance between tuple and list (list may even be faster).
tuples are more concise in terms of memory and faster for instantiation like a = (1,2,3), compared to a = [1,2,3].
'''
