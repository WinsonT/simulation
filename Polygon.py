from sys import path
#path.append(path[0]+"/Maps")  
#path.append(path[0]+"/NPC")  

import pygame
from pygame.locals import *
from math import pi, sin, cos
from functionlib import copy_list2

import Projection
import BC

pygame.init()

#CONSTANTS
GRAD_TOL = 0.001
THETA_TOL = 0.001 #this is a bit large. around 1 degrees!
OMEGA_TOL = 0.001

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
#Template                               Circle      Polygon  <---------------------------- Projection: in Polygon.__init__(...) to store .projections
#                                            |      |                                      BC: for broad collDetect. has .c, .r
#                                            v      v                                      
#                                    imgCircle      imgPolygon                            
#
#          PASSED TO TANGIBLE                             
#           Polygon.__init__(world, ...) [.world is only used to add/remove oneself (and one's BC) to/from world.entityList (world.AABB_list)]
#           .render(screen, alpha)
#           .process(event)
#           .trans()
#          ENTITIES HAVE
#           .idNum
#           .pos
#           .render()
#           .process()
#          TANGIBLES HAVE
#           .AABB
#           "shape type": Polygon or Circle
#          POLYGONS HAVE
#           .vert, .normal, .projections (for narrow phase collision detection)
#

#README
#A subclass that inherits Polygon can be defined for objects that need the quality that is classified as NO for Polygon. Such
#classes automatically inherits all the quality classified as YES in Polygon.

#DESCRIPTION
#A Polygon, with many sides. complex collision object. YES: Tangible, Image. NO: Brain, Animation, CollisionResponse
#Tangible does not mean Impenetrable. It is so if the collision response is rigid body non-penetration constraint,
#but Tangible here only means interacting with others through touch. Collision response can be sth. else, e.g., burn the object upon 
#contact. In fact, Polygon does not have a collision response. Hence, Polygon is only a template to be used by other classes
#that own a particular collision response. Also, Polygon does not have an animation. Animation is specific to object type, hence objects that do have animation
#should be defined in its own class, e.g., class Box has an animation for when the Box is broken.

#RULES on defining vertices0:
#1. vertices0 are defined CLOCKWISE in direction
#3. object must be convex (all corners have angles LESS than 180 degrees)
#4. preferrably oriented such that at theta=0 AABB/BC is the tightest
#5. no more than 32 points!

#PARAMETERS
#world: The World (a class) in which the Polygon is
#pos [x,y,theta]: The position of the Polygon, theta is the degree of CLOCKWISE rotation; 0 <= theta <= 2*pi
#vertices0 [[x1,y1],[x2,y2],...]: Coordinates of the vertices at theta == 0. Origin is also the Center of Mass (CM)
#image: file name of the image of the Polygon
#imgCM [x,y]: the coordinate of the center of mass (in pixel) in the image, using the top left as the origin.
#speed [vx,vy,w]: time derivatives of pos
class Polygon():
    def __init__ (self, world, pos, vertices0, image, imgCM, speed = [0,0,0]):
        #to implement in C++, use typename template since Entity needs World and World needs Entity
        self.world = world
        self.idNum = None #this will be determined by the world entityList
        
        #KINEMATICS
        self.pos = pos
        self.speed = speed
        
        #DYNAMICS
        self.massInv = 1.0
        self.inertiaInv = self.massInv/1000.0

    	'''GEOMETRY'''
        self.vert0 = vertices0
        self.vert = copy_list2(self.vert0)
        
        #number of vertices, n >= 3
        self.n = 0
        for i in self.vert0:
            self.n += 1

        #maps the index of vert (key) to the index of the corresponding normal (value), the one perpendicular to the line joining vert and nextVertex(vert)
        self.vertToNormal = [0 for i in xrange(self.n)]


        '''Bounding Circle'''
        #midBC0 is the sphere center in local coordinate
        self.midBC0, radius = computeBC(self.vert0)
        self.BC = BC.BC(self, radius)
        self.updateBC() #defines the sphere center in the world coordinate
        
        self.world.addBC(self.BC) #world.BC_list should be sorted (decreasing in r). insert at end, then isort. deletion is log n using binary search


        '''create a surface such that the middle of the surface is the center of mass'''
        #imgCM is the center of mass position relative to the top left of image. Now we make CM to also be the center of image
        imgWidth0 = imgCM[0]*2-1
        imgHeight0 = imgCM[1]*2-1
        self.image0 = pygame.Surface([imgWidth0, imgHeight0], pygame.SRCALPHA, 32)
        self.image0 = self.image0.convert_alpha()
        image0 = pygame.image.load(image).convert_alpha()
        self.image0.blit(image0, (0,0))

        self.updateImg(self.pos[2])  #this defines .image, .imgWidth, .imgHeight, and .lastImgTheta
        
    
    
#        self.active = 1
#        self.active_state = [(self.pos[0], self.pos[1], self.pos[2]), 0]
#        self.forces = []

        
        '''defining normal vectors in local coordinates, including max min projections of vertices onto each normal'''
        self.projections = []
        self.normal0 = []
        self.nn = 0 #the number of normals

        #when the binary digit of pointMask is 1, the point is to be skipped
        pointMask = 0
        for i in range(self.n):
            #if the to-be-added normal is a duplicate, skip
            if 1 & pointMask:
                continue #this statement can be used in cpp
                        
            j = self.nextVertex(i)

            a = self.vert0[j][0] - self.vert0[i][0]
            b = self.vert0[j][1] - self.vert0[i][1]
            
            c = (a*a + b*b)**0.5
            
            #add normal0 to the list. 
            self.normal0.append((b/c,-a/c))
            
            #mark vertToNormal
            self.vertToNormal[i] = self.nn
            
            #Now vertices are defined CLOCKWISE, so that [b/c, -a/c] is the unit vector pointing OUTWARDS. Hence, vert0[i] and [j] are MAX projections
            maxValue = -(self.vert0[i][1]*self.vert0[j][0]-self.vert0[i][0]*self.vert0[j][1])/c
            
            #to find the MIN, we traverse the vertices until a turning point is found. We shall traverse CLOCKWISE.
            #find the projection of the vertex after j
            k = self.nextVertex(j)
            minValue = self.vert0[k][0]*b/c - self.vert0[k][1]*a/c
            
            #special parameter. initially == 1<<1, points to j (so it's lagging behind k by 1). indicates which vertex can be skipped if there is one.
            check = 2
            
            #identifying if code does not go into the while loop (n = 3).
            value = minValue + 1000*GRAD_TOL #sentinel
            
            k_final = self.prevVertex(i)
            while (k != k_final):
                k = self.nextVertex(k)
                check <<= 1
            
                #The projection of vert0[k]
                value = self.vert0[k][0]*b/c - self.vert0[k][1]*a/c
                                
                if value > minValue - GRAD_TOL:
                    #minimum is the vertex before
                    k = self.prevVertex(k)
                    #since check is lagging behind k by 1, no need to do check >>= 1. now check points to k
                    break
                else:
                    minValue = value
                    value = minValue + 1000*GRAD_TOL #sentinel. indicates there are only 1 MIN vertex
                
            #If there are two MIN vertices.
            #the way the MAXes and MINs are stored in Projection, MAX2 is the next vertex if one travels CLOCKWISE from MAX1. likewise for MIN2 and MIN1.
            if value < minValue + GRAD_TOL:
                #prevents duplicate
                pointMask |= check
                
                #mark vertToNormal
                self.vertToNormal[k] = self.nn
                
                self.projections.append( Projection.Projection(maxValue, i, minValue, k, j, self.nextVertex(k)) )
                
            else:
                self.projections.append( Projection.Projection(maxValue, i, minValue, k, j, None) )

            self.nn += 1
            pointMask >>= 1
                
        #make it constant
        self.normal0 = tuple(self.normal0)
        
        '''normal vectors in world coordinates'''
        self.normal = []
        for normal in self.normal0:
            self.normal.append(normal)

        '''for .vert and .normal'''
        self.updateRotation() #this defines .lastVertTheta
        
        '''for Hgrid test'''
        self.nextObj = None

    def nextVertex(self, i):
        #0 <= i <= self.n-1
        return (i+1)%self.n
        
    def prevVertex(self, i):
        #0 <= i <= self.n-1
        return (i-1+self.n)%self.n
                
    def updateBC(self):
        theta = self.pos[2]
        c = cos(theta)
        s = sin(theta)
        self.BC.c = (self.pos[0] + c*self.midBC0[0] - s*self.midBC0[1],self.pos[1] + c*self.midBC0[1] + s*self.midBC0[0])
            
    #image update prior rendering
    def updateImg(self, theta):
        self.image = pygame.transform.rotate(self.image0, -180*theta/pi)
        self.imgWidth, self.imgHeight = self.image.get_size()
        self.lastImgTheta = theta
                                      
    #vertices and normals update, done only prior to separating axis test if lastVertTheta != self.pos[2] up to a precision of THETA_TOL
    def updateRotation(self): #would this require deletion (garbage collecting) in c++ to free memory?
        theta = self.pos[2]
        c = cos(theta)
        s = sin(theta)
        for i in range(self.n):
            vert0_ix = self.vert0[i][0]
            vert0_iy = self.vert0[i][1]
            self.vert[i] = ( vert0_ix*c - vert0_iy*s , \
                             vert0_ix*s + vert0_iy*c )
        for i in range(self.nn):
            normal0_ix = self.normal0[i][0]
            normal0_iy = self.normal0[i][1]
            self.normal[i] = ( normal0_ix*c - normal0_iy*s , \
                               normal0_ix*s + normal0_iy*c )
                               
        self.lastVertTheta = theta
        
    def render(self, surface, alpha = 0):
        x = self.pos[0] - (self.speed[0] * self.world.dt)*(1-alpha)
        y = self.pos[1] - (self.speed[1] * self.world.dt)*(1-alpha)
        theta = resetTheta(self.pos[2] - (self.speed[2] * self.world.dt)*(1-alpha))
                
        #note that this has the flip flop instability if theta oscillates over 0. one solution is to extend theta range to be from -2pi to 2pi
        if abs(theta - self.lastImgTheta) > THETA_TOL:
            self.updateImg(theta)
        
        #there is a rounding off of the coordinates to integer
        surface.blit(self.image,(x - self.imgWidth/2,y - self.imgHeight/2))
            
   
    def trans(self): #With euler's method, prevPos is not needed as it can simply be calculated using the velocity
        self.pos[0] = self.pos[0]+self.speed[0]*self.world.dt
        self.pos[1] = self.pos[1]+self.speed[1]*self.world.dt
        
        #note: image rotation is done prior to rendering, to avoid unnecessarily updating image every physics time step
        if abs(self.speed[2]) > OMEGA_TOL:
            self.pos[2] = resetTheta(self.pos[2]+self.speed[2]*self.world.dt)
        else:
            self.speed[2] = 0
            
        self.updateBC()

    '''For objects that can be destroyed,
    create a new class where in __init__ define
    self.health = 1
        def damaged(self, dmg = 0):
            self.health -= dmg
            if self.health <= 0:
                self.death()

        def death(self):
            self.world.removeEntity(self.id_num)
            self.world.removeAABB(self)

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
            
  


#keeping 0\-2*pi <= theta <= 2*pi. the extended limits is to stabilize theta oscillation
def resetTheta(theta):
    if theta > 2*pi:
        theta -= 2*pi
    elif theta < -2*pi:
        theta += 2*pi
    return theta

def computeBC(vertices):
    #Calculating the minimum bounding sphere (inaccurate simple algorithm)
    #find top left and bottom right coordinates of AABB using center of mass as origin.
    x = [vertices[0][0],vertices[0][0]]
    y = [vertices[0][1],vertices[0][1]]
    for i in vertices: #the first run is actually a waste
        if i[0] < x[0]:
            x[0] = i[0]
        elif i[0] > x[1]:
            x[1] = i[0]
        if i[1] < y[0]:
            y[0] = i[1]
        elif i[1] > y[1]:
            y[1] = i[1]
                    
    #use AABB center as sphere center
    midx = (x[0]+x[1])/2.0
    midy = (y[0]+y[1])/2.0
    
    rSqMax = vertices[0][0]**2 + vertices[0][1]**2
    for i in vertices:
        rSq = i[0]**2 + i[1]**2
        if rSq > rSqMax:
            rSqMax = rSq
    
    return ((midx, midy), rSqMax**0.5)
        
#returns the coordinates of the vertices of an n-sided regular polygon of side length == sideLength,
#with the center as the origin. Note that odd n causes BC to be non-optimal for current computeBC algorithm
def regularPolygon(n, sideLength): 
    theta = 2*pi/n
    c = cos(theta)
    s = sin(theta)
    
    #r is the radius of the common circle on which all vertices lie
    r = sideLength/(2*sin(theta/2.0))
    
    #the vertices are listed CLOCKWISE, as required by the class Polygon
    #let the first vertex be the bottom-most right
    vertices0 = [(sideLength/2.0, r* cos(theta/2.0))] 

    #we add the next coordinate, which is the previous rotated by theta CLOCKWISE. Remembering y-axis points downwards,
    for i in range(n-1):
        x = vertices0[i][0]*c - vertices0[i][1]*s
        y = vertices0[i][1]*c + vertices0[i][0]*s
        vertices0.append((x,y))

    return vertices0

'''
to tuple vertices0, use tuple(). i dont because for indexing, which is what vertices and normal in Polygon are only for,
there's no noticeable difference in performance between tuple and list (list may even be faster).
tuples are more concise in terms of memory and faster for instantiation like a = (1,2,3), compared to a = [1,2,3].
'''
