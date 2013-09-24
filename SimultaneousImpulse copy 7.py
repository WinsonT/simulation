#works better for large inertia. and small mass. minimum .mass_inv/.inertia_inv is about 1000. gravity 600 px/sec. physics calculations at 60 Hz

#REMEMBER to avoid doing anything based on number of frame. do everything based on time passed

#CURRENT COLLISION (with suggested solution in the bracket)
#jitter a bit
#inevitable slide on inclined plane
#stacking objects quiver
#frictionless stacking objects results in the inactive objects above to float (group touching objects, inactivating and activating them together. and im
#probably gonna change the (de)activation system. objects colliding with active objects are activated. objects losing contacts are also activated.)
#a point of contact is owned by two bodies. modifying this point of contact must induce change in the information stored in both bodies. bodies connected
#with one(or two) point(s) of contact need not undergo the usual collision detection routine. .contacts are NOT reset every frame. .forces are reset every frame.


#TO DO


#i am in the middle of revising activation/deactivation scheme
#shock system still bad. body penetrates too deep for many-bodies collision

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

#CONSTANTS
pgs_iter = 10 #ProjectedGaussSeidel (to solve collision) maximum number of iterations
#Window size
cam_width = 960
cam_height = 720
#Map size (in tiles). because map size is always a multiple of tile size.
map_width = 15 #960 px
map_height = 12 #768 px
# Tile
tile_w = 64#128
tile_h = 64#128
time_step = 1.0/60
fps = 60


trans_speed = 25*60
drag = 1.2
deactive_time = 2
skin = 2
#9.8 m/s**2=10.5 px/frame**2 with 1 m = 64 px and 60 fps
grav = 10*60
speed_limit = 50 #speed limit below which lower gravity is applied
dis = 43 #distance between hexagons

screen_1 = pygame.display.set_mode((cam_width, cam_height), 0, 32)
pygame.display.set_caption("Collision!")
fullscreen = False
#pygame.event.set_blocked((MOUSEMOTION))
#AFRAMEPASSED = USEREVENT + 1






#ADMINISTRATIVE FUNCTIONS


#copy_list2 creates an independent copy of "list2",
#where list2 is a list in a list, i.e [[a,b], [c,d], [e,f], ...] where a, b, c, d, and etc are numbers or string.
#a, b, c, d and etc must not be an object, or else the output of the function, ie the "copy", will not be an entirely independent copy of the input.
#independent means if "list2" is changed, the "copy" remains the same, and vice versa.
#meaning "list2" and "copy" do not refer to the same object at all.
#note that list2 can be tuple2
def copy_list2(list2):
    copy = []
    for i in list2:
        a = []
        for j in i:
            a.append(j)
        copy.append(a)
    return copy

#blankcopy creates a copy of list2, with each component replaced by a new value "blank" (blank must not be an object)
#note that list2 can be tuple2
def blankcopy_list2(list2, blank):
    copy = []
    for i in list2:
        a = [blank for j in i]
#        for j in i:
 #           a.append(blank)
        copy.append(a)
    return copy

#for a list in a list in a list
#note that list3 can be tuple3
def copy_list3(list3):
    copy = []
    for i in list3:
        a = []
        for j in i:
            b = []
            for k in j:
                b.append(k)
            a.append(b)
        copy.append(a)
    return copy

#convert list in a list into tuple in tuple
def tuple2(list2):
    a = []
    for i in list2:
        a.append(tuple(i))
    return tuple(a)

def tuple3(list3):
    a = []
    for i in list3:
        b = []
        for j in i:
            b.append(tuple(j))
        a.append(tuple(b))
    return tuple(a)

#TILES
def tile_coordinate(pos, w, h): #position, width, height
    x = pos[0]
    y = pos[1]




    tile_x = int(x/tile_w)
    tile_y = int(y/tile_h)
    return tile_x, tile_y

#####experimental

class Tile():
    def __init__(self, kind):
        self.kind()

    def a (self):
        pass

    def b (self):
        pass

    def c (self):
        pass

    def d (self):
        pass

    def g (self):
        pass

    def f (self):
        pass

#and circular surfaces if needed








#class Box():
 #   def __init__(self, x, y, width = 64, height = 64):
  #      self.pos = [x, y, 0]
   #     self.width = width
    #    self.height = height
     #   self.rect = [[x - width/2, y - height/2], [x + width/2, y + height/2]]
        

#Entity is for objects that do not collide. this is mainly (if not only) for decorations that do not require vertices and normal axes.
#they can still perform actions. for ghosts, i might want to use Polygon, because they still need to interact, which requires collision detection.
#examples are dust particles, smoke.

class Entity():
    def __init__(self, pos, image, world):
        self.world = world
        
        self.health = 1 #if .health < 0, entity is removed
    
        self.pos = pos #.pos is in the form [x, y, theta], x and y are at the center of the width and height of the self.image. theta is the rotation angle
        self.prev_pos = (pos[0], pos[1], pos[2]) #.prev_pos is the position at time .world.t - .world.dt
                             
        self.speed = [0,0,0] #rotational speed in rad/s. if > 0, clockwise. angular velocity into the screen, the z direction.
        
        if image != None:
            self.image0 = pygame.image.load(image).convert_alpha() #base image, will not be rotated
            self.update_img_rotation(self.pos[2])

        else:
            self.image0 = None
            self.image = None
    #        self.width0 = 0
     #       self.height0 = 0

        
        self.count_num = 0
        

#count() counts the time passed since it is first called. it resets when time passed exceeds val. the time passed is stored in self.count_num
#count2 count3 and so on can be defined individually at each object if more variables to keep track of time is needed
    def count(self, val, time_passed): #val is an integer
        if self.count_num + time_passed < val: self.count_num += time_passed
        else: self.count_num = 0

    def process(self, event):
#        self.trans()
        pass

    def damaged(self, dmg = 0):
        self.health -= dmg
        if self.health <= 0:
            self.world.remove_entity(self.id_num)

    def trans(self): #Entity's trans doesn't need integration (for now), as an Entity moves with a constant speed

        self.pos = [self.pos[0]+self.speed[0]*self.world.dt, self.pos[1]+self.speed[1]*self.world.dt, self.pos[2]+self.speed[2]*self.world.dt]

        #-2*pi <= .pos[2] <= 2*pi. the over stretched limit is to prevent pos[2] jumping from 0.0 to 6.283185307179, disrupting the deactivation algorithm
        if self.pos[2] >= 2*pi:
            self.pos[2] -= 2*pi
            self.prev_pos = (self.prev_pos[0], self.prev_pos[1], self.prev_pos[2]-2*pi) #to prevent stutter due to anti-(temporal)aliasing
        elif self.pos[2] <= -2*pi:
            self.pos[2] += 2*pi
            self.prev_pos = (self.prev_pos[0], self.prev_pos[1], self.prev_pos[2]+2*pi) #to prevent stutter due to anti-(temporal)aliasing
            
    def update_img_rotation(self, angle):
        self.image = pygame.transform.rotate(self.image0, -180*angle/pi) #negative angle because pygame rotate() is counterclockwise



#RULES on defining vertices0:
#1. maximum y should be equal to negative minimum y, meaning the center is height/2.0 from the top and bottom of the object.
#2. the same applies to x
#3. object must be convex (all corners have angles less than 180 degrees). corners are registered clockwise
#NOTE rule number 1 and 2 may change in the future, if custom center-of-mass and/or pivot feature is added. 
class Polygon(Entity):
    def __init__ (self, pos, vertices0, image):
        Entity.__init__(self, pos, image, world)
                
        self.mass_inv = 1.0
        self.inertia_inv = self.mass_inv/1000 #0.001

        self.vert0 = vertices0
        self.vert = copy_list2(self.vert0)

        
        #find xmin xmax ymin ymax
        x = [self.vert0[0][0],self.vert0[0][0]]
        y = [self.vert0[0][1],self.vert0[0][1]]
        for i in self.vert0:
            if i[0] < x[0]:
                x[0] = i[0]
            elif i[0] > x[1]:
                x[1] = i[0]
            if i[1] < y[0]:
                y[0] = i[1]
            elif i[1] > y[1]:
                y[1] = i[1]
                
        #initial bounding volume. for optimization, it is preferred to have it initially freely-oriented/OBB (compared with the current axis-aligned/AABB).
        #alternatively, the list of points could be rotated such that at theta=0 the tightest OBB == the AABB.
        #NOTE: this is independent of the image size. the image rendering is such that the center of the image is the center of mass.
        self.width0 = x[1] - x[0]
        self.height0 = y[1] - y[0]
        self.update_wh()
        
        self.active = 1
        self.active_state = [(self.pos[0], self.pos[1], self.pos[2]), 0]
        
        #.contacts[ent2] = (length, normal, point(s)). for "shock" step purpose. or collision with springs purpose
        self.contacts = {}

        self.forces = []


#count vertices
        self.n = 0
        for i in self.vert0:
            self.n += 1
#normal0 vectors
        self.normal0 = []
        self.nn = 0 #the number of normals
        for i in range(self.n):
            
            if i == self.n - 1: k = 0
            else: k = i + 1

            a = self.vert0[i][0] - self.vert0[k][0]
            b = self.vert0[i][1] - self.vert0[k][1]

            c = (a*a + b*b)**0.5

            d = 0 #d = 0 indicates that the added normal vector has not exist in the list
#important. the comparison is more intuitively written as j[0]/j[1] == -b/a. but the denominator can be very small (e.g, 1e-16)
#multiplication is used to avoid complication. notice the sign difference compared with the similar algorithm in World.collision().
#1000 * (...) means each component is evaluated to a precision of 0.001. the last statement has also been verified (though using only 1 test case. and reasoning)
            for j in self.normal0:
                if int(1000*(j[0]*a + j[1]*b)) == 0:
                    d = 1
                    break

            if d == 0:
                self.normal0.append([ b/c, -a/c ])
                self.nn += 1

#normal vectors
        self.normal = copy_list2(self.normal0)
        self.normal_poly_update()
        
#projections of vertices to self.normals. the ordering is based on .normal.
        self.proj = []
        index = 0
        for i in self.normal0:
            self.proj.append([])
            for j in range(self.n):
                self.proj[index].append((self.vert0[j][0]*i[0] + self.vert0[j][1]*i[1], j))
            self.proj[index].sort() #for tuples in a list, sorting is based on the first component of the tuple
            index += 1
        

    def process(self, event):
#        self.contacts = {}
#        self.forces = []
        
        
        #activation and deactivation still need some work out to prevent objects from floating in mid air

        if abs(self.pos[0] - self.active_state[0][0]) >= 1 or\
           abs(self.pos[1] - self.active_state[0][1]) >= 1 or\
           abs(self.pos[2] - self.active_state[0][2]) >= 0.1:
            if self.id_num == 0:
                print 'negated'
            self.active = 1
            self.active_state = [(self.pos[0], self.pos[1], self.pos[2]), 0]
                
        elif self.active:
            
            self.active_state[1] += self.world.dt
                
            if self.active_state[1] > deactive_time:
                if self.id_num == 0:
                    print '!!!!!!!!!!!!!!'
                self.active = 0
                self.speed = [0,0,0]
                self.active_state = [(self.pos[0], self.pos[1], self.pos[2]), 0]


        if self.active: pass

#            else: self.trans()
        self.register_tile() #omit this step if the object does not collide. Thus far only class Polygon can collide                

    def trans_old(self): #Entity's trans doesn't need integration (for now), as an Entity moves with a constant speed

        self.pos = [self.pos[0]+self.speed[0]*self.world.dt, self.pos[1]+self.speed[1]*self.world.dt, self.pos[2]+self.speed[2]*self.world.dt]

        #-2*pi <= .pos[2] <= 2*pi. the over stretched limit is to prevent pos[2] jumping from 0.0 to 6.283185307179, disrupting the deactivation algorithm
        if self.pos[2] >= 2*pi:
            self.pos[2] -= 2*pi
            self.prev_pos = (self.prev_pos[0], self.prev_pos[1], self.prev_pos[2]-2*pi) #to prevent stutter due to anti-(temporal)aliasing
        elif self.pos[2] <= -2*pi:
            self.pos[2] += 2*pi
            self.prev_pos = (self.prev_pos[0], self.prev_pos[1], self.prev_pos[2]+2*pi) #to prevent stutter due to anti-(temporal)aliasing

        if self.prev_pos[2] != self.pos[2]: #note: image rotation is done prior to blitting, so as to not unnecessarily update image every physics time step
            self.update_wh()
            self.normal_poly_update_needed = 1 #so that normal_poly is updated when checking collision


    #the RK4 integrator
    def trans(self):
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
            self.prev_pos = (self.prev_pos[0], self.prev_pos[1], self.prev_pos[2]-2*pi) #to prevent stutter due to anti-(temporal)aliasing
        elif self.pos[2] <= -2*pi:
            self.pos[2] += 2*pi
            self.prev_pos = (self.prev_pos[0], self.prev_pos[1], self.prev_pos[2]+2*pi) #to prevent stutter due to anti-(temporal)aliasing

        if self.prev_pos[2] != self.pos[2]: #note: image rotation is done prior to blitting, so as to not unnecessarily update image every physics time step
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
            
        return ((force[0]*self.mass_inv, force[1]*self.mass_inv, torque*self.inertia_inv))

    def update_wh(self):
        theta = self.pos[2]
        if theta < 0:
            theta += 2*pi
        if theta > pi:
            theta -= pi
        if theta > pi/2:
            theta = pi - theta

        self.width = cos(theta)*self.width0 + sin(theta)*self.height0
        self.height = cos(theta)*self.height0 + sin(theta)*self.width0



    def register_tile(self): #note that this assumes that no wall is present. when walls r present, they have to be included in the checking
        x_min = int((self.pos[0] - self.width/2)/tile_w)
        x_max = int((self.pos[0] + self.width/2)/tile_w)
        y_min = int((self.pos[1] - self.height/2)/tile_h)
        y_max = int((self.pos[1] + self.height/2)/tile_h)

        for i in range(x_min, x_max + 1): #x_max + 1 because the output of range(a,b) ends with b-1 instead of b
            for j in range(y_min, y_max + 1): 
                #for now (and probably forever), duplicate entries is not possible
                #the no collision with allied unit needs enables simpler coding though. improve on this part later.
#the downside of this is that a tile will be evaluated many times (by each object in the tile). if it is a downside. i think it's already the most efficient.
                for k in tiling[j][i]:

                    if type(k) == int and (self.active or self.world.entity_list[k].active): #collision with another object, only if one of them is active

#note that duplicates like (a,b) and (b,a) will never occur. this is because tiling is reset every frame and this tile_register is done based on id_num.
#so in (a, b), a is always larger than b
                        if not ((self.id_num, k) in self.world.coll_checklist) and not ((self.id_num, k) in self.world.contacts): #prevents duplicates.
                            self.world.coll_checklist.append((self.id_num, k))

                    else: #collision with tile. a tile is stored as a tuple in tiling. (a, b). a indicates if there is a corner. b is the tile ID.
                        pass

                tiling[j][i].append(self.id_num)


    def normal_poly_update(self):
        for i in range(self.n):
            self.vert[i] = ( self.vert0[i][0]*cos(self.pos[2]) - self.vert0[i][1]*sin(self.pos[2]) , \
                             self.vert0[i][0]*sin(self.pos[2]) + self.vert0[i][1]*cos(self.pos[2]) )
        for i in range(self.nn):
            self.normal[i] = ( self.normal0[i][0]*cos(self.pos[2]) - self.normal0[i][1]*sin(self.pos[2]) , \
                               self.normal0[i][0]*sin(self.pos[2]) + self.normal0[i][1]*cos(self.pos[2]) )

        self.normal_poly_update_needed = 0 #.normal_poly_update_needed = 1 indicates that normal and poly should be updated when checking collision




#I Broad phase (rough) collision detection
#for rough collision detection (broad phase), we need to identify if collision is imminent for an object without checking 1 by 1 against every other object present in the map. this can be done using a list of 0s representing the tiles of the map. The tiles (roughly) occupied by an object can then be labeled with its object id number (id number in world.entity_list). while labeling the tile value with an id number, the previous tile value must be 0, or else collision is imminent.

#e.g

#000000000 an object with id number 1 and 2 occupie the area labeled 1 and 2 respiectively
#000000220 in this case, no collision is possible
#011100220
#011100000
#011100000
#000000000


#II More exact collision detection: parallel axis theorem

#note for collision against 2 objects, find the longest projection axis, undo move until no longer collide against this, then move 
#perepndicular to this axis until collide with the second object (which at first cause a shorter projection axis)
#alternatively, find the longest projection axis, move parallel to this axis until no longer collide against this as usual, then move
#perpendicular to this axis until NO LONGER collide with the second object (which at first cause a shorter projection axis)

        
#input inherit = "self" if uniform_polygon is called for a class that inherits the class "Polygon". e.g, see init function of class "Player"
def uniform_polygon(pos, n, a, image, inherit = 0): 
    n = int(n)
    theta = 2*pi/n
    r = a/(2*sin(theta/2))
    vertices0 = [(-a/2, r* cos(theta/2))] 

    for i in range(n-1):
        c,d = vertices0[i][0]*cos(theta) - vertices0[i][1]*sin(theta),vertices0[i][0]*sin(theta) + vertices0[i][1]*cos(theta)
        vertices0.append((c,d))

    #to tuple it, use tuple(). i dont use it because for indexing, which is what vertices and normal are only for, there's no noticeable
    #difference in performance between tuple and list though (list may even be faster).
    #tuples are (idk how much) more concise in terms of memory. tuples are also faster for instantiation like a = (1,2,3), compare to a = [1,2,3].
    #vertices0 = tuple(vertices0)

    if inherit != 0: return Polygon.__init__(inherit, pos, vertices0, image)
    else: return Polygon(pos, vertices0, image)

class Player(Polygon):
    def __init__ (self):
        uniform_polygon([100, 100, 0], 4, 60, "box.png", self)
        self.acc_x = [0,0] #acc_x[0] is left accel, acc_x[1] is right accel
        self.acc_y = [0,0] #acc_y[0] is up accel, acc_y[1] is down accel

    def process(self, event):
        
#        self.forces = []
        
        if event.type == KEYDOWN:
            if event.key == K_UP:
                self.acc_y[0] = 1
            elif event.key == K_DOWN:
                self.acc_y[1] = 1
            elif event.key == K_LEFT:
                self.acc_x[0] = 1
            elif event.key == K_RIGHT:
                self.acc_x[1] = 1                    


        elif event.type == KEYUP:
            if event.key == K_UP:
                self.acc_y[0] = 0                    
            elif event.key == K_DOWN:
                self.acc_y[1] = 0                    
            elif event.key == K_LEFT:
                self.acc_x[0] = 0
            elif event.key == K_RIGHT:
                self.acc_x[1] = 0

#trans_speed is the magnitude of the total velocity. when moving diagonally, each component has a magnitude less than the variable trans_speed
        vx = (self.acc_x[1] - self.acc_x[0])
        vy = (self.acc_y[1] - self.acc_y[0])
        
        v = vx**2 + vy**2
        if v == 1: v = trans_speed
        elif v == 2: v = trans_speed/2**0.5

 #       self.speed = [v*vx, v*vy, self.speed[2]]
#        self.speed[0] /= drag
 #       self.speed[1] /= drag
        self.forces.append(Gravity((v*vx, v*vy), self.mass_inv))
  #      self.speed[0] += v*vx*world.dt
   #     self.speed[1] += v*vy*world.dt

        
#        self.trans()

        if self.prev_pos[2] != self.pos[2]: #note: image rotation is done prior to blitting, so as to not unnecessarily update image every physics time step
            self.update_wh()
            self.normal_poly_update_needed = 1 #so that normal_poly is updated when checking collision

        self.register_tile()

#between Smoke and Sys_obj, there might need to be class Particle. or class Particle can be set to be above class Entity. but that means it require image.
#which opacity is harder to play with. since surfaces with per pixel alpha cannot be faded quickly enough in pygame. except you wanna use colorkey.
class Smoke(Entity):
    def __init__(self, pos, speed):

#        Entity.__init__(self, pos, "Smoke.png", world)

#modified entity init
        self.world = world
        self.health = 1
#        self.pos0 = pos
        self.pos = pos
        self.speed = [0,0,0]
        #smoke image is handled by smoke gen.
#        self.image = pygame.image.load("Smoke.png").convert_alpha() #base image, will not be rotated
 #       self.width = self.image.get_width()
  #      self.height = self.image.get_height()
        self.count_num = 0
##
        self.speed = speed
        
        self.mass_inv = 0
        self.inertia_inv = 0

    def process(self, event):
#        self.trans()
        self.count_num += self.world.dt
        if self.count_num >= 5:
            self.damaged(1)

    def trans(self):
        self.pos = [self.pos[0]+self.speed[0]*self.world.dt, self.pos[1]+self.speed[1]*self.world.dt, 0]

#Big chunk smoke deleted. DO NOT use many transformation. VERY VERY slow. I think because of the blitting etc.
#what about character's arms and legs swing and all that? havent figure it yet

#sysobj, system object, is an object that do actions per frame, but do not collide (and has no image). e.g, smoke generator (class Smoke_gen)
#img collide class
# 1    1       Polygon
# 1    0       Entity
# 0    1       Polygon
# 0    0       Entity


#class Sys_obj():
 #   def __init__(self, world):
  #      self.world = world
    
   #     self.health = 1 #if self.health < 0, object is removed
        
    #    self.count_num = 0

#count() counts the time passed since it is first called. it resets when time passed exceeds val. the time passed is stored in self.count_num
#count2 count3 and so on can be defined individually at each object if more variables to keep track of time is needed
#    def count(self, val, time_passed): #val is an integer
 #       if self.count_num + time_passed < val: self.count_num += time_passed
  #      else: self.count_num = 0

   # def process(self, event):
    #    pass

#    def damaged(self, dmg = 0):
 #       self.health -= dmg
  #      if self.health <= 0:
   #         self.world.remove_sysobj(self.id_num)



#this NEEDS to be remade. the problem with this is that it cannot make a thick enough smoke/fog with a reasonble amount of particles.
#an idea is to make the image more opaque. but this would require pygame to be able to make it fade, or else it will look really bad.
#the key in making good smoke is shape dynamic, which is currently achieved by the massive amount of particles - just like real smoke.

#actually current smoke using 9% layer opacity in photoshop is quite okay (though 50 blits could be used for so many
#more useful purpose). a source of 10x10 looks decent with 50 particles.
#only if the whole image can be made to fade
class Smoke_gen(Entity):
    def __init__ (self, prop, lifetime = 0):
        pos = [0,0,0]
        Entity.__init__(self, pos, None, world)
        #x and y is the position of the top left of the box. in contrast to x and y in entities, which represent the center of their images.
        self.x, self.y, self.range_x, self.range_y = prop
        self.lifetime = lifetime
        
        #density is number of smoke particles per unit time (second) per unit area (in tens pixel square).
        self.density = 0.1

        #speed of smoke should be proportional to temperature. the hotter the faster. also, smoke due to cold objects is slow, though might be dense.
        #for now, the speed is fixed. also remember the direction of average velocity is against gravity
        self.v = 30 #60 for hot fire. with density 0.005. 20 for cold object. density 0.002
        
        temp = pygame.image.load("Smoke.png").convert_alpha()
        self.w = temp.get_width()
        self.h = temp.get_height()
        #images of a smoke particle with 1, 2, and 3x magnification
        self.image_list = (temp, pygame.transform.scale2x(temp), pygame.transform.scale(temp, (3*self.w, 3*self.h)))
        #for i in self.image_list:
         #   locked_image.append(i)
         
        self.mass_inv = 0
        self.inertia_inv = 0

    def process(self, event):
        if self.lifetime > 0:
            self.damaged(self.world.dt/self.lifetime)
        
        self.count_num += self.density * self.range_x * self.range_y * self.world.dt #coding takes decimal places into account

        while self.count_num >= 1:
            self.count_num -= 1

            #direction = random()*2*pi

            #Smoke(pos, speed)
            smoke = Smoke([randint(self.x, self.x + self.range_x), randint(self.y, self.y + self.range_y), 0],\
            [randint(-self.v, self.v), -randint(0, self.v), 0])
            #[self.v*(0.66+random()/3) * cos(direction),self.v*(0.66+random()/3) * sin(direction), 0])
            
            a = randint(1,2) #thus far the best combination. might change of course.
            smoke.image = self.image_list[a]
            smoke.width = (a+1)*self.w
            smoke.height = (a+1)*self.h
            
            self.world.add_entity(smoke)
            #self.world.add_sysobj(smoke)
 
    

#create a map (currently, the camera is the map)
class World():
    def __init__(self):

        self.scr = 0

        self.clock = pygame.time.Clock()
        self.current_t = 0 #the time up to which the game is trying to simulate. starts with time()
        self.frame_t = 0 #the time elapsed since the last frame
        self.t = 0 #the time up to which the game has simulated. starts with 0
        self.dt = time_step #time step for calculations
        self.remainder_t = 0

        self.entity_list = {}
        self.next_id_num = 0
#        self.sysobj_list = {}
 #       self.next_sysobj_id_num = 0

        self.coll_checklist = []
        self.contacts = {}
        self.shocks = {} #contacts with immovables
        self.immovables = []

#it is a dictionary, since i cannot use a list. this is because when an entity dies, when it is removed from the list, the index of other 
#entities will change, making it difficult to organize.
        
    def add_entity(self, entity):
        self.entity_list[self.next_id_num] = entity
        entity.id_num = self.next_id_num
        self.next_id_num += 1
        
#    def add_sysobj(self, sysobj):
 #       self.sysobj_list[self.next_sysobj_id_num] = sysobj
  #      sysobj.id_num = self.next_sysobj_id_num
   #     self.next_sysobj_id_num += 1

    def remove_entity(self, id_num):
        del self.entity_list[id_num]

#    def remove_sysobj(self, id_num):
 #       del self.sysobj_list[id_num]

    def process(self, event):

        d = 0

        for i in self.entity_list.values():

#            i.contacts = {}  #**
                
            if i.active:
                
                i.prev_pos = (i.pos[0], i.pos[1], i.pos[2])

#                if i.id_num == 0: print i.pos
                
                #gravity
                if i.mass_inv != 0:
                    if i.speed[1] > 0 and i.speed[1] < speed_limit:
                        i.forces.append(Gravity((0,grav*0.7), i.mass_inv))
                    else:
                        i.forces.append(Gravity((0,grav), i.mass_inv))
                #gravity
     #           if i.mass_inv != 0:
      #              if i.speed[1] > 0 and i.speed[1] < speed_limit:
       #                 i.speed[1] += grav*0.7*self.dt
        #            else:
         #               i.speed[1] += grav*self.dt

            else: d += 1
            
            i.process(event)
            i.trans()
            i.forces = []
                        

            
        

        print d
        
#collision routine
        for p in xrange(pos_iter):
            

            for i in self.contacts.values():
                self.contact_update(i, self.contacts)

            for i in self.shocks.values():
                self.contact_update(i, self.shocks)


            for ent1, ent2 in self.coll_checklist:
                self.coll_detect(self.entity_list[ent1], self.entity_list[ent2])                

            for i in self.contacts.values():
                self.pos_response(i)
    #        for i in self.shocks.values():
     #           self.pos_response(i)
                
            for i in self.immovables:

                self.preshock(i)

        for p in xrange(vel_iter):
            for i in self.contacts.values():
                self.vel_response(i)
            for i in self.shocks.values():
                self.vel_response(i)

        self.coll_checklist = []
        self.immovables = []
#        self.contacts = {}#**
 #       self.shocks = {}  #**

     #   for i in self.entity_list.values():
      #      if i.active:
       #         i.trans()
        #        i.forces = []


    def render(self):
#        screen_1.fill((0,0,0))#((255,255,255,255)) #alternately, u can blit BG here
#        screen_1.fill((255,255,255))
        screen_1.fill((155,155,155))

#anti-(temporal) aliasing
        alpha = self.remainder_t/self.dt
        #alpha = 0
        
#image system needs to be changed. an object needs to be able to use multiple images as its full image. e.g, a person has head, legs, hands, and body.
#also, macroscopic object might want to be divided into smaller ones, especially if a side of it is concave (my collision detection is for convex object).
#this might allow the merge of sysobj and entity class,
#because the main difference between the two is only (as far as i can remember) that sysobj do not have images. blitting will probably inside a method of
#class World. this also allows blitting after world.process()

        for i in self.entity_list.values():

            if i.image != None: i.update_img_rotation((i.pos[2]-i.prev_pos[2])*alpha + i.prev_pos[2] )
            if i.image != None: screen_1.blit(i.image,(int((i.pos[0]-i.prev_pos[0])*alpha + i.prev_pos[0]) - i.image.get_width()/2.0,\
                                                       int((i.pos[1]-i.prev_pos[1])*alpha + i.prev_pos[1]) - i.image.get_height()/2.0)  )

    
    
    def contact_update(self, contact, contacts):            
        

        
        ent1 = self.entity_list[contact[0][0]]
        ent2 = self.entity_list[contact[0][1]]
        proj_vec = contact[1]

        if ent1.normal_poly_update_needed: ent1.normal_poly_update()
        if ent2.normal_poly_update_needed: ent2.normal_poly_update()

#|     body2   |
#p3  ^         |
# \  | normal  |
#t1=====t2    /   ---> tangent
#|   \p1_|__p2
#| body1 |

        if proj_vec[1][0] == ent1:
            body1 = ent1
            body2 = ent2
            if proj_vec[1][2] == 1:
                limit = -ent1.proj[proj_vec[1][1]][0][0] #if ent2 max proj is less than this, no more collision
                normal = (-ent1.normal[proj_vec[1][1]][0], -ent1.normal[proj_vec[1][1]][1])
                tangent = (ent1.normal[proj_vec[1][1]][1], -ent1.normal[proj_vec[1][1]][0])
                t1 = ent1.proj[proj_vec[1][1]][0][1]
                t2 = ent1.proj[proj_vec[1][1]][1][1]

                
            else:
                limit = ent1.proj[proj_vec[1][1]][ent1.n-1][0]
                normal = ent1.normal[proj_vec[1][1]]
                tangent = (-ent1.normal[proj_vec[1][1]][1], ent1.normal[proj_vec[1][1]][0])
                t1 = ent1.proj[proj_vec[1][1]][ent1.n-1][1]
                t2 = ent1.proj[proj_vec[1][1]][ent1.n-2][1]
        
        #normal1 is the same as the evaluated proj_vec[1]. for shock step purpose
            normal1 = (-normal[0], -normal[1])
        
        else:
            body1 = ent2
            body2 = ent1
            if proj_vec[1][2] == 1:
                limit = ent2.proj[proj_vec[1][1]][ent2.n-1][0]
                normal = ent2.normal[proj_vec[1][1]]
                tangent = (-ent2.normal[proj_vec[1][1]][1], ent2.normal[proj_vec[1][1]][0])
                t1 = ent2.proj[proj_vec[1][1]][ent2.n-1][1]
                t2 = ent2.proj[proj_vec[1][1]][ent2.n-2][1]

            else:
                limit = -ent2.proj[proj_vec[1][1]][0][0]
                normal = (-ent2.normal[proj_vec[1][1]][0], -ent2.normal[proj_vec[1][1]][1])
                tangent = (ent2.normal[proj_vec[1][1]][1], -ent2.normal[proj_vec[1][1]][0])
                t1 = ent2.proj[proj_vec[1][1]][0][1]
                t2 = ent2.proj[proj_vec[1][1]][1][1]

        #normal1 is the same as the evaluated proj_vec[1]. for shock step purpose
            normal1 = normal


        point = body1.vert[t1]
        t1t = point[0]*tangent[0] + point[1]*tangent[1]

        point = body1.vert[t2]
        t2t = point[0]*tangent[0] + point[1]*tangent[1]
        
        if t2t < t1t:
            t1t, t2t = t2t, t1t
            t1, t2 = t2, t1
            
        pos_12n = normal[0]*body2.pos[0] + normal[1]*body2.pos[1] - normal[0]*body1.pos[0] - normal[1]*body1.pos[1]
        pos_12t = tangent[0]*body2.pos[0] + tangent[1]*body2.pos[1] - tangent[0]*body1.pos[0] - tangent[1]*body1.pos[1]
        
        #limit == 0 is now at the center of body2
        limit -= pos_12n
        
        #single point of contact

        if proj_vec[2][1] == 0:
            p1 = proj_vec[2][0][1]
            if p1 == body2.n-1:
                p2 = p1-1
                p3 = 0
            elif p1 == 0:
                p2 = body2.n - 1
                p3 = p1+1
            else:
                p2 = p1-1
                p3 = p1+1
        #2 points of contact
        else:
            j = []
            for i in proj_vec[2]:
                if i[0] == body2:
                    j.append(i[1])
            p1 = j[0]
            p2 = j[1]
            p3 = None
            

        point = body2.vert[p1]
        p1n = point[0]*normal[0] + point[1]*normal[1]

        point = body2.vert[p2]
        p2n = point[0]*normal[0] + point[1]*normal[1]
        
        #we want p1n < p2n.
        if p1n > p2n:
            p1, p2 = p2, p1
            p1n, p2n = p2n, p1n

        #after this, p3 and p3n will not be used
        if p3 != None:
            point = body2.vert[p3]
            p3n = point[0]*normal[0] + point[1]*normal[1]
            if p3n < p1n:
                p2 = p1
                p2n = p1n
                p1 = p3
                p1n = p3n

            elif p3n < p2n:
                p2 = p3
                p2n = p3n

        #if no longer in contact
        if p1n > limit:

            del contacts[contact[0]]
            del body1.contacts[body2]
            del body2.contacts[body1]

            return

        #else if only 1 point of contact
        elif p2n - p1n > 0.1 or p2n > limit:
            point = body2.vert[p1]
            p1t = point[0]*tangent[0] + point[1]*tangent[1]
            if t1 < p1t and p1t < t2:
                
                if ent1.mass_inv == 0:
                    if not(ent1 in self.immovables): self.immovables.append(ent1)
                elif ent2.mass_inv == 0:
                    if not(ent2 in self.immovables): self.immovables.append(ent2)
                
                proj_vec = (limit - p1n, proj_vec[1], ((body2, p1),0))
                
                n1 = (proj_vec[0], normal1)#, proj_vec[2])
                n2 = (proj_vec[0], (-normal1[0], -normal1[1]))#, proj_vec[2])
                
                contacts[contact[0]] = (contact[0], proj_vec)

                ent1.contacts[ent2] = n1
                ent2.contacts[ent1] = n2
                
            else:
                #optional: add the ent1 ent2 pair to .coll_checklist to process the collision directly
                self.coll_checklist.append((ent1.id_num,ent2.id_num))
                
                del contacts[contact[0]]
                del body1.contacts[body2]
                del body2.contacts[body1]

                return

        #else there are 2 points of contact
        else:      
            point = body2.vert[p1]
            p1t = point[0]*tangent[0] + point[1]*tangent[1]
            point = body2.vert[p2]
            p2t = point[0]*tangent[0] + point[1]*tangent[1]
            
            a = [(t1t, (body1, t1)), (t2t, (body1, t2)), (p1t, (body2, p1)), (p2t, (body2, p2))]
            a.sort()

            if ent1.mass_inv == 0:
                if not(ent1 in self.immovables): self.immovables.append(ent1)
            elif ent2.mass_inv == 0:
                if not(ent2 in self.immovables): self.immovables.append(ent2)

            proj_vec = (limit - p1n, proj_vec[1], (a[1][1], a[2][1], a[0][1], a[3][1]))

            n1 = (proj_vec[0], normal1)#, proj_vec[2])
            n2 = (proj_vec[0], (-normal1[0], -normal1[1]))#, proj_vec[2])

            contacts[contact[0]] = (contact[0], proj_vec)

            ent1.contacts[ent2] = n1
            ent2.contacts[ent1] = n2

    
        
    def coll_detect(self, ent1, ent2):
    
        
        if ent1.normal_poly_update_needed: ent1.normal_poly_update()
        if ent2.normal_poly_update_needed: ent2.normal_poly_update()
         
        #broader phase: AABB (use x and y axis). useful if many objects are cramming in the same tile, typically if the objects are small.
        if ent1.width/2.0 < ent2.pos[0] - ent2.width/2.0 - ent1.pos[0] or\
           ent2.width/2.0 < ent1.pos[0] - ent1.width/2.0 - ent2.pos[0] or\
           ent1.height/2.0 < ent2.pos[1] - ent2.height/2.0 - ent1.pos[1] or\
           ent2.height/2.0 < ent1.pos[1] - ent1.height/2.0 - ent2.pos[1]:
            return 0

        
        #contains indices of normals of ent2 that are equal to one of ent1's normals.
        checked_list = []
        
        proj_vec = ("",["",""], ("",""))

        #iterate through all normal vectors
        
        #SUGGESTION:
        #instead of preparing and sorting the projections (each one of them) before checking if the entities are overlapping,
        #straightaway only check if the min/max projections of the other entities is penetrating.
        #and instead of projecting all vertices and find the max/min, project 1, and check the next vertex against the desired normal vector
        #looking for a turning point in the projection value to find min/max projection, also checking if point penetrate, making use of whether
        #we are currently on the min or max side of the other polygon
        #this also reminds me to fix ent.proj. instead of keeping all vertices, only store min and max projections
        #instead of checking if ent1 and ent2 have common normals, straightaway pick 1 normal and check for collision.
        #to prevent redundance, just make sure that the current normal vector is not the same as any of the checked normal.
        for i in range(ent1.nn+ent2.nn):
        
            if i < ent1.nn:
                k = ent1.normal[i]
                
                ent1_proj = ent1.proj[i]
            
                d = 0
                j = 0
            
                while j < ent2.nn:
#the comparison is more intuitively written as j[1]/j[0] == i[1]/i[0]. but the denominator can be very small (e.g, 1e-16)
#multiplication is used to avoid complication. 1000 * (...) means each component is evaluated to a precision of 0.001
#the last statement has also been verified (though using only 1 test case. and reasoning)
                    if not(j in checked_list) and int(1000*(ent2.normal[j][1]*k[0] - k[1]*ent2.normal[j][0])) == 0:
                        checked_list.append(j)
                        d = 1
                        #if the two normals are anti-parallel.
                        #wouldnt it be easier to just dot the two vectors? parallel if > 0
                        if int( abs(ent2.normal[j][1]-k[1]) + abs(ent2.normal[j][0]-k[0]) ) > 0:

                            ent2_proj = []
                            for m in range(ent2.n):
                                ent2_proj.append((-ent2.proj[j][ent2.n-(m+1)][0], ent2.proj[j][ent2.n-(m+1)][1]))
                    #else they are parallel
                        else:
                            ent2_proj = ent2.proj[j]
                        break #in C, to not implement break, use additional condition on the for loop to get out as soon as d != ""
                    j += 1


                if d == 0:
                    ent2_proj = []
                    for m in range(ent2.n):
                        ent2_proj.append( (ent2.vert[m][0]*k[0] + ent2.vert[m][1]*k[1], m) )
                    ent2_proj.sort()

            elif not((i-ent1.nn) in checked_list):
                k = ent2.normal[i-ent1.nn]
                ent2_proj = ent2.proj[i-ent1.nn]
                
                ent1_proj = []
                for m in range(ent1.n):
                    ent1_proj.append( (ent1.vert[m][0]*k[0] + ent1.vert[m][1]*k[1], m) )
                ent1_proj.sort()
            
            
            else:
                continue
                                

            pos_proj12 = k[0]*ent2.pos[0] + k[1]*ent2.pos[1]-k[0]*ent1.pos[0] - k[1]*ent1.pos[1]


#if not overlapping
            if ent1_proj[ent1.n - 1][0] < ent2_proj[0][0]  +pos_proj12 or\
               ent2_proj[ent2.n - 1][0] < ent1_proj[0][0]  -pos_proj12:
                return 0

#else it is overlapping. then if ent2 is on the "right" where "right" is the direction vector k is pointing
            elif ent1_proj[ent1.n - 1][0] - ent2_proj[0][0] < ent2_proj[ent2.n - 1][0] - ent1_proj[0][0] +2*pos_proj12:
                length = ent1_proj[ent1.n - 1][0] - ent2_proj[0][0] -pos_proj12# - 1
                k2 = k
                k = (-k[0], -k[1])
                if i >= ent1.nn:
                    direction = (ent2, i-ent1.nn, -1)
                else:
                    direction = (ent1, i, -1)

#find the point of contact
                if abs(ent1_proj[ent1.n - 1][0] - ent1_proj[ent1.n - 2][0]) > 0.1: #point is a vertice of ent1
                    point = ((ent1, ent1_proj[ent1.n - 1][1]),0) #the point[1] == 0 indicates there's only 1 point of contact 
                elif abs(ent2_proj[0][0] - ent2_proj[1][0]) > 0.1: #point is a vertice of ent2
                    point = ((ent2, ent2_proj[0][1]),0)
                else: #collision is between flat surfaces. pick a point. by using random, it simulates the situation where the force is uniformly distributed.
                    a = []
                    for m, n in (\
                              ((ent1.vert[ent1_proj[ent1.n-1][1]][0]+ent1.pos[0],ent1.vert[ent1_proj[ent1.n-1][1]][1]+ent1.pos[1]), (ent1,ent1_proj[ent1.n-1][1])),\
                              ((ent1.vert[ent1_proj[ent1.n-2][1]][0]+ent1.pos[0],ent1.vert[ent1_proj[ent1.n-2][1]][1]+ent1.pos[1]), (ent1,ent1_proj[ent1.n-2][1])),\
                              ((ent2.vert[ent2_proj[0][1]][0] + ent2.pos[0], ent2.vert[ent2_proj[0][1]][1] + ent2.pos[1]), (ent2, ent2_proj[0][1])),\
                              ((ent2.vert[ent2_proj[1][1]][0] + ent2.pos[0], ent2.vert[ent2_proj[1][1]][1] + ent2.pos[1]), (ent2, ent2_proj[1][1])) ):
                        a.append((m[0]*k[1] - m[1]*k[0], n))
                    a.sort()
                    point = (a[1][1], a[2][1], a[0][1], a[3][1]) #the first two are the points of contact. the last two are the remaining points involved
#                    if random() < 0.5: point = (2, (a[1][1],a[2][1]))
 #                   else: point = (2, (a[2][1], a[1][1]))

                
#if ent2 is on the left
            else:
                length = ent2_proj[ent2.n - 1][0] - ent1_proj[0][0] + pos_proj12# + 1
                k2 = (-k[0], -k[1])
                if i >= ent1.nn:
                    direction = (ent2, i-ent1.nn, 1)
                else:
                    direction = (ent1, i, 1)


#find the point of contact
                if abs(ent1_proj[0][0] - ent1_proj[1][0]) > 0.1: #point is a vertice of ent1
                    point = ((ent1, ent1_proj[0][1]),0)
                elif abs(ent2_proj[ent2.n - 1][0] - ent2_proj[ent2.n - 2][0]) > 0.1: #point is a vertice of ent2
                    point = ((ent2, ent2_proj[ent2.n - 1][1]),0)
                else: #collision is between flat surfaces. pick a point. by using random, it simulates the situation where the force is uniformly distributed. 
                    a = []
                    for m, n in (
                              ((ent2.vert[ent2_proj[ent2.n-1][1]][0]+ent2.pos[0],ent2.vert[ent2_proj[ent2.n-1][1]][1]+ent2.pos[1]), (ent2,ent2_proj[ent2.n-1][1])),\
                              ((ent2.vert[ent2_proj[ent2.n-2][1]][0]+ent2.pos[0],ent2.vert[ent2_proj[ent2.n-2][1]][1]+ent2.pos[1]), (ent2,ent2_proj[ent2.n-2][1])),\
                              ((ent1.vert[ent1_proj[0][1]][0] + ent1.pos[0], ent1.vert[ent1_proj[0][1]][1] + ent1.pos[1]), (ent1, ent1_proj[0][1])),\
                              ((ent1.vert[ent1_proj[1][1]][0] + ent1.pos[0], ent1.vert[ent1_proj[1][1]][1] + ent1.pos[1]), (ent1, ent1_proj[1][1])) ):
                        a.append((m[0]*k[1] - m[1]*k[0], n))
                    a.sort()
                    point = (a[1][1], a[2][1], a[0][1], a[3][1])
#                    if random() < 0.5: point = (2, (a[1][1],a[2][1]))
 #                   else: point = (2, (a[2][1], a[1][1]))


#only replace proj_vec if abs(length) is smaller than proj_vec[0], since for all k in axes, we want the minimum overlap.
#pushing ent1 in the direction of proj_vec[1] separates the objects
            if proj_vec[0] == "" or proj_vec[0] > length:
                proj_vec = (length, direction, point)
                proj_vec2 = (length, (direction[0], direction[1], -direction[2]), point)
                
                n1 = (length, k)#, point)
                n2 = (length, k2)#, point)

        
        ent1.contacts[ent2] = n1
        ent2.contacts[ent1] = n2

        
        #contact or shock (collision against an immovable). ent1 and ent2 cannot be both immovables
        #it is important to note that ent1.id_num is always larger than ent2.id_num.
        if ent1.mass_inv == 0:
            self.shocks[(ent1.id_num, ent2.id_num)] = ((ent1.id_num, ent2.id_num), proj_vec)
            if not(ent1 in self.immovables): self.immovables.append(ent1)

        elif ent2.mass_inv == 0:
            self.shocks[(ent1.id_num, ent2.id_num)] = ((ent1.id_num, ent2.id_num), proj_vec)
            if not(ent2 in self.immovables): self.immovables.append(ent2)            

        else:
            self.contacts[(ent1.id_num, ent2.id_num)] = ((ent1.id_num, ent2.id_num), proj_vec)

        
        return 0


    def pos_response(self, contact):
        ent1 = self.entity_list[contact[0][0]]
        ent2 = self.entity_list[contact[0][1]]
        
        proj_vec = contact[1]
        normal = (proj_vec[1][0].normal[proj_vec[1][1]][0]*proj_vec[1][2], proj_vec[1][0].normal[proj_vec[1][1]][1]*proj_vec[1][2])
        
            #world.add_entity(Smoke_gen((int(proj_vec[2][0]-5), int(proj_vec[2][1]-5), 10, 10), 0.1))

#for sticky substance, set e = 0, and no need to update position to remove overlap when colliding! :D

#we are projecting the bodies not until they completely separate (note the proj_vec[0]-skin). this is done to continuously record contact between objects
#to reduce jittering and allow more stacking. we do not modify position if proj_vec <= 1 or else objects will stick to each other.

#to move the object, undo their movement (instead of directly moving the object. this method causes inevitable sliding on inclined plane). but undo movement
#is equivalent to sweep test... will be back to this matter.
        if proj_vec[0] > skin:
            ent1.pos[0] += (proj_vec[0]-skin)*normal[0]*ent1.mass_inv/(ent1.mass_inv + ent2.mass_inv)
            ent1.pos[1] += (proj_vec[0]-skin)*normal[1]*ent1.mass_inv/(ent1.mass_inv + ent2.mass_inv)
            ent2.pos[0] -= (proj_vec[0]-skin)*normal[0]*ent2.mass_inv/(ent1.mass_inv + ent2.mass_inv)
            ent2.pos[1] -= (proj_vec[0]-skin)*normal[1]*ent2.mass_inv/(ent1.mass_inv + ent2.mass_inv)

    def preshock(self, immovable):
    
        checked = []
        
        for ent in immovable.contacts:
            overlap = immovable.contacts[ent]
            if ent.mass_inv != 0 and overlap[0] > skin:
                ent.pos[0] -= (overlap[0]-skin)*overlap[1][0]
                ent.pos[1] -= (overlap[0]-skin)*overlap[1][1]
                
                checked.append(ent)
        for ent in immovable.contacts:
            overlap = immovable.contacts[ent]
            if ent.mass_inv != 0 and overlap[0] > skin:
                self.shock(ent, overlap, checked)
                
#        immovable.contacts = {}


    def shock(self, immovable, overlap, checked):
        #normal inv is the direction opposite to the normal force ent feels
        #overlap is how much ent should be moved, in the format of proj_vec. overlap[2] is not used though, so it is present only in the first depth of recursion
        #overlap[0] is length, overlap[1] is the direction opposite to how much ent should be moved

        for ent in immovable.contacts:
            normal_inv = immovable.contacts[ent][1]
            if ent.mass_inv != 0 and not(ent in checked) and overlap[0] > skin:
                length = (overlap[0]-skin)*(normal_inv[0]*overlap[1][0] + normal_inv[1]*overlap[1][1])
                if length > 0:
                    ent.pos[0] -= length*normal_inv[0]
                    ent.pos[1] -= length*normal_inv[1]
            
                    checked.append(ent)

        for ent in immovable.contacts:
            normal_inv = immovable.contacts[ent][1]
            if ent.mass_inv != 0 and not(ent in checked) and overlap[0] > skin and (normal_inv[0]*overlap[1][0] + normal_inv[1]*overlap[1][1]) > 0:
                self.shock(ent, (length+skin, normal_inv), checked)#length + 1

#        immovable.contacts = {}


            


    def vel_response(self, contact):

        ent1 = self.entity_list[contact[0][0]]
        ent2 = self.entity_list[contact[0][1]]
        
        proj_vec = contact[1]
        normal = (proj_vec[1][0].normal[proj_vec[1][1]][0]*proj_vec[1][2], proj_vec[1][0].normal[proj_vec[1][1]][1]*proj_vec[1][2])


        #coefficient of restitution/bounciness. perhaps use values from both entities and multiply them to calculate e.

        #NOTE: small friction ruins stacking using this method 

        f = 0.5#0.2 
        e = 0.5#0.5

        if proj_vec[2][1] == 0:
            poc = (proj_vec[2][0][0].vert[proj_vec[2][0][1]][0] + proj_vec[2][0][0].pos[0], proj_vec[2][0][0].vert[proj_vec[2][0][1]][1] + proj_vec[2][0][0].pos[1])
        else:

            poc1 = (proj_vec[2][0][0].vert[proj_vec[2][0][1]][0] + proj_vec[2][0][0].pos[0], proj_vec[2][0][0].vert[proj_vec[2][0][1]][1] + proj_vec[2][0][0].pos[1])
            poc2 = (proj_vec[2][1][0].vert[proj_vec[2][1][1]][0] + proj_vec[2][1][0].pos[0], proj_vec[2][1][0].vert[proj_vec[2][1][1]][1] + proj_vec[2][1][0].pos[1])
            
            j = []
            
            for poc in (poc1, poc2):

            # r_1p is the vector from the center of ent1 to p, poc, the point of contact. 
                r_1p = ((poc[0] - ent1.pos[0]), (poc[1] - ent1.pos[1]))
                r_2p = ((poc[0] - ent2.pos[0]), (poc[1] - ent2.pos[1]))

            #speeds of point of contact. v_1p = v_1 + w_1 cross r_1p. v_1 = center of mass velocity. w_1 = angular velocity
                v_1p = (ent1.speed[0] - ent1.speed[2]*r_1p[1], ent1.speed[1] + ent1.speed[2]*r_1p[0])
                v_2p = (ent2.speed[0] - ent2.speed[2]*r_2p[1], ent2.speed[1] + ent2.speed[2]*r_2p[0])

            #v_apbpn = (v_2p - v_1p) dot n > 0.
                v_apbpn = (v_2p[0]-v_1p[0])*normal[0]+(v_2p[1]-v_1p[1])*normal[1]
                if v_apbpn <= 0:
                    v_apbpn = 0
                elif v_apbpn < 5*grav*self.dt:
                    e = 0

                r_1pf = -r_1p[0] * normal[1] + r_1p[1] * normal[0]
                r_2pf = -r_2p[0] * normal[1] + r_2p[1] * normal[0]
                m_inv1 = (ent1.mass_inv + ent2.mass_inv + ent1.inertia_inv*r_1pf**2 + ent2.inertia_inv*r_2pf**2) # > 0

                j.append( (1+e)*(v_apbpn)/m_inv1)
            if j[0] == 0 and j[1] == 0:
                return #there's a return here!!!!!!!!!!!!!!!!!!!!!!!!
            else:
                poc = [(poc1[0]*j[0] + poc2[0]*j[1])/(j[0]+j[1]), (poc1[1]*j[0] + poc2[1]*j[1])/(j[0]+j[1]) ]

            
        #let n be the direction of normal force felt by ent1. n = normal. let f be the unit vector n rotated 90 degrees clockwise. f is given by
        #(-normal[1], normal[0]). n cross f = z, the unit vector into the screen. z is also the direction of positive angular velocities.

        # r_1p is the vector from the center of ent1 to p, poc, the point of contact. 
        r_1p = ((poc[0] - ent1.pos[0]), (poc[1] - ent1.pos[1]))
        r_2p = ((poc[0] - ent2.pos[0]), (poc[1] - ent2.pos[1]))


        #speeds of point of contact. v_1p = v_1 + w_1 cross r_1p. v_1 = center of mass velocity. w_1 = angular velocity
        v_1p = (ent1.speed[0] - ent1.speed[2]*r_1p[1], ent1.speed[1] + ent1.speed[2]*r_1p[0])
        v_2p = (ent2.speed[0] - ent2.speed[2]*r_2p[1], ent2.speed[1] + ent2.speed[2]*r_2p[0])


        #the following algorithm is based on an unrealistic model of friction. the direction of real friction depends on the value of relative velocity.
        #However in this model, the normal force also determines the direction of the friction (in reality, this only happens if the relative velocity
        #is 0). The friction, assumed to be constant throughout the collision, is such that it minimizes sliding after the collision response.
        #update: friction now affects normal impulse such that v_apbpn (see below) after collision is e * v_apbpn before collision. this is to avoid adding
        #        energy to the system.

        #v_apbpn = (v_2p - v_1p) dot n > 0.
        v_apbpn = (v_2p[0]-v_1p[0])*normal[0]+(v_2p[1]-v_1p[1])*normal[1]
        #v_apbpf = (v_2p - v_1p) dot f. 
        v_apbpf = -(v_2p[0]-v_1p[0])*normal[1]+(v_2p[1]-v_1p[1])*normal[0]

        #if v_apbpn <= 0, objects are trying to move away from one another. in this case, there is no collision.
        #do not immediately return, because there might be two points of contact. v_apbpn = 0 results in j = k = 0 (see below)
        if v_apbpn <= 0:
            return #note the return here!!!!!!!!!!!!
        elif v_apbpn < 5*grav*self.dt:
            e = 0

        #r_1pf = r_1p dot f.
        r_1pf = -r_1p[0] * normal[1] + r_1p[1] * normal[0]
        r_2pf = -r_2p[0] * normal[1] + r_2p[1] * normal[0]
        #r_1pn = r_1p dot n.
        r_1pn = r_1p[0] * normal[0] + r_1p[1] * normal[1]
        r_2pn = r_2p[0] * normal[0] + r_2p[1] * normal[1]

        #m_inv1 and m_inv2 equals to 0 only if ent1 and ent2 are both immovablem, i.e. their .mass_inv and .inertia_inv are both 0.
        m_inv1 = (ent1.mass_inv + ent2.mass_inv + ent1.inertia_inv*r_1pf**2 + ent2.inertia_inv*r_2pf**2) # > 0
        m_inv2 = (ent1.mass_inv + ent2.mass_inv + ent1.inertia_inv*r_1pn**2 + ent2.inertia_inv*r_2pn**2) # > 0
        m_inv3 = (r_1pn*r_1pf*ent1.inertia_inv + r_2pn*r_2pf*ent2.inertia_inv)

        #calculate the required frictional impulse
        k = ((1+e)*v_apbpn*m_inv3 + v_apbpf*m_inv1)/(m_inv1*m_inv2 - m_inv3**2)

        #compare with the maximum frictional impulse. notice the direction of l depends on k, the required friction, which depends on the normal impulse.
        if k != 0:
            l = f*(1+e)*v_apbpn/(m_inv1*abs(k)/k - m_inv3*f) #this method allows l to be in the direction opposite to k, which is very very unrealistic
            if abs(k) > abs(l):
                k = abs(k*l)/k #so this measure is needed to ensure l to point in the k direction 
#                k = l

        #j is the impulse on ent1 in the direction of n. since we know n is the correct direction for impulse on ent1, j should always be > 0
        j = ((1+e)*v_apbpn+k*m_inv3)/m_inv1        

        
        #update speed change due to normal force
        #ent1.speed += ent1.mass_inv * j*n.
        ent1.speed[0] += normal[0]*ent1.mass_inv*j 
        ent1.speed[1] += normal[1]*ent1.mass_inv*j
        #update speed change due to dynamic friction 
        ent1.speed[0] -= normal[1]*ent1.mass_inv*k
        ent1.speed[1] += normal[0]*ent1.mass_inv*k

        #update speed change due to normal force
        ent2.speed[0] -= normal[0]*ent2.mass_inv*j
        ent2.speed[1] -= normal[1]*ent2.mass_inv*j
        #update speed change due to dynamic friction
        ent2.speed[0] += normal[1]*ent2.mass_inv*k
        ent2.speed[1] -= normal[0]*ent2.mass_inv*k

        
        #ent1.speed[2] += (r_1p cross n) * ent1.inertia_inv * j. 
        ent1.speed[2] -= r_1pf*ent1.inertia_inv*j
        ent2.speed[2] += r_2pf*ent2.inertia_inv*j
        #ent1.speed[2] += (r_1p cross f) * ent1.inertia_inv * k. 
        ent1.speed[2] += r_1pn*ent1.inertia_inv*k
        ent2.speed[2] -= r_2pn*ent2.inertia_inv*k
        
        if ent1.active == 0:
            if abs(ent1.speed[0]) < 1 and abs(ent1.speed[1]) < 1 and abs(ent1.speed[2]) < 0.1:
                ent1.speed = [0,0,0]
                ent1.pos[0] = ent1.active_state[0][0]
                ent1.pos[1] = ent1.active_state[0][1]
                ent1.pos[2] = ent1.active_state[0][2]
            else:
                if ent2.id_num == 0:
                    print 'negated2'
                ent1.active = 1
        if ent2.active == 0:
            if abs(ent2.speed[0]) < 1 and abs(ent2.speed[1]) < 1 and abs(ent2.speed[2]) < 0.1:
                ent2.speed = [0,0,0]
                ent2.pos[0] = ent2.active_state[0][0]
                ent2.pos[1] = ent2.active_state[0][1]
                ent2.pos[2] = ent2.active_state[0][2]
            else:
                if ent2.id_num == 0:
                    print 'negated3', ent2.speed
                ent2.active = 1
            
        
        v_1p = (ent1.speed[0] - ent1.speed[2]*(poc[1] - ent1.pos[1]), ent1.speed[1] + ent1.speed[2]*(poc[0] - ent1.pos[0]))
        v_2p = (ent2.speed[0] - ent2.speed[2]*(poc[1] - ent2.pos[1]), ent2.speed[1] + ent2.speed[2]*(poc[0] - ent2.pos[0]))


class Gravity():
    def __init__(self, accel, mass_inv):
        self.f = (accel[0]/mass_inv, accel[1]/mass_inv)
        
    def calculate(self, x, v, t):
        return ((self.f[0], self.f[1]), 0)

world = World()


for i in range(8):
    hexagon = uniform_polygon([600, 600 - i*dis*2, 0], 6, 50, "hexagon50v3.png")
    hexagon.speed = [0, 0, 0]#1.6]
    world.add_entity(hexagon)
#    world.add_sysobj(hexagon)
    del hexagon
for i in range(0):
    hexagon = uniform_polygon([500, 600 - i*dis, 0], 6, 25, "hexagon25.png")
    hexagon.speed = [0, 0, 0]#1.6]
    world.add_entity(hexagon)
#    world.add_sysobj(hexagon)
    del hexagon
for i in range(0):
    hexagon = uniform_polygon([700, 600 - i*dis, 0], 6, 25, "hexagon25.png")
    hexagon.speed = [0, 0, 0]#1.6]
    world.add_entity(hexagon)
#    world.add_sysobj(hexagon)
    del hexagon
for i in range(0):
    square = uniform_polygon([500, 300 + i*dis, 0], 4, 20, "box20.png")
    square.speed = [0, 0, 0]#1.6]
    world.add_entity(square)
#    world.add_sysobj(square)
    del square
for i in range(0):
    square = uniform_polygon([700, 300 + i*dis, 0], 4, 20, "box20.png")
    square.speed = [0, 0, 0]#1.6]
    world.add_entity(square)
#    world.add_sysobj(square)
    del square

p = Player()
p.mass_inv /= 10.0
p.inertia_inv /= 10.0
world.add_entity(p)
del p


#world.add_sysobj(Smoke_gen((300, 500, 10, 10)))
#world.add_sysobj(Smoke_gen(500, 100, 400, 100))


#entity_list = [hexagon, p1]


tiling0 = []

for i in range(map_height):
    tiling0.append([])
    for j in range(map_width):
        tiling0[i].append([])

tiling = copy_list3(tiling0)

ground = Polygon([275,600,0], ([-175,-51], [175,49], [175,51], [-175,51]), None)
ground.mass_inv = 0
ground.inertia_inv = 0
world.add_entity(ground)
del ground
ground = Polygon([626,649,0], ([-175,-2], [175,-2], [175,2], [-175,2]), None)
ground.mass_inv = 0
ground.inertia_inv = 0
world.add_entity(ground)
del ground
ground = Polygon([804,300,0], ([-2,-350], [2,-338], [2,338], [-2,350]), None)
ground.mass_inv = 0
ground.inertia_inv = 0
world.add_entity(ground)
del ground

#world = World()
#world.new_map()

#world.scr = screen_1
#del screen_1

#images that is accessed hundreds of time in a frame.
#locked_image = []


#################################
world.current_t = time()

#class State():
#    def __init__(self,):
        

while True:
    #sets max fps
    world.frame_t = world.clock.tick(fps)/1000.0
    
    if world.frame_t > 0.25:
        world.frame_t = 0.25

    world.remainder_t += world.frame_t

    while world.remainder_t >= world.dt: #the codes inside this while loop NEEDS to be executed in a time significantly LESS than dt or game slows down!
        #event handler
        event1 = pygame.event.poll()
        
        #system events
        if event1.type == KEYDOWN and event1.key == K_z:
            world.add_sysobj(Smoke_gen((100, 500, 10, 10), 10))
        if event1.type == QUIT: 
            exit()
        if event1.type == KEYDOWN and event1.key == K_TAB:
            fullscreen = not fullscreen
            if fullscreen: screen1 = pygame.display.set_mode((cam_width, cam_height), FULLSCREEN, 32) # DOUBLEBUF | HWSURFACE | 
            else: screen1 = pygame.display.set_mode((cam_width, cam_height), 0, 32)

        #the pysics processes.
        world.t += world.dt
        world.remainder_t -= world.dt

        tiling = copy_list3(tiling0)

        #prevState = currentState
        world.process(event1) #update currentState

#remember to handle temporal aliasing:
# NOTE that it's not as simple as positional coordinates. animation needs to be based on TIME instead of FRAMES from now on, though i guess i alr use TIME.
#      but my point is, for example let the walk animation consist of 9 frames indicated by the parameter "walk" = 0,1,2,…8, with each frame supposed to be 0.1s.
#      when doing "currentState*alpha + prevState*(1-alpha)", the "walk" also needs to be recalculated. maybe define the method .walk for walk animation, and
#      do the same for all other possible animations, and store the time elapsed as well as the time per frame in a dictionary. but for now let's not talk about
#      animation and do rigid body dynamics
# OH also consider what happens when an object is created/deleted. i guess for creation there is no need for intrapolation since by right it is created only 
# at the currentState and the renderer is rendering state at some time BEFORE currentState. the same for deletion.
# AND consider running the rendering on different threads though i dont see why is this necessary?


#alpha = world.remainder_t/world.dt
#state = currentState*alpha + prevState*(1-alpha)
#
#world.render(state)

#the display update
    world.render()
    pygame.display.update()

    print world.frame_t

    