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

    def process(self, event):
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
            self.prev_pos = (self.prev_pos[0], self.prev_pos[1], self.prev_pos[2]-2*pi) #to prevent stutter due to ANTI-(temporal)aliasing which requires prevpos
        elif self.pos[2] <= -2*pi:
            self.pos[2] += 2*pi
            self.prev_pos = (self.prev_pos[0], self.prev_pos[1], self.prev_pos[2]+2*pi) #to prevent stutter due to ANTI-(temporal)aliasing which requires prevpos
            
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

    def trans(self): #Entity's trans doesn't need integration (for now), as an Entity moves with a constant speed

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
        if theta < 0: #0 < theta < 2*pi
            theta += 2*pi
        if theta > pi: #rotation by pi doesnt change width/height, so theta and theta-pi are equivalent
            theta -= pi
        if theta > pi/2: #finally rotation by theta = pi/2 + alpha is equivalent to pi/2-alpha = pi - theta. effectively, in the end 0 < theta < pi/2
            theta = pi - theta

        
        #algorith needed to be able to have COM that is not in the middle. this can be implemented simply by finding the line that connects COM and the "middle point"
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

    if inherit != 0: Polygon.__init__(inherit, pos, vertices0, image)
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

