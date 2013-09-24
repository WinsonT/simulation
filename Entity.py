from sys import path
#path.append(path[0]+"/Maps")  
#path.append(path[0]+"/NPC")  

import pygame
import os
from pygame.locals import *
from sys import exit
from time import time
from math import pi, sin, cos

pygame.init()

#Entity itself has no image, no brain, and will not collide. It serves as a template for 
class Entity():
    def __init__(self, world, pos):
        self.world = world
        
#        self.health = 1 #if .health < 0, entity is removed
    
        self.pos = pos 
        self.speed = [0,0,0] #
        
    def process(self, event): #subclasses that inherit Entity can fill processes with the Task class
        if self.health <= 0:
            self.death()    #in anticipation of death animation. if there is one, death() can be redefined in the corresponding Classes
    #(note that the Entity.death() is actually still there. see pygame ebook pg 148 GameEntity.render() event though Ant has .render())



    def damaged(self, dmg = 0):
        self.health -= dmg
        
    def death():
        self.world.remove_entity(self.id_num)

    def render():
        pass

    def trans(self): #semi-implicit Euler, using speed at t+dt to integrate position
        self.pos = [self.pos[0]+self.speed[0]*self.world.dt, self.pos[1]+self.speed[1]*self.world.dt, self.pos[2]+self.speed[2]*self.world.dt]

        #-2*pi <= .pos[2] <= 2*pi. the over stretched limit is to prevent pos[2] jumping from 0.0 to 6.283185307179, disrupting the deactivation algorithm
#        if self.pos[2] >= 2*pi:
#            self.pos[2] -= 2*pi
#        elif self.pos[2] <= -2*pi:
#            self.pos[2] += 2*pi
            
#    def update_img_rotation(self, angle):
#        self.image = pygame.transform.rotate(self.image0, -180*angle/pi) #pygame rotation is counter-clockwise

