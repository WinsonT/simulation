from sys import path
#path.append(path[0]+"/Maps")  
#path.append(path[0]+"/NPC")  

from WorldTemplate import WorldTemplate
from Polygon import Polygon, regularPolygon

#CONSTANTS
GRAVITY = 2400

#World1 is a place for performing various tests like rigid body simulation.
class World1(WorldTemplate): #WorldTemplate is the parent class containing the properties that all worlds have .entityList, .dt, .clock, etc.
# 
#                           |
#                           |
#                           |
#                           |
#                           |
#                           |
# |\______                  | 
# |       \___              |
# |___________\             |
#             =============

    def __init__(self, screen):
        WorldTemplate.__init__(self, screen)
        
        self.gravity = GRAVITY
        
        side = 50
        stack = 5
        
        #some hexagon
        self.addEntity(Polygon(self, [300,600 + (50-side),0], regularPolygon(6, side), "hexagon50.png", (50,43), [0, 0, 8]))
        self.addEntity(Polygon(self, [593,600 + (50-side),0], regularPolygon(6, side), "hexagon50.png", (50,43), [0, 0, 0]))
        self.addEntity(Polygon(self, [400,647,0], ((-400,-3),(400,-3),(400,3),(-400,3)), "hexagon50.png", (50,43), [0, 0, 0]))
        self.addEntity(Polygon(self, [700,635,0], ((-3,-3),(3,-3),(3,3),(3,3)), "hexagon50.png", (50,43), [0, 0, 0]))
        self.addEntity(Polygon(self, [100,635,0], ((-3,-3),(3,-3),(3,3),(3,3)), "hexagon50.png", (50,43), [0, 0, 0]))
        for i in xrange(stack-1):
            self.addEntity(Polygon(self, [593,600 + (50-side) - (i+1)*(side/5+side*(3)**0.5),0], regularPolygon(6, side+i*10), "hexagon50.png", (50,43), [0, 0, 0]))        
#        for i in xrange(stack-1):
#            self.addEntity(Polygon(self, [300,600 + (50-side) - (i+1)*(side/5+side*(3)**0.5),0], regularPolygon(6, side), "hexagon50.png", (50,43), [0, 0, 0]))        
        
        self.entityList[3].massInv = 0
        self.entityList[3].inertiaInv = 0
        self.entityList[4].massInv = 0
        self.entityList[4].inertiaInv = 0
        self.entityList[5].massInv = 0
        self.entityList[5].inertiaInv = 0

#        self.entityList[2].inertiaInv = 0
#        self.entityList[1].inertiaInv = 0
    def process(self, event):
        WorldTemplate.preProcess(self, event)
        
        '''force and constraint update'''
        #add gravity
        
        #constraint between body 1 and 2
#        F = (self.entityList[2].pos[0]-self.entityList[1].pos[0],self.entityList[2].pos[1]-self.entityList[1].pos[1],0,\
#             self.entityList[1].pos[0]-self.entityList[2].pos[0],self.entityList[1].pos[1]-self.entityList[2].pos[1],0)
#        self.entityList[1].speed[0] += F[0]*self.dt/10
#        self.entityList[1].speed[1] += F[1]*self.dt/10
#        self.entityList[2].speed[0] += F[3]*self.dt/10
#        self.entityList[2].speed[1] += F[4]*self.dt/10
        for i in self.entityList:
            if self.entityList[i].massInv != 0: self.entityList[i].speed[1] += self.gravity*self.dt
        
        
        WorldTemplate.postProcess(self, event)

    def render(self):
        self.screen.fill((155, 155, 155))
        WorldTemplate.render(self)
    
