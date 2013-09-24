from sys import path
#path.append(path[0]+"/Maps")  
#path.append(path[0]+"/NPC")  

import pygame
import os
from pygame.locals import*
from sys import exit
from time import time
from random import random, randint, randrange

import World1

pygame.init()

#CONSTANTS
FPS = 50
CAM_WIDTH = 720
CAM_HEIGHT = 960

screen_1 = pygame.display.set_mode((CAM_HEIGHT, CAM_WIDTH), 0, 32)
pygame.display.set_caption("Simultaneous Impulse!")
fullscreen = False
#pygame.event.set_blocked((MOUSEMOTION))
#AFRAMEPASSED = USEREVENT + 1


def main():
    #initialize tick
    world = World1.World1(screen_1) #load World1, a world for rigid body simulation
    
    world.clock.tick(FPS) #tick before entering main loop
    while True:
        #limits max fps
        world.frame_t = world.clock.tick(FPS)/1000.0
    
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
                return 0
            if event1.type == KEYDOWN and event1.key == K_TAB:
                fullscreen = not fullscreen
                if fullscreen: screen1 = pygame.display.set_mode((cam_width, cam_height), FULLSCREEN, 32) # DOUBLEBUF | HWSURFACE | 
                else: screen1 = pygame.display.set_mode((cam_width, cam_height), 0, 32)

            #the pysics processes.
            world.t += world.dt
            world.remainder_t -= world.dt

#            tiling = copy_list3(tiling0)

            #prevState = currentState
            world.process(event1) #update currentState

#remember to handle temporal aliasing:
# NOTE that it's not as simple as positional coordinates. animation needs to be based on TIME instead of FRAMES from now on, though i guess i alr use TIME.
#      but my point is, for example let the walk animation consist of 9 frames indicated by the parameter "walk" = 0,1,2,8, with each frame supposed to be 0.1s.
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

#        print world.frame_t
    

if __name__ == "__main__":
    exit(main())    