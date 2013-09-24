from sys import path
#path.append(path[0]+"/Maps")  
#path.append(path[0]+"/NPC")  

import pygame
from time import time

#COMMENTS
#So quite OK now. Only how to stabilize stack even more? Currently a stack of only 3 hexagons is shaking (gravity of 2400). 1 way is to increase PGS iteration
#But couldnt contact caching take care of this??
#I think i know. After trying position-based error correction, i found that objects dont penetrate even at ERP == 0 (although they jitter if RESTITUTION > 0),
#so that the error-correction is totally useless. This is because the velocity constraint alone (given current model of collision detection, especially the two
#points of collision part) is able to keep objects from penetrating when SKIN is large enough. This brings us to the core of the problem. This means the role of
#Baumgarte's stabilization here is to bring 2 objects CLOSER together, stabilizing at close-to-0 separation. But the object never is stable (separation varies
#a little). So the Baumgarte is in fact what makes the stack unstable. All the more so if restitution > 0. Solution must be either in collision detection or
# contact caching. Or modifying parameters? also note that for large iterations stack is stable and separation goes to 0 (^-12 or so). And that larger coeff
# of friction makes the stack swings more.
#
#
#
#Next up, inclined plane test

#CONSTANTS

#FLOAT
CIRCLE_TO_CELL_RATIO = 1.0/4
GRAD_TOL = 0.001 #if too large, objects drift due to inaccuracy, especially if frictionless
THETA_TOL = 0.001
ERP = 0.1 #needs to be adjusted depending on restitution
RESTITUTION = 0. #needs to be tuned together with ERP
FRICTION_COEFF = 0.2
RELAXATION = 1. #0 < RELAXATION < 2, but 1.000 seems to be the best.

#INTEGERS
MIN_CELL_SIZE = 64
CELL_TO_CELL_RATIO = 2
HGRID_MAX_LEVELS = 10#5
NUM_BUCKETS = 1024
SKIN = 3
PGS_ITERATION = 10
#ER_ITERATION = 10


#Two kinds of game entity: Tangible and Decoration
#   (based on whether they interact/collide with others through touch)
#This is so that collision detection does not have to ask every time "who does not collide?"
#The division is not based on whether they have image because renderer does not have to ask "who does not have an image?"
#as rendering invokes a .render() method specific to each object. Those with no image simply have a .render() that does nothing
#
#Everything has: -.process()
#                -.trans()
#                -.render()
#                -.pos() (actually not yet used)
#
#Tangible: -interact with others through touch. tested in collision detection
#          -might have no image, in which case .render() does nothing
#          -broad phase (hgrid test) collision detection requires .BC. Also .nextObj, .timeStamp, but whether these belong to Tangible or BC 
#               hasnt been decided yet
#                   > also, object must be small enough to be contained in the largest hgrid
#          -For collision groups (members of each group being invisible to each other), before group 1 is inserted to hgrid, they are queried against ground.
#               then when inserting group1, no need to query for collision. next before group 2 is inserted, each is queried against objects already
#               inserted in the hgrid. Then group 2 is inserted without querying. then before group 3 is inserted, each is queried. etc.
#          -have "shape type": Polygon or Circle
#          -narrow phase collision detection for
#               Polygon: .vert, .normal, .axes.
#               Circle: .r
#
#Decoration: -does not interact with others via touch (hence is undetectable by game entities).
#            -only have a pre-programmed animation, whether a repeated one (e.g., flickering candle) or not (e.g, flickering candle that will extinguish)
#            -animation is specific per object. So there is no parent class Decoration. To make candles, create class Candle and define the animation specificly for it.
#
#
#
#

'''
#returns true if the minimum of aabb1 is < minimum of aabb2
def xIsLess(aabb1, aabb2):
    return aabb1.min[0] < aabb2.min[0]
def yIsLess(aabb1, aabb2):
    return aabb1.min[1] < aabb2.min[1]
#returns true if the radius of bc1 is < radius of bc2
def rIsGreater(bc1, bc2):
    return bc1.r > bc2.r

#a function used for quicksort. put all values "cmp"-er than list[pivotIdx] on the left side. left(right) is the left-(right-)most index of the list
def partition(list, left, right, pivotIdx, cmp):
    pivot = list[pivotIdx]
    list[pivotIdx], list[right] = list[right], pivot #move pivot to end
    storeIdx = left
    for i in xrange(left, right):
        if cmp(list[i], pivot):
            list[i], list[storeIdx] = list[storeIdx], list[i]
            storeIdx += 1
    list[storeIdx], list[right] = pivot, list[storeIdx] #move pivot to end
    return storeIdx

#quick sort. cmp is the comparator function. uses insertion sort whenever array size <= 10.
#left will be the cmp-est
def qsort(list, left, right, cmp):
    if right-left < 10:
        pivotIdx = (left + right)/2 #rounding off. NEVER choose left/right as this re-orders an already sorted list.
        
        #Get lists of bigger and smaller items and final position of pivot
        pivotNewIdx = partition(list, left, right, pivotIdx, cmp)
 
        #Recursively sort elements "cmp"-er than the pivot
        quicksort(list, left, pivotNewIdx - 1)
 
        #Recursively sort elements at least as "not cmp" as the pivot
        quicksort(list, pivotNewIdx + 1, right)
    else:
        isort(list, left, right, cmp)
        
#insertion sort. left will be the cmp-est
def isort(list, left, right, cmp):
    #start with the second element
    for i in xrange(left + 1, right):
        pivot = list[i]
        pivotTestIdx = i
        
        #keep moving test index to the left until either list[pivotTestIdx-1] is "cmp"-er than pivot or pivotTestIdx == left
        while (pivotTestIdx > left and cmp(pivot, list[pivotTestIdx - 1])): #SHORT CIRCUIT: second expression invalid for pivotTestIdx == 0
            #and while doing so, replace list[pivotTestIdx] with the element to its left
            list[pivotTestIdx] = list[pivotTestIdx-1]
            pivotTestIdx -= 1
            
        #put the pivot into place
        list[pivotTestIdx] = pivot
'''
     
def ComputeHashBucketIdx(x, y, z):
    h1 = int(0x8da6b343%NUM_BUCKETS) # Large multiplicative constants
    h2 = int(0xd8163841%NUM_BUCKETS) # here arbitrarily chosen primes
    h3 = int(0xcb1ab31f%NUM_BUCKETS) #
    
    n = (h1 * x + h2 * y + h3 * z) % NUM_BUCKETS
    if (n < 0): n += NUM_BUCKETS
    return n
        
            
            
class Hgrid():
    def __init__(self):
        #Time stamp (timeStamp must also be reset when tick is reset)
        self.tick = 0
        self.timeStamp = [0 for i in xrange(NUM_BUCKETS)]
        
        self.occupiedLevelsMask = 0
        self.objectBucket = [None for i in xrange(NUM_BUCKETS)]
        
    def resetHgrid(self):
        self.tick = 0
        self.timeStamp = [0 for i in xrange(NUM_BUCKETS)]
        self.occupiedLevelsMask = 0
        self.objectBucket = [None for i in xrange(NUM_BUCKETS)]
    
    

#RULES
#-World calls addEntity when adding an object. Tangible calls addBC in its initialization
#-Object calls removeBC (if relevant) and removeEntity upon death
#
class WorldTemplate(): #WorldTemplate is the one containing the properties that all worlds have, like entity list, time, timpe_step, clock, etc.
    def __init__(self, screen):
        #This is the main screen, the window screen, where entities will blit themselves.
        self.screen = screen #there are more screens for every World. There might be BG1, BG2, etc.

        self.clock = pygame.time.Clock()
        self.frame_t = 0 #the time elapsed since the last frame
        self.t = 0 #the time up to which the game has simulated. starts with 0
        self.dt = 1.0/60 #time step for calculations, a constant
        self.remainder_t = 0

        #create dictionaries
        #world is responsible for adding entities into this dictionary
        self.entityList = {} #entityList contains everything, objects that are Tangible and are not, that have images and do not.
        self.next_idNum = 1 #reserve idNum == 0 for immovable? useful for constraint solver
        
        self.BC_list = [] #sorted list in decreasing radius of every BC of each Tangible. Tangible is responsible for adding their BC into this list
        self.numberOfBC = 0
        
        self.grid = Hgrid()


        #constraints are re-registered into solver every frame
        #but world can keep a list of (non-collision) Constraint class (e.g, Joints), which are kept (instead of be reset every frame).
        #each constraint has a method (.register, or .process) which are called each frame to re-register into solver. some checks can also be made. e.g.,
        #a check to see if a destructor should be called to remove oneself from world's Constraint list. Also, the lambda values for these constraint are saved to be
        #used by the solver as initial guess.
        
        #FORMULATION
        #A constraint on 2 bodies is of the form J V = zeta, V is the velocity of the two bodies. J is the jacobian = derivative of the position constraint Cn.
        #zeta is the desired value of J V. for example, for collision constraint, zeta = elasticity*initial relative speed - beta*Cn, the first term corresponds to
        #the elasticity, the second Baumgarte stabilization, with beta = ERP/dt, ERP < 1 (experimentally determined). A motor can also be formulated as constraint
        #with zeta such that it does the desired action. Jmap is to store the sparse J into a dense form (see paper). The force on the bodies will be J^T lambda0,
        #hence the aim of the solver is to determine lambda0. The minimum and maximum allowed values of lambda for each constraint is stored in lambdaMinMax.
        #contactIDList is used to cache the identifier of current collisions for querying if contact is previously present. If yes, the previous lambda value is used
        #as initial guess. The stored lambda values itself is in storedLambdaValues. The number of constraints present is counted in numOfConstraint.

        #COLLISION SOLVER VARIABLES. see IterativeDynamics.pdf by Erin Catto
        self.zeta = [] #an s times 1 column vector
        self.Jmap = []
        self.J = []
        self.lambdaMinMax = []
        self.contactIDList = []
        self.numOfConstraint = 0
                
        #STORED LAMBDA VALUES
        #For collision constraint, key is of the form (obj1.idNum, obj2.idNum, obj1collisionVertex, obj2collisionVertex), obj1.idNum being smaller than obj2's
        #For joints, ... (to be announced)
        self.storedLambdaValues = {}

        #SO FAR NOT GOOD
        #objects don't sink into another even at ERP 0, so this position based error correction does nothing at all. This is bcs SKIN == 3, and 
        #with 2 points of collisions the velocity constraints alone managed to prevent bodies from penetrating. If SKIN == 0, objects jitter bcs
        #now the velocity constraint is not good enough to handle constraint. So it's better to go back to Baumgarte's stabilization. Keep in mind though,
        #this means that the role of the ERP is actually to get the objects CLOSER together if SKIN >= 3.
#        #POSITION BASED ERROR CORRECTION - SIMULTANEOUS METHOD
#        self.zeta2 = [] #an s times 1 column vector
#        self.Jmap2 = []
#        self.J2 = []
#        self.kappaMinMax = []
#        self.contactIDList2 = []
#        self.numOfConstraint2 = 0
#        self.storedKappaValues = {}
        
		#create a list of group and contact points? for sleeping. but later.
		

    def addEntity(self, entity):
        self.entityList[self.next_idNum] = entity
        entity.idNum = self.next_idNum
        self.next_idNum += 1

    def removeEntity(self, idNum):
        del self.entityList[idNum]


    def addBC(self, BC):
        self.BC_list.append(BC)
        
        #insertion sort only the last element
        idx = self.numberOfBC
        #keep moving index to the left until either list[idx-1] is not smaller than BC or idx == 0
        while (idx > 0 and self.BC_list[idx - 1].r < BC.r): #SHORT CIRCUIT: second expression invalid for pivotTestIdx == 0
            #and while doing so, replace list[idx] with the element to its left
            self.BC_list[idx] = self.BC_list[idx-1]
            idx -= 1
        #put the pivot into place
        self.BC_list[idx] = BC

        self.numberOfBC += 1

    def removeBC(self, object):
        self.numberOfBC -= 1
        
        #index to rightmost element
        right = self.numberOfBC
        #next to be checked
        idx = right - 1
        #while the rightmost is not the object
        while (self.BC_list[right].object.idNum != object.idNum): #object must be inside or error occurred
                #swap the next to be checked with the rightmost
                self.BC_list[right], self.BC_list[idx] = self.BC_list[idx], self.BC_list[right]
                idx -= 1
        
        del self.BC_list[right]


    def preProcess(self, event):
        self.grid.resetHgrid()
        self.zeta = []
        self.Jmap = []
        self.J = []
        self.lambdaMinMax = []
        self.contactIDList = []
        self.numOfConstraint = 0
#        self.zeta2 = []
#        self.Jmap2 = []
#        self.J2 = []
#        self.kappaMinMax = []
#        self.contactIDList2 = []
#        self.numOfConstraint2 = 0

    def postProcess(self, event):
        timestart = time()
        
        '''force and constraint update (cont'd)'''
        #pre to calling WorldTemplate.postProcess(), world.process() should have applied GENERAL external forces such as gravity,
        #and registered any constraint present in the system.
        
        #collision detection, includes registering collision constraints, and perhaps external forces by terrain on entities
        for circle in self.BC_list:
            self.registerToHGrid(circle.object)

        '''velocity update (and constraints from characters)'''
        #entities process. entities should NOT modify position.
        #they only propose velocities, which include velocity change computation due to external forces
        #AND perhaps additional constraints
        for entity in self.entityList.values():
            entity.process(event)
        
        '''velocity check and error correction'''
        #simultaneous-impulse based constraint solver
        self.solveConstraint()
                
        '''position update'''
        #modify position, using the velocities approved by constraint handler
        for entity in self.entityList.values():
            entity.trans()
            
        print "time elapsed in calculation", time() - timestart
            
    def render(self):
        #screen reset should be handled by each world. worldTemplate only calls render for every entity
    
        #parameter for anti-(temporal) aliasing
        alpha = self.remainder_t/self.dt
        
        for entity in self.entityList.values():
            entity.render(self.screen, alpha)
            

    '''    
    def sortSweep(self): #prefers AABB to have min and max
        n = self.numberOfAABB
        
        #x-axis test. objects are assumed to have greater variance in the x direction
        qsort(self.AABB_list, 0, n-1, xIsLess)
        
        #sweep forward
        for i in xrange(0,n):
            for j in xrange(i+1,n):
            #stop if AABB j is outside AABB i
                if self.AABB_list[j].min[0] > self.AABB_list[i].max[0]:
                    break
                else:
                    #check for y overlap
                    if self.AABB_list[i].min[1] > self.AABB_list[j].max[1] or self.AABB_list[i].max[1] < self.AABB_list[j].min[1]:
                        #check pair
                        collDetect(self.AABB_list[i], self.AABB_list[j])
                        
        #This is OK for object against object collision test, but not for collision with ground/tiles.
    '''


    #TO DO: handle collision group (see notes on Tangible above)
    #add object to hgrid, at the same time checking for collision
    def registerToHGrid(self, obj):
        grid = self.grid
        
        #for time stamp
        grid.tick += 1
        
        '''Find lowest level where object fully fits inside cell, taking RATIO into account'''
        startLevel = 0
        cellSize = MIN_CELL_SIZE
        r = obj.BC.r
        diameter = 2*r
        
        #CIRCLE_TO_CELL_RATIO < 1, hence max overlap is 4 tiles
        while (cellSize * CIRCLE_TO_CELL_RATIO < diameter):
            cellSize *= CELL_TO_CELL_RATIO
            startLevel += 1

        #Assert if object is larger than largest grid cell
        assert(startLevel < HGRID_MAX_LEVELS)
        
        # Find all cells overlapped by object, extended by the max radius can be contained in the level. and SKIN
        delta = cellSize * CIRCLE_TO_CELL_RATIO/2.0 + r + SKIN
        cx = obj.BC.c[0]
        cy = obj.BC.c[1]
        x1 = int((cx-delta)/ cellSize) #note rounding
        y1 = int((cy-delta)/ cellSize)
        x2 = int((cx+delta)/ cellSize)
        y2 = int((cy+delta)/ cellSize)
                
        #and for each of the cell 
        for i in xrange(x1,x2+1):
            for j in xrange(y1,y2+1):
                #starting with this level
                mask = grid.occupiedLevelsMask >> startLevel
                
                '''Check collision'''
                level = startLevel
                cellX = i
                cellY = j
                # with objects located in this cell or above
                while (level < HGRID_MAX_LEVELS):
                    #if no more object
                    if mask == 0:
                        break
                    #if level is empty
                    if mask & 1 == 0:
                        cellX = cellX / CELL_TO_CELL_RATIO #note rounding
                        cellY = cellY / CELL_TO_CELL_RATIO
                        mask >>= 1
                        level += 1
                        continue
                        
                    #compute cell hash ID
                    bucket = ComputeHashBucketIdx(cellX, cellY, level)
                                        
                    #if cell has not been checked (essential as higher level cells maybe checked several times)
                    if grid.timeStamp[bucket] != grid.tick:
                        grid.timeStamp[bucket] = grid.tick
                                                
                        #while there is still object to check
                        obj2 = grid.objectBucket[bucket]
                        
                        
                        while obj2 != None:
                            #if obj2 is truly nearby
                            if (obj2.BC.c[0]-cx)**2 + (obj2.BC.c[1]-cy)**2 < (obj2.BC.r + r + SKIN)**2:
                                #test obj against obj2. also, if one is immovable, set this to be obj (this is the convention required)
                                self.separatingAxisTest(obj, obj2)
                                minOverlap = obj.BC.r + obj2.BC.r + 4*SKIN - ((obj.BC.c[0] - obj2.BC.c[0])**2 + (obj.BC.c[1] - obj2.BC.c[1])**2)**0.5
                            obj2 = obj2.nextObj

                    cellX = cellX / CELL_TO_CELL_RATIO #note rounding
                    cellY = cellY / CELL_TO_CELL_RATIO
                    mask >>= 1
                    level += 1
                    
        '''and add object to head of linked list'''
        bucket = ComputeHashBucketIdx(int(cx/cellSize), int(cy/cellSize), startLevel)
        obj.nextObj = grid.objectBucket[bucket]
        grid.objectBucket[bucket] = obj
        
        # Mark this level as having one more object. Also indicate level is in use
        grid.occupiedLevelsMask |= (1 << startLevel)
                
   
        
    def separatingAxisTest(self, obj1, obj2):
#        timestart = time()
        
        if obj1.massInv == 0 and obj2.massInv == 0: return 0
        
        '''work in world coordinate'''
        #instead of local ones. this is so that when doing subsequent tests for obj1/obj2, the vertices and normals are already correctly rotated 
        if abs(obj1.lastVertTheta - obj1.pos[2]) > THETA_TOL: obj1.updateRotation()
        if abs(obj2.lastVertTheta - obj2.pos[2]) > THETA_TOL: obj2.updateRotation()
                
        #code is faster if obj1 have more normals??????????????????????????????????????????????????
#        if obj1.nn < obj2.nn:
#            obj1, obj2 = obj2, obj1

        if obj1.idNum > obj2.idNum:
            obj1, obj2 = obj2, obj1
        
        #the vector connecting the center of mass of obj1 and obj2
        obj12X = obj2.pos[0] - obj1.pos[0]
        obj12Y = obj2.pos[1] - obj1.pos[1]

        '''initialize output parameters'''
        minOverlap = obj1.BC.r + obj2.BC.r + 4*SKIN - ((obj1.BC.c[0] - obj2.BC.c[0])**2 + (obj1.BC.c[1] - obj2.BC.c[1])**2)**0.5
        #pushing obj2 along minAxis resolve collision
        minAxis = (obj12X, obj12Y)
        #Collision point(s) in world coordinate (in terms of orientation), using respective center of mass as origin. Note there can be up to two collision points
        obj1cp1 = None 
        obj2cp1 = None
        obj1cp2 = None
        obj2cp2 = None
        obj1collidingVertex1 = None
        obj2collidingVertex1 = None
        obj1collidingVertex2 = None
        obj2collidingVertex2 = None
        
        
        #QUESTION TO PONDER
        #Wait a minute! minAxis will NOT differ from the initialized vector by more than 90 degree, provided that the centers of mass are contained within
        #each polygon AND that collision is detected before any 'significant' overlap is made! less robust (due to the last requirement), but much simpler!
        
        
        #if bit is 1, means axis is a duplicate. applies only for normals of obj2
        duplicateAxisMask = 0

        '''for each axis of obj1'''
        for obj1AxisID in xrange(obj1.nn):
            #this will be minAxis if collision is at obj1 MAX. if at obj1 MIN, minAxis is (-nx,-ny)
            nx = obj1.normal[obj1AxisID][0]
            ny = obj1.normal[obj1AxisID][1]
                        
            #projection of the vector connecting the two bodies
            obj12N = obj12X*nx + obj12Y*ny

            
            #projection of first vertex of obj2
            k = 0
            minVal = obj2.vert[k][0]*nx + obj2.vert[k][1]*ny
            maxVal = minVal
            
            #2nd vertex
            k = 1
            value = obj2.vert[k][0]*nx + obj2.vert[k][1]*ny
            
            deltaMax = obj12N-obj1.projections[obj1AxisID].min #when added with the MAX projection of obj2, gives overlap
            deltaMin = -obj12N+obj1.projections[obj1AxisID].max #when subtracted by the MIN projection of obj2, gives overlap
            
            '''Case 1, if the same, (the lucky case where things can be resolved quickly)'''
            if abs(value - maxVal) < GRAD_TOL:
                #mark the axis in obj2 as duplicate
                duplicateAxisMask |= 1
            
                k = 2
                value = obj2.vert[k][0]*nx + obj2.vert[k][1]*ny
                
                #if k == 0 and 1 are MIN (obj2.normal is anti-parallel to (nx,ny))
                if value > minVal:
                    k = 0
                    #code is exactly like JUMPING TO *) BELOW, with k == 0, step == 1 and minVal == minVal2################################################
                    overlap = deltaMin - minVal

                    #if non overlapping, stop now
                    if SKIN + overlap < 0: 
#                        print 1, overlap
                        return 0 #-----------------------------------RETURN----------------
                    
                    projection = obj2.projections[k]
                    
                    #EXPLOITING INFORMATION FROM obj2.projections
                    #find second overlap, obj2's MIN with obj1's MIN
                    overlap2 = deltaMax - projection.min #== deltaMax + maxVal

                    #if non overlapping, stop now
                    if SKIN + overlap2 < 0:
#                        print 2, overlap
                        return 0 #-----------------------------------RETURN----------
                    
                    #if it's about obj1 MIN and obj2 MIN
                    if overlap2 < overlap and overlap2 < minOverlap:
                        minOverlap = overlap2
                        #NOTE REDEFINITION OF NORMAL
                        nx = -nx
                        ny = -ny
                        minAxis = (nx, ny)
                        
                        #if there is only one obj2 MIN
                        if projection.min2 == None:
                            #if there are two obj1 MINs
                            if obj1.projections[obj1AxisID].min2 != None:
                                #contact point is obj2 MIN
                                obj2cp1 = obj2.vert[projection.min1]
                                obj1cp1 = (obj2cp1[0]+nx*minOverlap+obj12X, obj2cp1[1]+ny*minOverlap+obj12Y)
                                obj1cp2 = None
                                obj2cp2 = None
                                obj1collidingVertex1 = None
                                obj2collidingVertex1 = projection.min1
                                obj1collidingVertex2 = None
                                obj2collidingVertex2 = None
                                
                            #else, do nothing since this surely isn't the separating axis. 
                            #we're done with this axis
                            continue
                            
                        #else, if there is only one obj1 MIN
                        elif obj1.projections[obj1AxisID].min2 == None:
                            #contact point is obj1 MIN
                            obj1cp1 = obj1.vert[obj1.projections[obj1AxisID].min1]
                            obj2cp1 = (obj1cp1[0]-nx*minOverlap-obj12X, obj1cp1[1]-ny*minOverlap-obj12Y)
                            obj1cp2 = None
                            obj2cp2 = None
                            obj1collidingVertex1 = obj1.projections[obj1AxisID].min1
                            obj2collidingVertex1 = None
                            obj1collidingVertex2 = None
                            obj2collidingVertex2 = None
                            
                            #we're done with this axis
                            continue
                            
                        else:
                            #the points in local coordinates
                            #projection of obj1's min2 on t is > obj1's min1, projection of obj2's min1 is > obj2's min2
                            #t for tangent == the minAxis rotated 90degree wrt z (CLOCKWISE).
                            obj1minID = obj1.projections[obj1AxisID].min1
                            obj1maxID = obj1.projections[obj1AxisID].min2
                            obj2minID = projection.min2
                            obj2maxID = projection.min1
                            
                            #... to 2 COLLISION POINTS below
                        

                    #else if it's about obj1 MAX and obj2 MAX
                    elif overlap < minOverlap:
                        minOverlap = overlap
                        minAxis = (nx, ny)
                        
                        #the points in local coordinates
                        #projection of obj1's max2 on t is > obj1's max1, projection of obj2's max1 is > obj2's max2,
                        #t is tangent == minAxis rotated 90degree wrt z (CLOCKWISE).
                        obj1minID = obj1.projections[obj1AxisID].max1
                        obj1maxID = obj1.projections[obj1AxisID].max2
                        obj2minID = projection.max2
                        obj2maxID = k

                        #... to 2 COLLISION POINTS below
                        
                        
                    #else minOverlap is smaller. don't update collision points, just mark axis as duplicate.
                    else:
                        #we're done with this axis
                        continue
                    
                    ###################################################################################################################################
                    
                    
                #else k == 0 and 1 are MAXes (obj2.normal is parallel to (nx,ny))
                else:                    
                    #overlap of obj1 MIN and obj2 MAX
                    overlap2 = deltaMax + maxVal
                    
                    #if non-overlapping, stop now
                    if SKIN + overlap2 < 0:
#                        print 3, overlap
                        return 0 #-----------------------------------RETURN----------

                    projection = obj2.projections[0]

                    #EXPLOITING INFORMATION FROM obj2.projections
                    #find overlap, obj2's MIN with obj1's MAX
                    overlap = deltaMin - projection.min
                
                    #if non overlapping, stop now
                    if SKIN + overlap < 0:
#                        print 4, overlap
                        return 0 #-----------------------------------RETURN----------

                    #if it's about obj1 MIN and obj2 MAX (COPIED-PASTED from below, from obj1 MIN and obj2 MAX section)
                    if overlap2 < overlap and overlap2 < minOverlap:
                        minOverlap = overlap2
                        #NOTE REDEFINITION OF NORMAL
                        nx = -nx
                        ny = -ny
                        minAxis = (nx, ny)

                        #if there is only one obj1 MIN, then that's the point
                        if obj1.projections[obj1AxisID].min2 == None:
                            obj1cp1 = obj1.vert[obj1.projections[obj1AxisID].min1]
                            obj2cp1 = (obj1cp1[0]-nx*minOverlap-obj12X, obj1cp1[1]-ny*minOverlap-obj12Y)
                            obj1cp2 = None
                            obj2cp2 = None
                            obj1collidingVertex1 = obj1.projections[obj1AxisID].min1
                            obj2collidingVertex1 = None
                            obj1collidingVertex2 = None
                            obj2collidingVertex2 = None
                            
                            #we're done with this axis
                            continue
                            
                        #else, two points of collision
                        #the points in local coordinates
                        #projection of obj1's min2 on t is > obj1's min1, projection of obj2's max1 is > obj2's max2
                        #t for tangent == the minAxis rotated 90degree wrt z (CLOCKWISE).
                        obj1minID = obj1.projections[obj1AxisID].min1
                        obj1maxID = obj1.projections[obj1AxisID].min2
                        obj2minID = projection.max2
                        obj2maxID = projection.max1

                            
                        #... to 2 COLLISION POINTS below
                    
                    #else if it's about obj1 MAX and obj2 MIN (COPIED-PASTED from below, from obj1 MAX and obj2 MIN section)
                    elif overlap < minOverlap:
                        minOverlap = overlap
                        minAxis = (nx, ny)
                        
                        #the points in local coordinates
                        #projection of obj1's max2 on t is > obj1's max1, projection of obj2's min1 is > obj2's min2,
                        #t is tangent == minAxis rotated 90degree wrt z (CLOCKWISE).
                        obj1minID = obj1.projections[obj1AxisID].max1
                        obj1maxID = obj1.projections[obj1AxisID].max2
                        obj2minID = projection.min2
                        obj2maxID = k
                                                        
                        #... to 2 COLLISION POINTS below
                                        
                    #else minOverlap is smaller. don't update collision points, just mark axis as duplicate.
                    else:
                        #we're done with this axis
                        continue

    
                #FOR 2 COLLISION POINTS
                #find and register the points
                obj1min = obj1.vert[obj1minID]
                obj1max = obj1.vert[obj1maxID]
                obj2min = obj2.vert[obj2minID]
                obj2max = obj2.vert[obj2maxID]
                                                                
                #t for tangent == minAxis rotated 90degree wrt z (CLOCKWISE).
                tx = -ny
                ty = nx
                obj12T = obj12X*tx + obj12Y*ty

                #obj1min projection < obj2min projection
                if obj1min[0]*tx + obj1min[1]*ty < obj2min[0]*tx + obj2min[1]*ty + obj12T:
                    #obj1max projection > obj2max projection
                    if obj1max[0]*tx + obj1max[1]*ty < obj2max[0]*tx + obj2max[1]*ty + obj12T:
                        #the two collision points are obj2's vertices
                        obj1cp1 = (obj2min[0]+nx*minOverlap+obj12X, obj2min[1]+ny*minOverlap+obj12Y)
                        obj1cp2 = (obj2max[0]+nx*minOverlap+obj12X, obj2max[1]+ny*minOverlap+obj12Y)
                        obj2cp1 = obj2min
                        obj2cp2 = obj2max
                        obj1collidingVertex1 = None
                        obj2collidingVertex1 = obj2minID
                        obj1collidingVertex2 = None
                        obj2collidingVertex2 = obj2maxID
                    else:
                        #the two collision points are obj2min and obj1max
                        obj1cp1 = (obj2min[0]+nx*minOverlap+obj12X,obj2min[1]+ny*minOverlap+obj12Y)
                        obj1cp2 = obj1max
                        obj2cp1 = obj2min
                        obj2cp2 = (obj1max[0]-nx*minOverlap-obj12X,obj1max[1]-ny*minOverlap-obj12Y)
                        obj1collidingVertex1 = None
                        obj2collidingVertex1 = obj2minID
                        obj1collidingVertex2 = obj1maxID
                        obj2collidingVertex2 = None
                else:
                    #obj1max projection > obj2max projection
                    if obj1max[0]*tx + obj1max[1]*ty < obj2max[0]*tx + obj2max[1]*ty + obj12T:
                        #the two collision points are obj1min and obj2max
                        obj1cp1 = obj1min
                        obj1cp2 = (obj2max[0]+nx*minOverlap+obj12X,obj2max[1]+ny*minOverlap+obj12Y)
                        obj2cp1 = (obj1min[0]-nx*minOverlap-obj12X,obj1min[1]-ny*minOverlap-obj12Y)
                        obj2cp2 = obj2max
                        obj1collidingVertex1 = obj1minID
                        obj2collidingVertex1 = None
                        obj1collidingVertex2 = None
                        obj2collidingVertex2 = obj2maxID
                    else:
                        #the two collision points are obj1's vertices
                        obj1cp1 = obj1min
                        obj1cp2 = obj1max
                        obj2cp1 = (obj1min[0]-nx*minOverlap-obj12X,obj1min[1]-ny*minOverlap-obj12Y)
                        obj2cp2 = (obj1max[0]-nx*minOverlap-obj12X, obj1max[1]-ny*minOverlap-obj12Y)
                        obj1collidingVertex1 = obj1minID
                        obj2collidingVertex1 = None
                        obj1collidingVertex2 = obj1maxID
                        obj2collidingVertex2 = None

                #we're done with this axis
                continue

                
            '''Case2, else'''
            #if to find MIN needs to traverse backward
            if value > maxVal:
                maxVal = value
                k = obj2.n
                step = -1
                k_final = 2
                #to find MAX if necessary
                j = 1
                j_final = obj2.n-1
            #else, traverse forward to find MIN
            else:
                minVal = value
                #k == 1
                step = 1
                k_final = obj2.n-1
                #to find MAX if necessary
                j = obj2.n
                j_final = 2
            
            #traverse to find MIN
            while (k != k_final):
                k += step #there is no need for k = (k+step+obj2.n)%obj2.n
                minVal2 = obj2.vert[k][0]*nx + obj2.vert[k][1]*ny
                if minVal2 > minVal - GRAD_TOL:
                    k = (k-step)%obj2.n #there is no need for (k-step+obj2.n)%obj2.n
                    break
                else:
                    minVal = minVal2
                    minVal2 += 1000*GRAD_TOL #sentinel

            #*) after minVal is found
            overlap = deltaMin - minVal
                        
            #if non overlapping, stop now
            if SKIN + overlap < 0:
#                print 5, overlap
                return 0 #-----------------------------------RETURN----------------

                  
            '''If there are 2 MIN vertices:'''
            # This section deals with the possibility that there are 2 contact points, 
            # while at the same time marking the corresponding axis of obj2 as "checked"
            # the result is that when checking obj2's axis, there can only be 1 contact point
            if minVal2 < minVal + GRAD_TOL:
                
                #k should be such that if one travels CLOCKWISE from it, the next vertex is the other MIN vertex
                #this is so that k is either a .max1 or .min1, but not .max2 or .min2
                if step == -1:
                    k = (k+step+obj2.n)%obj2.n
                
                projectionID = obj2.vertToNormal[k]
                projection = obj2.projections[projectionID]
                
                #mark the axis in obj2 as duplicate
                duplicateAxisMask |= 1 << projectionID

                #EXPLOITING INFORMATION FROM obj2.projections
                #if k is a MAX (note that here projection.max should be == -minVal, obj2.normal is antiparallel to (nx,ny))
                if projection.max1 == k:
                    #find second overlap, obj2's MIN with obj1's MIN
                    overlap2 = deltaMax - projection.min
                    
                    #if non overlapping, stop now
                    if SKIN + overlap2 < 0:
#                        print 6, overlap
                        return 0 #-----------------------------------RETURN----------
                    
                    #if it's about obj1 MIN and obj2 MIN
                    if overlap2 < overlap and overlap2 < minOverlap:
                        minOverlap = overlap2
                        #NOTE REDEFINITION OF NORMAL
                        nx = -nx
                        ny = -ny
                        minAxis = (nx, ny)
                                                
                        #if there is only one obj2 MIN
                        if projection.min2 == None:
                            #if there are two obj1 MINs
                            if obj1.projections[obj1AxisID].min2 != None:
                                #contact point is obj2 MIN
                                obj2cp1 = obj2.vert[projection.min1]
                                obj1cp1 = (obj2cp1[0]+nx*minOverlap+obj12X, obj2cp1[1]+ny*minOverlap+obj12Y)
                                obj1cp2 = None
                                obj2cp2 = None
                                obj1collidingVertex1 = None
                                obj2collidingVertex1 = projection.min1
                                obj1collidingVertex2 = None
                                obj2collidingVertex2 = None

                                
                            #else, do nothing since this surely isn't the separating axis. 
                            #we're done with this axis
                            continue
                            
                        #else, if there is only one obj1 MIN
                        elif obj1.projections[obj1AxisID].min2 == None:
                            #contact point is obj1 MIN
                            obj1cp1 = obj1.vert[obj1.projections[obj1AxisID].min1]
                            obj2cp1 = (obj1cp1[0]-nx*minOverlap-obj12X, obj1cp1[1]-ny*minOverlap-obj12Y)
                            obj1cp2 = None
                            obj2cp2 = None
                            obj1collidingVertex1 = obj1.projections[obj1AxisID].min1
                            obj2collidingVertex1 = None
                            obj1collidingVertex2 = None
                            obj2collidingVertex2 = None
                            
                            
                            #we're done with this axis
                            continue
                            
                        else:
                            #the points in local coordinates
                            #projection of obj1's min2 on t is > obj1's min1, projection of obj2's min1 is > obj2's min2
                            #t for tangent == the minAxis rotated 90degree wrt z (CLOCKWISE).
                            obj1minID = obj1.projections[obj1AxisID].min1
                            obj1maxID = obj1.projections[obj1AxisID].min2
                            obj2minID = projection.min2
                            obj2maxID = projection.min1
                            
                            #... to 2 COLLISION POINTS below
                        

                    #else if it's about obj1 MAX and obj2 MAX (so obj2 MAX (wrt obj2 normal) is MIN wrt obj1 normal as in this case they are antiparallel)
                    elif overlap < minOverlap:
                        minOverlap = overlap
                        minAxis = (nx, ny)
                        
                        #the points in local coordinates
                        #projection of obj1's max2 on t is > obj1's max1, projection of obj2's max1 is > obj2's max2,
                        #t is tangent == minAxis rotated 90degree wrt z (CLOCKWISE).
                        obj1minID = obj1.projections[obj1AxisID].max1
                        obj1maxID = obj1.projections[obj1AxisID].max2
                        obj2minID = projection.max2
                        obj2maxID = k

                        #... to 2 COLLISION POINTS below
                        
                        
                    #else minOverlap is smaller. don't update collision points, just mark axis as duplicate.
                    else:
                        #we're done with this axis
                        continue
                    
                    
                #else k is a MIN (note that here projection.min should be == minVal, obj2.normal is parallel to (nx,ny))
                #projection.min1 == k
                else:
                    #find second overlap, obj2's MAX with obj1's MIN
                    overlap2 = deltaMax + projection.max
                    
                    #if non overlapping, stop now
                    if SKIN + overlap2 < 0:
#                        print 7, overlap
                        return 0 #-----------------------------------RETURN----------
                    
                    #if it's about obj1 MIN and obj2 MAX
                    if overlap2 < overlap and overlap2 < minOverlap:
                        minOverlap = overlap2
                        #NOTE REDEFINITION OF NORMAL
                        nx = -nx
                        ny = -ny
                        minAxis = (nx, ny)
                        
                        #if obj1 has only 1 MIN, then this obj1 MIN is the point of collision
                        if obj1.projections[obj1AxisID].min2 == None:
                            obj1cp1 = obj1.vert[obj1.projections[obj1AxisID].min1]
                            obj2cp1 = (obj1cp1[0]-nx*minOverlap-obj12X, obj1cp1[1]-ny*minOverlap-obj12Y)
                            obj1cp2 = None
                            obj2cp2 = None
                            obj1collidingVertex1 = obj1.projections[obj1AxisID].min1
                            obj2collidingVertex1 = None
                            obj1collidingVertex2 = None
                            obj2collidingVertex2 = None
                            
                            #we're done with this axis
                            continue
                            
                        #else, two points of collision
                        #the points in local coordinates
                        #projection of obj1's min2 on t is > obj1's min1, projection of obj2's max1 is > obj2's max2
                        #t for tangent == the minAxis rotated 90degree wrt z (CLOCKWISE).
                        obj1minID = obj1.projections[obj1AxisID].min1
                        obj1maxID = obj1.projections[obj1AxisID].min2
                        obj2minID = projection.max2
                        obj2maxID = projection.max1
                            
                        #... to 2 COLLISION POINTS below
                                                                
                        
                    #else if it's about obj1 MAX and obj2 MIN
                    elif overlap < minOverlap:
                        minOverlap = overlap
                        minAxis = (nx, ny)

                        #the points in local coordinates
                        #projection of obj1's max2 on t is > obj1's max1, projection of obj2's min1 is > obj2's min2,
                        #t is tangent == minAxis rotated 90degree wrt z (CLOCKWISE).
                        obj1minID = obj1.projections[obj1AxisID].max1
                        obj1maxID = obj1.projections[obj1AxisID].max2
                        obj2minID = projection.min2
                        obj2maxID = k
                                                        
                        #... to 2 COLLISION POINTS below
                        
                        
                    #else minOverlap is smaller. don't update collision points, just mark axis as duplicate
                    else:
                        #we're done with this axis
                        continue

                
                #FOR 2 COLLISION POINTS
                #find and register the points
                obj1min = obj1.vert[obj1minID]
                obj1max = obj1.vert[obj1maxID]
                obj2min = obj2.vert[obj2minID]
                obj2max = obj2.vert[obj2maxID]
                
                #t for tangent == minAxis rotated 90degree wrt z (CLOCKWISE).
                tx = -ny
                ty = nx
                obj12T = obj12X*tx + obj12Y*ty

                #obj1min projection < obj2min projection
                if obj1min[0]*tx + obj1min[1]*ty < obj2min[0]*tx + obj2min[1]*ty + obj12T:
                    #obj1max projection > obj2max projection
                    if obj1max[0]*tx + obj1max[1]*ty < obj2max[0]*tx + obj2max[1]*ty + obj12T:
                        #the two collision points are obj2's vertices
                        obj1cp1 = (obj2min[0]+nx*minOverlap+obj12X, obj2min[1]+ny*minOverlap+obj12Y)
                        obj1cp2 = (obj2max[0]+nx*minOverlap+obj12X, obj2max[1]+ny*minOverlap+obj12Y)
                        obj2cp1 = obj2min
                        obj2cp2 = obj2max
                        obj1collidingVertex1 = None
                        obj2collidingVertex1 = obj2minID
                        obj1collidingVertex2 = None
                        obj2collidingVertex2 = obj2maxID
                    else:
                        #the two collision points are obj2min and obj1max
                        obj1cp1 = (obj2min[0]+nx*minOverlap+obj12X,obj2min[1]+ny*minOverlap+obj12Y)
                        obj1cp2 = obj1max
                        obj2cp1 = obj2min
                        obj2cp2 = (obj1max[0]-nx*minOverlap-obj12X,obj1max[1]-ny*minOverlap-obj12Y)
                        obj1collidingVertex1 = None
                        obj2collidingVertex1 = obj2minID
                        obj1collidingVertex2 = obj1maxID
                        obj2collidingVertex2 = None
                else:
                    #obj1max projection > obj2max projection
                    if obj1max[0]*tx + obj1max[1]*ty < obj2max[0]*tx + obj2max[1]*ty + obj12T:
                        #the two collision points are obj1min and obj2max
                        obj1cp1 = obj1min
                        obj1cp2 = (obj2max[0]+nx*minOverlap+obj12X,obj2max[1]+ny*minOverlap+obj12Y)
                        obj2cp1 = (obj1min[0]-nx*minOverlap-obj12X,obj1min[1]-ny*minOverlap-obj12Y)
                        obj2cp2 = obj2max
                        obj1collidingVertex1 = obj1minID
                        obj2collidingVertex1 = None
                        obj1collidingVertex2 = None
                        obj2collidingVertex2 = obj2maxID
                    else:
                        #the two collision points are obj1's vertices
                        obj1cp1 = obj1min
                        obj1cp2 = obj1max
                        obj2cp1 = (obj1min[0]-nx*minOverlap-obj12X,obj1min[1]-ny*minOverlap-obj12Y)
                        obj2cp2 = (obj1max[0]-nx*minOverlap-obj12X, obj1max[1]-ny*minOverlap-obj12Y)
                        obj1collidingVertex1 = obj1minID
                        obj2collidingVertex1 = None
                        obj1collidingVertex2 = obj1maxID
                        obj2collidingVertex2 = None
                            
                #we're done with this axis
                continue
                        

            '''else, there is only 1 MIN vertex of obj2'''
            #if there is only one obj1 MIN
            if obj1.projections[obj1AxisID].min2 == None:
                #then collision point is k
                if overlap < minOverlap:
                    minOverlap = overlap
                    minAxis = (nx, ny)
                    
                    obj2cp1 = obj2.vert[k]
                    obj1cp1 = (obj2cp1[0]+nx*minOverlap+obj12X, obj2cp1[1]+ny*minOverlap+obj12Y)
                    obj1cp2 = None
                    obj2cp2 = None
                    obj1collidingVertex1 = None
                    obj2collidingVertex1 = k
                    obj1collidingVertex2 = None
                    obj2collidingVertex2 = None

                    
                #done with this axis
                continue

                                    
            #else there are two obj1 MINs. Hence, check obj2 MAX 
            #traverse to find MAX
            while (j != j_final):
                j -= step #there is no need for j = (j+step+obj2.n)%obj2.n
                maxVal2 = obj2.vert[j][0]*nx + obj2.vert[j][1]*ny
                if maxVal2 < maxVal + GRAD_TOL:
                    j = (j+step)%obj2.n #there is no need for (j+step+obj2.n)%obj2.n
                    break
                else:
                    maxVal = maxVal2
                    maxVal2 -= 1000*GRAD_TOL #sentinel

            overlap2 = deltaMax + maxVal
            
            #if non overlapping, stop now
            if SKIN + overlap2 < 0:
#                print 8, overlap
                return 0 #-----------------------------------RETURN----------------

            #if it's about obj1 MIN and obj2 MAX
            if overlap2 < overlap and overlap2 < minOverlap:
                minOverlap = overlap2
                #NOTE REDEFINITION OF NORMAL
                nx = -nx
                ny = -ny
                minAxis = (nx, ny)
                
                #if there are two MAXes
                if maxVal2 > maxVal - GRAD_TOL:
                    #one MAX is j, the other is k
                    k = (j-step+obj2.n)%obj2.n
                    
                    #swap if necessary to have k be the next vertex if one travels CLOCKWISE from j
                    if step == 1:
                        k, j = j, k
                    
                    #Find collision points
                    #j and k must be MAXes, since there is only one MIN when projected wrt obj1's normal
                    #The points in local coordinates
                    #Projection of obj1's min2 on t is > obj1's min1, projection of obj2's max1 is > obj2's max2
                    #t for tangent == the minAxis rotated 90degree wrt z (CLOCKWISE).
                    obj1min = obj1.vert[obj1.projections[obj1AxisID].min1]
                    obj1max = obj1.vert[obj1.projections[obj1AxisID].min2]
                    obj2min = obj2.vert[k]
                    obj2max = obj2.vert[j]
                    
                    #FOR 2 COLLISION POINTS (copied-pasted from above. exactly the same)
                    #find and register the points
                    #t for tangent == minAxis rotated 90degree wrt z (CLOCKWISE).
                    tx = -ny
                    ty = nx
                    obj12T = obj12X*tx + obj12Y*ty

                    #obj1min projection < obj2min projection
                    if obj1min[0]*tx + obj1min[1]*ty < obj2min[0]*tx + obj2min[1]*ty + obj12T:
                        #obj1max projection > obj2max projection
                        if obj1max[0]*tx + obj1max[1]*ty < obj2max[0]*tx + obj2max[1]*ty + obj12T:
                            #the two collision points are obj2's vertices
                            obj1cp1 = (obj2min[0]+nx*minOverlap+obj12X, obj2min[1]+ny*minOverlap+obj12Y)
                            obj1cp2 = (obj2max[0]+nx*minOverlap+obj12X, obj2max[1]+ny*minOverlap+obj12Y)
                            obj2cp1 = obj2min
                            obj2cp2 = obj2max
                            obj1collidingVertex1 = None
                            obj2collidingVertex1 = k
                            obj1collidingVertex2 = None
                            obj2collidingVertex2 = j
                        else:
                            #the two collision points are obj2min and obj1max
                            obj1cp1 = (obj2min[0]+nx*minOverlap+obj12X,obj2min[1]+ny*minOverlap+obj12Y)
                            obj1cp2 = obj1max
                            obj2cp1 = obj2min
                            obj2cp2 = (obj1max[0]-nx*minOverlap-obj12X,obj1max[1]-ny*minOverlap-obj12Y)
                            obj1collidingVertex1 = None
                            obj2collidingVertex1 = k
                            obj1collidingVertex2 = obj1.projections[obj1AxisID].min2
                            obj2collidingVertex2 = None
                    else:
                        #obj1max projection > obj2max projection
                        if obj1max[0]*tx + obj1max[1]*ty < obj2max[0]*tx + obj2max[1]*ty + obj12T:
                            #the two collision points are obj1min and obj2max
                            obj1cp1 = obj1min
                            obj1cp2 = (obj2max[0]+nx*minOverlap+obj12X,obj2max[1]+ny*minOverlap+obj12Y)
                            obj2cp1 = (obj1min[0]-nx*minOverlap-obj12X,obj1min[1]-ny*minOverlap-obj12Y)
                            obj2cp2 = obj2max
                            obj1collidingVertex1 = obj1.projections[obj1AxisID].min1
                            obj2collidingVertex1 = None
                            obj1collidingVertex2 = None
                            obj2collidingVertex2 = j
                        else:
                            #the two collision points are obj1's vertices
                            obj1cp1 = obj1min
                            obj1cp2 = obj1max
                            obj2cp1 = (obj1min[0]-nx*minOverlap-obj12X,obj1min[1]-ny*minOverlap-obj12Y)
                            obj2cp2 = (obj1max[0]-nx*minOverlap-obj12X, obj1max[1]-ny*minOverlap-obj12Y)
                            obj1collidingVertex1 = obj1.projections[obj1AxisID].min1
                            obj2collidingVertex1 = None
                            obj1collidingVertex2 = obj1.projections[obj1AxisID].min2
                            obj2collidingVertex2 = None
                    
                    #mark the axis in obj2 as duplicate
                    duplicateAxisMask |= 1 << obj2.vertToNormal[j]
                        
                    #we're done with this axis
                    continue
                    
                #else, point of collision is the one MAX, j
                obj2cp1 = obj2.vert[j]
                obj1cp1 = (obj2cp1[0]+nx*minOverlap+obj12X, obj2cp1[1]+ny*minOverlap+obj12Y)
                obj1cp2 = None
                obj2cp2 = None
                obj1collidingVertex1 = None
                obj2collidingVertex1 = j
                obj1collidingVertex2 = None
                obj2collidingVertex2 = None
                
                #Done with this axis
            
            #if it's about obj1 MAX and obj2 MIN
            elif overlap < minOverlap:
                minOverlap = overlap
                minAxis = (nx, ny)
                
                obj2cp1 = obj2.vert[k]
                obj1cp1 = (obj2cp1[0]+nx*minOverlap+obj12X, obj2cp1[1]+ny*minOverlap+obj12Y)
                obj1cp2 = None
                obj2cp2 = None
                obj1collidingVertex1 = None
                obj2collidingVertex1 = k
                obj1collidingVertex2 = None
                obj2collidingVertex2 = None

            #else, we're done with this axis
        #end of [for obj1AxisID in xrange(obj1.nn)] loop

        '''for each axis in obj2'''
        #obj2 may have 2 MINs, which requires checking. but obj1 is guaranteed to only have 1 MAX and 1 MIN
        #NOTE: Here, when I said obj1 MIN, it's minimum projection along minAxis candidate. when I said obj2 MIN, it's MIN as stored in obj2.projections (which means
        #it is based on projection along obj2 normal, which is ANTI-parallel to minAxis candidate. The same applies for MAX.
        for obj2AxisID in xrange(obj2.nn):
            #if axis is marked, skip
            if duplicateAxisMask & 1:
                continue
            
            #this will be minAxis if collision is at obj2 MAX. if at obj2 MIN, minAxis is (-nx,-ny)
            nx = -obj2.normal[obj2AxisID][0]
            ny = -obj2.normal[obj2AxisID][1]
                                        
            #projection of the vector connecting the two bodies
            obj12N = obj12X*nx + obj12Y*ny
            
            #projection of first vertex of obj1
            k = 0
            minVal = obj1.vert[k][0]*nx + obj1.vert[k][1]*ny
            maxVal = minVal
            
            #2nd vertex
            k = 1
            value = obj1.vert[k][0]*nx + obj1.vert[k][1]*ny
            
            deltaMax = -obj12N+obj2.projections[obj2AxisID].max #when added with the MAX projection of obj1, gives overlap
            deltaMin = obj12N-obj2.projections[obj2AxisID].min #when subtracted by the MIN projection of obj1, gives overlap
            
            '''if to find MAX needs to traverse forward'''
            if value > maxVal:
                maxVal = value
                j = 1
                step = 1
                j_final = obj1.n-1
                #to find MIN if necessary
                k = obj1.n
                k_final = 2
            #else, traverse backward to find MAX
            else:
                minVal = value
                j = obj1.n
                step = -1
                j_final = 2
                #to find MIN if necessary
                #k == 1
                k_final = obj1.n-1

            #traverse to find MAX
            while (j != j_final):
                j += step #there is no need for j = (j+step+obj2.n)%obj2.n
                maxVal2 = obj1.vert[j][0]*nx + obj1.vert[j][1]*ny
                if maxVal2 < maxVal: #usually should be maxVal2 < maxVal + GRAD_TOL, but it is guaranteed obj1 has only 1 MAX
                    j = (j-step)%obj1.n #there is no need for (j-step+obj2.n)%obj2.n
                    break
                else:
                    maxVal = maxVal2

            overlap = deltaMax + maxVal
            
            #if non overlapping, stop now
            if SKIN + overlap < 0:
#                print 9, overlap
                return 0 #-----------------------------------RETURN----------------

            '''if there is only one obj2 MIN'''
            if obj2.projections[obj2AxisID].min2 == None:
                #then possible collision point is j
                if overlap < minOverlap:
                    minOverlap = overlap
                    minAxis = (nx, ny)
                    
                    obj1cp1 = obj1.vert[j]
                    obj2cp1 = (obj1cp1[0]-nx*minOverlap-obj12X, obj1cp1[1]-ny*minOverlap-obj12Y)
                    obj1cp2 = None
                    obj2cp2 = None
                    obj1collidingVertex1 = j
                    obj2collidingVertex1 = None
                    obj1collidingVertex2 = None
                    obj2collidingVertex2 = None
                    
                #done with this axis
                continue
                
            '''else, there are 2 obj2 MINs. Hence check obj1 MIN'''
            #traverse to find MIN
            while (k != k_final):
                k -= step #there is no need for k = (k-step+obj1.n)%obj1.n
                minVal2 = obj1.vert[k][0]*nx + obj1.vert[k][1]*ny
                if minVal2 > minVal: #usually should be minVal2 > minVal - GRAD_TOL, but it is guaranteed obj1 has only 1 MIN
                    k = (k+step)%obj1.n #there is no need for (k+step+obj1.n)%obj1.n
                    break
                else:
                    minVal = minVal2

            overlap2 = deltaMin - minVal
            
            #if non overlapping, stop now
            if SKIN + overlap2 < 0:
#                print 10, overlap
                return 0 #-----------------------------------RETURN----------------

            #if it's about obj2 MIN and obj1 MIN
            if overlap2 < overlap and overlap2 < minOverlap:
                minOverlap = overlap2
                #NOTE REDEFINITION OF NORMAL
                nx = -nx
                ny = -ny
                minAxis = (nx, ny)
                
                #point of collision is the one MIN, k
                obj1cp1 = obj1.vert[k]
                obj2cp1 = (obj1cp1[0]-nx*minOverlap-obj12X, obj1cp1[1]-ny*minOverlap-obj12Y)
                obj1cp2 = None
                obj2cp2 = None
                obj1collidingVertex1 = k
                obj2collidingVertex1 = None
                obj1collidingVertex2 = None
                obj2collidingVertex2 = None
                
                #Done with this axis
                
            #if it's about obj2 MAX and obj1 MIN
            elif overlap < minOverlap:
                minOverlap = overlap
                minAxis = (nx, ny)
                    
                obj1cp1 = obj1.vert[j]
                obj2cp1 = (obj1cp1[0]-nx*minOverlap-obj12X, obj1cp1[1]-ny*minOverlap-obj12Y)
                obj1cp2 = None
                obj2cp2 = None
                obj1collidingVertex1 = j
                obj2collidingVertex1 = None
                obj1collidingVertex2 = None
                obj2collidingVertex2 = None

            #else, we're done with this axis
            
            duplicateAxisMask >>= 1
        #end of [for obj2AxisID in xrange(obj2.nn)] loop
        
#        #OUTPUT PARAMETERS
#        minOverlap
#        #pushing obj2 along minAxis resolve collision
#        minAxis
#        #Collision point(s) in local coordinate. Note there can be up to two collision points
#        obj1cp1
#        obj2cp1 
#        obj1cp2
#        obj2cp2       
#        obj1collidingVertex1
#        obj2collidingVertex1
#        obj1collidingVertex2
#        obj2collidingVertex2
        #to be used by constrain solver
        
        
        '''register normal constraint'''
        nx = minAxis[0]
        ny = minAxis[1]
#        print obj1.idNum, obj2.idNum, minOverlap

        #Jacobian
        J = (-nx, -ny, obj1cp1[1]*nx-obj1cp1[0]*ny,\
              nx, ny, obj2cp1[0]*ny-obj2cp1[1]*nx)
        #Relative velocity of contact point in the normal direction, negative means approaching
        vrel = J[0]*obj1.speed[0] + J[1]*obj1.speed[1] + J[2]*obj1.speed[2] + J[3]*obj2.speed[0] + J[4]*obj2.speed[1] + J[5]*obj2.speed[2]
        
        self.Jmap.append((obj1.idNum, obj2.idNum))
        self.J.append(J)
        self.zeta.append(minOverlap*ERP/self.dt-RESTITUTION*vrel)
        self.lambdaMinMax.append((0, None))

#        self.Jmap2.append((obj1.idNum, obj2.idNum))
#        self.J2.append(J)
#        self.zeta2.append(minOverlap*ERP/self.dt)
#        self.kappaMinMax.append((0, None))
        
        #construct contactID. obj1 by convention has smaller idNum
#        if obj1.idNum > obj2.idNum:
#            #format of collision constraint contactID: contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex, obj2collidingVertex)
#            contactID = (obj2.idNum, obj1.idNum, obj2collidingVertex1, obj1collidingVertex1)
#        else:
        contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex1, obj2collidingVertex1)
                
        #save contactID of constraint present
        self.contactIDList.append(contactID)
#        self.contactIDList2.append(contactID)
        #increase counter
        self.numOfConstraint += 1
#        self.numOfConstraint2 += 1
             
        #Friction (no error correction)
        J = (ny, -nx, -obj1cp1[1]*ny-obj1cp1[0]*nx,\
             -ny, nx, obj2cp1[1]*ny + obj2cp1[0]*nx)
        self.Jmap.append((obj1.idNum, obj2.idNum))
        self.J.append(J)
        self.zeta.append(0)
        #if the normal constraint contactID is saved
        if contactID in self.storedLambdaValues:
            normal = self.storedLambdaValues[contactID]
            self.lambdaMinMax.append((-normal*FRICTION_COEFF, normal*FRICTION_COEFF))
        else:
            if obj1.massInv == 0 or (obj1.pos[1] > obj2.pos[1] and obj2.massInv != 0):
                normal = self.gravity/obj2.massInv
            else:
                normal = self.gravity/obj1.massInv
            self.lambdaMinMax.append((-normal*FRICTION_COEFF, normal*FRICTION_COEFF))
        #rename contactID. obj1 by convention has smaller idNum
#        if obj1.idNum > obj2.idNum:
#            #format of friction constraint contactID: contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex, obj2collidingVertex, 'f')
#            contactID = (obj2.idNum, obj1.idNum, obj2collidingVertex1, obj1collidingVertex1, 'f')
#        else:
        contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex1, obj2collidingVertex1, 'f')
            
        #save contactID of constraint present
        self.contactIDList.append(contactID)
        #increase counter
        self.numOfConstraint += 1
             
        #THIS CODE IS WRONG! COLLISION DETECTION NEEDS ONLY TO DETERMINE 1 POINT THEN DO THE CORRECT VERSION OF THIS USING CLIPPING (SEE CATTO'S GCD PRES)
        #find whether there is a second contact point if none exist yet
        if obj1cp2 == None:
            #if the collision is at obj1 vertices
            if obj2collidingVertex1 ==  None:
                #possible points
                a = obj1.nextVertex(obj1collidingVertex1)
                b = obj1.prevVertex(obj1collidingVertex1)
                #we would like to take largest one with largest projection along minAxis
                proja = obj1.vert[a][0]*nx + obj1.vert[a][1]*ny
                projb = obj1.vert[b][0]*nx + obj1.vert[b][1]*ny
                
                #make a the point to be considered
                if projb > proja:
                    proja = projb
                    a = b
                
                overlap = proja - nx*obj1cp1[0] - ny*obj1cp1[1] + minOverlap

                #if point a is inside obj2.
                if overlap + SKIN > 0:
                    obj1cp2 = (obj1.vert[a][0],obj1.vert[a][1])
                    obj2cp2 = (obj1cp2[0] - minAxis[0]*overlap - obj12X, obj1cp2[1] - minAxis[1]*overlap - obj12Y)
                    obj1collidingVertex2 = a

                    #Jacobian
                    J = (-nx, -ny, obj1cp2[1]*nx-obj1cp2[0]*ny,\
                          nx, ny, obj2cp2[0]*ny-obj2cp2[1]*nx)
                    #Relative velocity of contact point in the normal direction, negative means approaching
                    vrel = J[0]*obj1.speed[0] + J[1]*obj1.speed[1] + J[2]*obj1.speed[2] + J[3]*obj2.speed[0] + J[4]*obj2.speed[1] + J[5]*obj2.speed[2]
                    
                    self.Jmap.append((obj1.idNum, obj2.idNum))
                    self.J.append(J)
                    self.zeta.append(overlap*ERP/self.dt-RESTITUTION*vrel)
                    self.lambdaMinMax.append((0, None))

#                    self.Jmap2.append((obj1.idNum, obj2.idNum))
#                    self.J2.append(J)
#                    self.zeta2.append(overlap*ERP/self.dt)
#                    self.kappaMinMax.append((0, None))

                    #construct contactID. obj1 by convention has smaller idNum
#                    if obj1.idNum > obj2.idNum:
#                        #format: contactID = (obj1.idNum, obj2.idNum, )
#                        contactID = (obj2.idNum, obj1.idNum, obj2collidingVertex2, obj1collidingVertex2)
#                    else:
                    contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex2, obj2collidingVertex2)
                    
                    #save contactID of constraint present
                    self.contactIDList.append(contactID)
#                    self.contactIDList2.append(contactID)
                    #increase counter
                    self.numOfConstraint += 1
#                    self.numOfConstraint2 += 1
                    
                    #Friction
                    J = (ny, -nx, -obj1cp2[1]*ny-obj1cp2[0]*nx,\
                         -ny, nx, obj2cp2[1]*ny + obj2cp2[0]*nx)
                    self.Jmap.append((obj1.idNum, obj2.idNum))
                    self.J.append(J)
                    self.zeta.append(0)
                    #if the normal constraint contactID is saved
                    if contactID in self.storedLambdaValues:
                        normal = self.storedLambdaValues[contactID]
                        self.lambdaMinMax.append((-normal*FRICTION_COEFF, normal*FRICTION_COEFF))
                    else:
                        if obj1.massInv == 0 or (obj1.pos[1] > obj2.pos[1] and obj2.massInv != 0):
                            normal = self.gravity/obj2.massInv
                        else:
                            normal = self.gravity/obj1.massInv
                        self.lambdaMinMax.append((-normal*FRICTION_COEFF, normal*FRICTION_COEFF))
                    #rename contactID. obj1 by convention has smaller idNum
#                    if obj1.idNum > obj2.idNum:
#                        #format of friction constraint contactID: contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex, obj2collidingVertex, 'f')
#                        contactID = (obj2.idNum, obj1.idNum, obj2collidingVertex2, obj1collidingVertex2, 'f')
#                    else:
                    contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex2, obj2collidingVertex2, 'f')
                        
                    #save contactID of constraint present
                    self.contactIDList.append(contactID)
                    #increase counter
                    self.numOfConstraint += 1
                
            #else obj1collidingVertex1 == None
            else:
                #possible points
                a = obj2.nextVertex(obj2collidingVertex1)
                b = obj2.prevVertex(obj2collidingVertex1)
                #we would like to take largest one with largest projection along minAxis
                proja = obj2.vert[a][0]*nx + obj2.vert[a][1]*ny
                projb = obj2.vert[b][0]*nx + obj2.vert[b][1]*ny
                
                #make a the point to be considered
                if projb < proja:
                    proja = projb
                    a = b
                
                overlap = -proja + nx*obj2cp1[0] + ny*obj2cp1[1] + minOverlap

                #if point a is inside obj1.
                if overlap + SKIN > 0:
                    obj2cp2 = (obj2.vert[a][0],obj2.vert[a][1])
                    obj1cp2 = (obj2cp2[0] + minAxis[0]*overlap + obj12X, obj2cp2[1] + minAxis[1]*overlap + obj12Y)
                    obj2collidingVertex2 = a

                    #Jacobian
                    J = (-nx, -ny, obj1cp2[1]*nx-obj1cp2[0]*ny,\
                          nx, ny, obj2cp2[0]*ny-obj2cp2[1]*nx)
                    #Relative velocity of contact point in the normal direction, negative means approaching
                    vrel = J[0]*obj1.speed[0] + J[1]*obj1.speed[1] + J[2]*obj1.speed[2] + J[3]*obj2.speed[0] + J[4]*obj2.speed[1] + J[5]*obj2.speed[2]
                    
                    self.Jmap.append((obj1.idNum, obj2.idNum))
                    self.J.append(J)
                    self.zeta.append(overlap*ERP/self.dt-RESTITUTION*vrel)
                    self.lambdaMinMax.append((0, None))

#                    self.Jmap2.append((obj1.idNum, obj2.idNum))
#                    self.J2.append(J)
#                    self.zeta2.append(overlap*ERP/self.dt)
#                    self.kappaMinMax.append((0, None))

                    #construct contactID. obj1 by convention has smaller idNum
#                    if obj1.idNum > obj2.idNum:
#                        #format: contactID = (obj1.idNum, obj2.idNum, )
#                        contactID = (obj2.idNum, obj1.idNum, obj2collidingVertex2, obj1collidingVertex2)
#                    else:
                    contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex2, obj2collidingVertex2)
                    
                    #save contactID of constraint present
                    self.contactIDList.append(contactID)
#                    self.contactIDList2.append(contactID)
                    #increase counter
                    self.numOfConstraint += 1
#                    self.numOfConstraint2 += 1
                    
                    #Friction
                    J = (ny, -nx, -obj1cp2[1]*ny-obj1cp2[0]*nx,\
                         -ny, nx, obj2cp2[1]*ny + obj2cp2[0]*nx)
                    self.Jmap.append((obj1.idNum, obj2.idNum))
                    self.J.append(J)
                    self.zeta.append(0)
                    #if the normal constraint contactID is saved
                    if contactID in self.storedLambdaValues:
                        normal = self.storedLambdaValues[contactID]
                        self.lambdaMinMax.append((-normal*FRICTION_COEFF, normal*FRICTION_COEFF))
                    else:
                        if obj1.massInv == 0 or (obj1.pos[1] > obj2.pos[1] and obj2.massInv != 0):
                            normal = self.gravity/obj2.massInv
                        else:
                            normal = self.gravity/obj1.massInv
                        self.lambdaMinMax.append((-normal*FRICTION_COEFF, normal*FRICTION_COEFF))
                    #rename contactID. obj1 by convention has smaller idNum
#                    if obj1.idNum > obj2.idNum:
#                        #format of friction constraint contactID: contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex, obj2collidingVertex, 'f')
#                        contactID = (obj2.idNum, obj1.idNum, obj2collidingVertex2, obj1collidingVertex2, 'f')
#                    else:
                    contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex2, obj2collidingVertex2, 'f')
                        
                    #save contactID of constraint present
                    self.contactIDList.append(contactID)
                    #increase counter
                    self.numOfConstraint += 1
                
                
        #else, register the second contact point
        else:
#            print '2 points'
            #Jacobian
            J = (-nx, -ny, obj1cp2[1]*nx-obj1cp2[0]*ny,\
                  nx, ny, obj2cp2[0]*ny-obj2cp2[1]*nx)
            #Relative velocity of contact point in the normal direction, negative means approaching
            vrel = J[0]*obj1.speed[0] + J[1]*obj1.speed[1] + J[2]*obj1.speed[2] + J[3]*obj2.speed[0] + J[4]*obj2.speed[1] + J[5]*obj2.speed[2]
            
            self.Jmap.append((obj1.idNum, obj2.idNum))
            self.J.append(J)
            self.zeta.append(minOverlap*ERP/self.dt-RESTITUTION*vrel)
            self.lambdaMinMax.append((0, None))

#            self.Jmap2.append((obj1.idNum, obj2.idNum))
#            self.J2.append(J)
#            self.zeta2.append(minOverlap*ERP/self.dt)
#            self.kappaMinMax.append((0, None))

            #construct contactID. obj1 by convention has smaller idNum
#            if obj1.idNum > obj2.idNum:
#                #format: contactID = (obj1.idNum, obj2.idNum, )
#                contactID = (obj2.idNum, obj1.idNum, obj2collidingVertex2, obj1collidingVertex2)
#            else:
            contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex2, obj2collidingVertex2)
            
            #save contactID of constraint present
            self.contactIDList.append(contactID)
#            self.contactIDList2.append(contactID)
            #increase counter
            self.numOfConstraint += 1
#            self.numOfConstraint2 += 1
            
            #Friction
            J = (ny, -nx, -obj1cp2[1]*ny-obj1cp2[0]*nx,\
                 -ny, nx, obj2cp2[1]*ny + obj2cp2[0]*nx)
            self.Jmap.append((obj1.idNum, obj2.idNum))
            self.J.append(J)
            self.zeta.append(0)
            #if the normal constraint contactID is saved
            if contactID in self.storedLambdaValues:
                normal = self.storedLambdaValues[contactID]
                self.lambdaMinMax.append((-normal*FRICTION_COEFF, normal*FRICTION_COEFF))
            else:
                if obj1.massInv == 0 or (obj1.pos[1] > obj2.pos[1] and obj2.massInv != 0):
                    normal = self.gravity/obj2.massInv
                else:
                    normal = self.gravity/obj1.massInv
                self.lambdaMinMax.append((-normal*FRICTION_COEFF, normal*FRICTION_COEFF))
            #rename contactID. obj1 by convention has smaller idNum
#            if obj1.idNum > obj2.idNum:
#                #format of friction constraint contactID: contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex, obj2collidingVertex, 'f')
#                contactID = (obj2.idNum, obj1.idNum, obj2collidingVertex2, obj1collidingVertex2, 'f')
#            else:
            contactID = (obj1.idNum, obj2.idNum, obj1collidingVertex2, obj2collidingVertex2, 'f')
                
            #save contactID of constraint present
            self.contactIDList.append(contactID)
            #increase counter
            self.numOfConstraint += 1


#        print "time elapsed in calculation", time()-timestart

    def solveConstraint(self):
        '''Solve for lambda'''
        #given the constraint equation J M_inv J^T lambda = (zeta - J V)/dt
        
        '''query the cache of lambdaValues'''
        lambda0 = []
        
        #if contact is in store
        for contactID in self.contactIDList:
            if contactID in self.storedLambdaValues:
                lambda0.append(self.storedLambdaValues[contactID])
            else:
                lambda0.append(0)
#        print "bfr", lambda0

        '''PGS Solver'''
        
        '''initialization. has been checked with MATLAB'''
        #compute eta = (zeta - J V)/dt, while simultaneously
        eta = []
        #computing d, diagonal elements of J M_inv J^T, and
        d = []
        #a = M_inv J^T lambda0, a dictionary, its keys being idNums and its values 3x1 vectors
        a = {}
        for i in xrange(self.numOfConstraint):
            a[self.Jmap[i][0]] = [0,0,0]
            a[self.Jmap[i][1]] = [0,0,0]
        for i in xrange(self.numOfConstraint):
            b1 = self.Jmap[i][0]
            b2 = self.Jmap[i][1]
            sum = 0
            sum2 = 0
            #only b1 can be == 0 (immovable, hence don't have massInv, speed, etc.). 
            if b1 != 0:
                #computing a
                temp = self.entityList[b1].massInv * lambda0[i]
                a[b1][0] += self.J[i][0] * temp
                a[b1][1] += self.J[i][1] * temp
                a[b1][2] += self.J[i][2] * self.entityList[b1].inertiaInv * lambda0[i]
                #computing eta
                temp = self.entityList[b1].speed
                sum += self.J[i][0]*temp[0] + self.J[i][1]*temp[1] +self.J[i][2]*temp[2]
                #computing d
                sum2 += (self.J[i][0]**2+self.J[i][1]**2)*self.entityList[b1].massInv + self.J[i][2]**2*self.entityList[b1].inertiaInv
            #computing a
            temp = self.entityList[b2].massInv * lambda0[i]
            a[b2][0] += self.J[i][3] * temp
            a[b2][1] += self.J[i][4] * temp
            a[b2][2] += self.J[i][5] * self.entityList[b2].inertiaInv * lambda0[i]
            #computing eta
            temp = self.entityList[b2].speed
            sum += self.J[i][3]*temp[0] + self.J[i][4]*temp[1] +self.J[i][5]*temp[2]
            eta.append((self.zeta[i] - sum)/self.dt)
            #computing d
            sum2 += (self.J[i][3]**2+self.J[i][4]**2)*self.entityList[b2].massInv + self.J[i][5]**2*self.entityList[b2].inertiaInv
            d.append(sum2)
            
        '''iteration'''
        for iter in xrange(PGS_ITERATION):
            for i in xrange(self.numOfConstraint):
                b1 = self.Jmap[i][0]
                b2 = self.Jmap[i][1]
                dlambdai = RELAXATION*(eta[i] - self.J[i][0]*a[b1][0] - self.J[i][1]*a[b1][1] - self.J[i][2]*a[b1][2]\
                                   - self.J[i][3]*a[b2][0] - self.J[i][4]*a[b2][1] - self.J[i][5]*a[b2][2])/d[i]
                lambdai0 = lambda0[i]
                lambdamin = self.lambdaMinMax[i][0]
                lambdamax = self.lambdaMinMax[i][1]
                if lambdamax == None:
                    if lambdamin == None:
                        lambda0[i] += dlambdai
                    else:
                        lambda0[i] = max(lambdamin, lambdai0 + dlambdai)
                elif lambdamin == None:
                    lambda0[i] = min(lambdai0 + dlambdai, lambdamax)
                else:
                    lambda0[i] = max(lambdamin, min(lambdai0 + dlambdai, lambdamax))
                dlambdai = lambda0[i] - lambdai0
                a[b1][0] += dlambdai * self.J[i][0] * self.entityList[b1].massInv
                a[b1][1] += dlambdai * self.J[i][1] * self.entityList[b1].massInv
                a[b1][2] += dlambdai * self.J[i][2] * self.entityList[b1].inertiaInv
                a[b2][0] += dlambdai * self.J[i][3] * self.entityList[b2].massInv
                a[b2][1] += dlambdai * self.J[i][4] * self.entityList[b2].massInv
                a[b2][2] += dlambdai * self.J[i][5] * self.entityList[b2].inertiaInv
            
            print lambda0
            
        #store lambda0 values
        self.storedLambdaValues = {}
        for i in xrange(self.numOfConstraint):
            self.storedLambdaValues[self.contactIDList[i]] = lambda0[i]

#        print "aft", lambda0

        #apply force = J^T lambda0
        for i in xrange(self.numOfConstraint):
            b1 = self.Jmap[i][0]
            b2 = self.Jmap[i][1]
            lambdaT = lambda0[i]*self.dt
            #only b1 can be == 0 (immovable, hence don't have massInv, speed, etc.). 
            if b1 != 0:
                temp = self.entityList[b1] ##########################################################
                temp.speed[0] += self.J[i][0] * temp.massInv * lambdaT
                temp.speed[1] += self.J[i][1] * temp.massInv * lambdaT
                temp.speed[2] += self.J[i][2] * temp.inertiaInv * lambdaT
            temp = self.entityList[b2]
            temp.speed[0] += self.J[i][3] * temp.massInv * lambdaT
            temp.speed[1] += self.J[i][4] * temp.massInv * lambdaT
            temp.speed[2] += self.J[i][5] * temp.inertiaInv * lambdaT

#        for i in self.entityList:
#            print self.entityList[i].speed
#            
#            
#        '''POSITION-BASED ERROR CORRECTION'''
#        '''Solve for kappa'''
#        #given the constraint equation J2 M_inv J2^T kappa = zeta/dt
#        
#        '''query the cache of kappaValues, re-using the variable lambda0'''
#        lambda0 = []
#
#        #if contact is in store
#        for contactID in self.contactIDList2:
#            if contactID in self.storedKappaValues:
#                lambda0.append(self.storedKappaValues[contactID])
#            else:
#                lambda0.append(0)
#
#        print 'bfr kappa', lambda0
#
#        '''initialization. Reconstruct a and d'''
#        #computing d, diagonal elements of J2 M_inv J2^T, and
#        d = []
#        #a = M_inv J^T lambda0, a dictionary, its keys being idNums and its values 3x1 vectors
#        a = {}
#        for i in xrange(self.numOfConstraint2):
#            a[self.Jmap2[i][0]] = [0,0,0]
#            a[self.Jmap2[i][1]] = [0,0,0]
#        for i in xrange(self.numOfConstraint2):
#            b1 = self.Jmap2[i][0]
#            b2 = self.Jmap2[i][1]
#            sum = 0
#            #only b1 can be == 0 (immovable, hence don't have massInv, speed, etc.). 
#            if b1 != 0:
#                #computing a
#                temp = self.entityList[b1].massInv * lambda0[i]
#                a[b1][0] += self.J2[i][0] * temp
#                a[b1][1] += self.J2[i][1] * temp
#                a[b1][2] += self.J2[i][2] * self.entityList[b1].inertiaInv * lambda0[i]
#                #computing d
#                sum += (self.J2[i][0]**2+self.J2[i][1]**2)*self.entityList[b1].massInv + self.J2[i][2]**2*self.entityList[b1].inertiaInv
#            #computing a
#            temp = self.entityList[b2].massInv * lambda0[i]
#            a[b2][0] += self.J2[i][3] * temp
#            a[b2][1] += self.J2[i][4] * temp
#            a[b2][2] += self.J2[i][5] * self.entityList[b2].inertiaInv * lambda0[i]
#            #computing d
#            sum += (self.J2[i][3]**2+self.J2[i][4]**2)*self.entityList[b2].massInv + self.J2[i][5]**2*self.entityList[b2].inertiaInv
#            d.append(sum)
#
#        '''ErrorReduction iteration'''
#        for iter in xrange(ER_ITERATION):
#            for i in xrange(self.numOfConstraint2):
#                b1 = self.Jmap2[i][0]
#                b2 = self.Jmap2[i][1]
#                dlambdai = RELAXATION*(self.zeta2[i]/self.dt - self.J2[i][0]*a[b1][0] - self.J2[i][1]*a[b1][1] - self.J2[i][2]*a[b1][2]\
#                                   - self.J2[i][3]*a[b2][0] - self.J2[i][4]*a[b2][1] - self.J2[i][5]*a[b2][2])/d[i]
#                lambdai0 = lambda0[i]
#                print dlambdai
#                lambdamin = self.kappaMinMax[i][0]
#                lambdamax = self.kappaMinMax[i][1]
#                if lambdamax == None:
#                    if lambdamin == None:
#                        lambda0[i] += dlambdai
#                    else:
#                        lambda0[i] = max(lambdamin, lambdai0 + dlambdai)
#                elif lambdamin == None:
#                    lambda0[i] = min(lambdai0 + dlambdai, lambdamax)
#                else:
#                    lambda0[i] = max(lambdamin, min(lambdai0 + dlambdai, lambdamax))
#                dlambdai = lambda0[i] - lambdai0
#                a[b1][0] += dlambdai * self.J2[i][0] * self.entityList[b1].massInv
#                a[b1][1] += dlambdai * self.J2[i][1] * self.entityList[b1].massInv
#                a[b1][2] += dlambdai * self.J2[i][2] * self.entityList[b1].inertiaInv
#                a[b2][0] += dlambdai * self.J2[i][3] * self.entityList[b2].massInv
#                a[b2][1] += dlambdai * self.J2[i][4] * self.entityList[b2].massInv
#                a[b2][2] += dlambdai * self.J2[i][5] * self.entityList[b2].inertiaInv
#
#        #store kappa0 values
#        self.storedKappaValues = {}
#        for i in xrange(self.numOfConstraint2):
#            self.storedKappaValues[self.contactIDList2[i]] = lambda0[i]
#
#        print 'aft kappa', lambda0
#
#        #apply correction
#        for i in xrange(self.numOfConstraint2):
#            b1 = self.Jmap2[i][0]
#            b2 = self.Jmap2[i][1]
#            lambdaT = lambda0[i]*self.dt*self.dt
#            #only b1 can be == 0 (immovable, hence don't have massInv, speed, etc.). 
#            if b1 != 0:
#                temp = self.entityList[b1] ##########################################################
#                temp.pos[0] += self.J2[i][0] * temp.massInv * lambdaT
#                temp.pos[1] += self.J2[i][1] * temp.massInv * lambdaT
#                temp.pos[2] += self.J2[i][2] * temp.inertiaInv * lambdaT
#            temp = self.entityList[b2]
#            temp.pos[0] += self.J2[i][3] * temp.massInv * lambdaT
#            temp.pos[1] += self.J2[i][4] * temp.massInv * lambdaT
#            temp.pos[2] += self.J2[i][5] * temp.inertiaInv * lambdaT





    def pos_response(self, contact):
        ent1 = self.entityList[contact[0][0]]
        ent2 = self.entityList[contact[0][1]]
        
        proj_vec = contact[1]
        normal = (proj_vec[1][0].normal[proj_vec[1][1]][0]*proj_vec[1][2], proj_vec[1][0].normal[proj_vec[1][1]][1]*proj_vec[1][2])
        
            #world.addEntity(Smoke_gen((int(proj_vec[2][0]-5), int(proj_vec[2][1]-5), 10, 10), 0.1))

#for sticky substance, set e = 0, and no need to update position to remove overlap when colliding! :D

#we are projecting the bodies not until they completely separate (note the proj_vec[0]-skin). this is done to continuously record contact between objects
#to reduce jittering and allow more stacking. we do not modify position if proj_vec <= 1 or else objects will stick to each other.

#to move the object, undo their movement (instead of directly moving the object. this method causes inevitable sliding on inclined plane). but undo movement
#is equivalent to sweep test... will be back to this matter.
        if proj_vec[0] > skin:
            ent1.pos[0] += (proj_vec[0]-skin)*normal[0]*ent1.massInv/(ent1.massInv + ent2.massInv)
            ent1.pos[1] += (proj_vec[0]-skin)*normal[1]*ent1.massInv/(ent1.massInv + ent2.massInv)
            ent2.pos[0] -= (proj_vec[0]-skin)*normal[0]*ent2.massInv/(ent1.massInv + ent2.massInv)
            ent2.pos[1] -= (proj_vec[0]-skin)*normal[1]*ent2.massInv/(ent1.massInv + ent2.massInv)

    def vel_response(self, contact):

        ent1 = self.entityList[contact[0][0]]
        ent2 = self.entityList[contact[0][1]]
        
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
                m_inv1 = (ent1.massInv + ent2.massInv + ent1.inertiaInv*r_1pf**2 + ent2.inertiaInv*r_2pf**2) # > 0

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

        #m_inv1 and m_inv2 equals to 0 only if ent1 and ent2 are both immovablem, i.e. their .massInv and .inertiaInv are both 0.
        m_inv1 = (ent1.massInv + ent2.massInv + ent1.inertiaInv*r_1pf**2 + ent2.inertiaInv*r_2pf**2) # > 0
        m_inv2 = (ent1.massInv + ent2.massInv + ent1.inertiaInv*r_1pn**2 + ent2.inertiaInv*r_2pn**2) # > 0
        m_inv3 = (r_1pn*r_1pf*ent1.inertiaInv + r_2pn*r_2pf*ent2.inertiaInv)

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
        #ent1.speed += ent1.massInv * j*n.
        ent1.speed[0] += normal[0]*ent1.massInv*j 
        ent1.speed[1] += normal[1]*ent1.massInv*j
        #update speed change due to dynamic friction 
        ent1.speed[0] -= normal[1]*ent1.massInv*k
        ent1.speed[1] += normal[0]*ent1.massInv*k

        #update speed change due to normal force
        ent2.speed[0] -= normal[0]*ent2.massInv*j
        ent2.speed[1] -= normal[1]*ent2.massInv*j
        #update speed change due to dynamic friction
        ent2.speed[0] += normal[1]*ent2.massInv*k
        ent2.speed[1] -= normal[0]*ent2.massInv*k

        
        #ent1.speed[2] += (r_1p cross n) * ent1.inertiaInv * j. 
        ent1.speed[2] -= r_1pf*ent1.inertiaInv*j
        ent2.speed[2] += r_2pf*ent2.inertiaInv*j
        #ent1.speed[2] += (r_1p cross f) * ent1.inertiaInv * k. 
        ent1.speed[2] += r_1pn*ent1.inertiaInv*k
        ent2.speed[2] -= r_2pn*ent2.inertiaInv*k
        
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

