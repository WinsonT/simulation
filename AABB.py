class AABB():
    def __init__(self, object, min = (0,0), max = (0,0)):
        self.object = object
        self.min = min
        self.max = max


'''in Polygon

    in __init__(...):
        ...
        
        #find top left and bottom right coordinates of AABB using center of mass as origin.
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
        self.AABB = AABB(self)
        self.updateAABB()
        
        self.world.addAABB(self.AABB)


    #AABB size update every time .pos[2] is updated
    def updateAABB(self):
        theta = self.pos[2]
        c = cos(theta)
        s = sin(theta)
        AABBx = self.pos[0] + c*self.midAABB0[0] - s*self.midAABB0[1]
        AABBy = self.pos[1] + c*self.midAABB0[1] + s*self.midAABB0[0]
        
        if theta > pi: #rotation by pi doesnt change width/height, so theta and theta-pi are equivalent
            theta -= pi
        if theta > pi/2.0: #finally rotation by pi/2 + alpha is equivalent to pi/2-alpha. effectively, in the end 0 < theta < pi/2
            theta = pi - theta
    
        c = cos(theta)
        s = sin(theta)
        halfWidth = (c*self.width0 + s*self.height0)/2
        halfHeight = (c*self.height0 + s*self.width0)/2
        
        self.AABB.min = (AABBx - halfWidth,AABBy - halfHeight)
        self.AABB.max = (AABBx + halfWidth,AABBy + halfHeight)
'''