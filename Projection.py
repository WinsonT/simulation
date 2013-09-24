class Projection():
    def __init__(self, max, pt1, min, pt3, pt2 = None, pt4 = None):
        self.max = max   #when vertices are projected onto a normal, the max value
        self.max1 = pt1  #index to the max vertex
        self.max2 = pt2  #or if there are 2 max vertices (None by default)
        self.min = min   #the min value
        self.min1 = pt3  #index to the min vertex
        self.min2 = pt4  #or if there are 2 min vertices (None by default)