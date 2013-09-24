
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
        a = [j for j in i]
        copy.append(a)
    return copy

#blankcopy creates a copy of list2, with each component replaced by a new value "blank" (blank must not be an object)
#note that list2 can be tuple2
def blankcopy_list2(list2, blank):
    copy = []
    for i in list2:
        a = [blank for j in i]
        copy.append(a)
    return copy

#for a list in a list in a list
#note that list3 can be tuple3
def copy_list3(list3):
    copy = []
    for i in list3:
        a = []
        for j in i:
            b = [k for k in j]
            a.append(b)
        copy.append(a)
    return copy

#convert list in a list into tuple in tuple
def tuple2(list2):
    a = [tuple(i) for i in list2]
    return tuple(a)

def tuple3(list3):
    a = []
    for i in list3:
        b = [tuple(j) for j in i]
        a.append(tuple(b))
    return tuple(a)

    
