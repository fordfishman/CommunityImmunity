## Ford Fishman

"""
Functions and classes used across classes
"""
from Enums import Type, Mutation
import numpy as np


# def switch(type:Type):
#     """
#     """
#     return type.name


def generateName(type:Type):
    """Make names for different objects of different classes"""
    # name = type.value + str(num)
    nums = range(0,10)
    idList = np.random.choice(nums, size=10, replace=True)
    nameEnd = "".join( map(str,idList) )
    name = "%s%s" % (type.value,nameEnd)
    return name


def dispatch_on_value(func):
    """
    Value-dispatch function decorator.
    
    Transforms a function into a value-dispatch function,
    which can have different behaviors based on the value of the first argument.
    http://hackwrite.com/posts/learn-about-python-decorators-by-writing-a-function-dispatcher/
    """
    
    registry = {}

    def dispatch(value):

        try:
            return registry[value]
        except KeyError:
            return func

    def register(value, func=None):
       
        if func is None:
            return lambda f: register(value, f)
        
        registry[value] = func
        
        return func

    def wrapper(*args, **kw):
        return dispatch(args[1])(*args, **kw)

    wrapper.register = register
    wrapper.dispatch = dispatch
    wrapper.registry = registry

    return wrapper

def runProcess(func):
    """
    runs a function based on an expected number of times
    The actual number of times the process is run depends on 
    a binomial process, where n is the number of interactions, 
    and p is the probability the process is run per interaction
    """
    def wrapper(*args,**kw):

        n = 1.0 # default expected number of events
        p = kw.get("p",0.0)

        for num in args[1:len(args)]:

            # try: 
            #     float(num)

            # except TypeError:

            #     print("Arg {} not tranformable to float. Treating as 1.0.".format(num))
            #     print()
            #     num = 1.0 # if arg is not a number, it won't modify n

            # else: 
            #     num = float(num)

            # finally: 
            #     print(num)
            n *= num 
        item_list = []
        if n > 0 and p>0:
            
            numEvents = np.random.binomial(n,p) 
            
            item_list = [func(args[0],**kw) for i in range(numEvents)]
            
        return item_list

    return wrapper

@runProcess
def testFunc(*args, b=""):
    return b

    

        
# print(testFunc(2,3, b = "Asd"))

# NUCLEOTIDES = ("A","C","G","T")

# class Test():

#     @dispatch_on_value
#     def mutate(self, mutation, genome:str):
#         pass # only runs if a mutation occurs that's not in thhe class
        

#     @mutate.register(Mutation.SNP)
#     def _(self, mutation, genome:str):

#         nt = np.random.choice(NUCLEOTIDES) # the new nucleotide

#         gLength = len(genome) # genome length
#         i = np.random.choice( range(0, gLength) )
#         # make genome into a list to change position
#         genomeList = list(genome) 
#         genomeList[i] = nt
#         newGenome = "".join(genomeList)

#         return newGenome

#     @mutate.register(Mutation.DELETION)
#     def _(self, mutation, genome:str):
        
#         gLength = len(genome) # genome length
#         i = np.random.choice( range(0, gLength) )
#         # make genome into a list to delete position
#         genomeList = list(genome)
#         genomeList.pop(i)
#         newGenome = "".join(genomeList)

#         return newGenome

# import timeit
# a = list(range(1000000))
# c = list()
# t0 = timeit.default_timer()
# for i in a:
#     x = i *4
#     c.append(x**2)
#     # c.append(i*3)
# t1 = timeit.default_timer()
# print(t1-t0)

# # a = [0]
# c = list()

# t0 = timeit.default_timer()
# # b = [i*2 for i in a]
# c = [(i*3)**2 for i in a]

# t1 = timeit.default_timer()
# print(t1-t0)
