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


def generateName(type:Type, num):
    """Make names for different objects of different classes"""
    # name = type.value + str(num)
    nums = range(0,10)
    idList = np.random.choice(nums, size=10, replace=True)
    name = type.value + "".join( map(str,idList) )
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
        p = kw["p"]
        item_list = list()

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
        if n > 0:
            
            numEvents = np.random.binomial(n,p) 
            # print(str(lam) + " -> " + str(numEvents))
            # print(numEvents)
            # i = 0 # for iteration
            
            for i in range(0,numEvents):
                
                item_list.append( func(args[0],**kw, num=i) )
                # print(value)
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

a = range(0,10)
N=10**6
P=10**7
ab=10**-7
m=10**-5
lam = N*P*ab*m
n = N*P*ab
# def run(**kw):
#     print(kw["p"])
# run(p=1)
# print(np.random.poisson(lam =lam, size = 10))
# print(np.random.binomial(n=n,p=m,size=10))