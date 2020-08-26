## Ford Fishman

"""
Functions and classes used across classes
"""
from Enums import Type, Mutation
import numpy as np
import pandas as pd

def generateName(type:Type):
    """Make names for different objects of different classes"""
    # name = type.value + str(num)
    nums = range(0,10)
    idList = np.random.choice(nums, size=10, replace=True)
    nameEnd = "".join( map(str,idList) )
    name = "%s%s" % (type.value,nameEnd)
    return name

def initRecord():
    """
    returns a dataframe that records information to this phage population
    """
    record = pd.DataFrame(None, columns=["timestep", "name", "pop", "dpop", "dpop_pop","type", "spacers"])
    
    return record

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

    

        
import pandas as pd
# import numpy as np

# d = {
#     "a":{"1":True,"2":True,"3":False},
#     "b":{"1":True,"2":True,"3":False},
#     "c":{"1":True,"2":True,"3":False}
# }

# s = pd.Series(data = {"a":1,"b":2,"c":3})
# print(s)
# s["a"] = "beeg"
# print(s["a"])

# def func(val):
#     b = {"5":"five","6":"six","4":"four"}
#     n = val.index
#     # print(val.where(val is dict,"not dict"))


import time
from tqdm import tqdm



# bar = pb.ProgressBar(maxval=100, widgets=widgets).start()

for i in tqdm(range(100)):
    time.sleep(0.1)
    # bar.update(i)


df = pd.DataFrame(None, columns=("a","b","c"))
df.loc[0] = df.columns.map(pd.Series({"a":1,"b":3,"c":3}))
df.loc[1] = df.columns.map(pd.Series({"b":29,"a":3,"c":3}))

# print(df)

# df2 = pd.DataFrame(None, columns=("a","b","c"))
# df2.loc[0] = df.columns.map(pd.Series({"a":1,"b":3,"c":3}))

# df2 = df2.append(df, ignore_index = True )
# print(df2)
# df1 = df.apply(lambda col:func(col),axis=0)
# print(df)
# import argparse

# parser = argparse.ArgumentParser()

# parser.add_argument('-p')
# parser.add_argument('-pS')
# parser.add_argument('--adsp')
# parser.add_argument('-o','--output', default='comim', type=str, help="Desired name of output path")

# arg = parser.parse_args(["-pS","2"])
# print(arg.pS)

# a = initRecord()

# a.loc[0] = [1, 'name', 1,2,"N"]
# print(a)