## Ford Fishman
## Custom errors for sims

class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class FunctionError(Error):
    """Errors raised by an incorrect function call"""
    pass

# class apples():
#     def __init__(self):
#         self.thing = set()
    
#     def setThing(self, string):
#         self.thing.add(string)

# from copy import copy, deepcopy

# a = apples()
# b = deepcopy(a) 
# b.setThing("ad")
# print(a.thing)




