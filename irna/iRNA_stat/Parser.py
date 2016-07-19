"""
 @brief: General treatment of dictionnary using bisect
 @author: Amine Ghozlane
 @version: 1.0 
"""
import bisect
from operator import itemgetter

class Parser:
    
    def __init__(self):
        """
        Instanciate Parser object
        """
        pass
    
    def sortdict(self,dictionnary,criteria):
        """
        Sort the dictionnary
        @param dictionnary: Dictionnary list
        @param criteria: Sort criteria
        """
        return sorted(dictionnary,key=itemgetter(criteria))
    
    def getData(self,input_list,name):
        """
        Search name in input list
        @param input_list: List
        @param name: Search criteria
        """
        #Searching the node with its name
        i=bisect.bisect_left(input_list,name)
        #Object has been found
        if(i!=len(input_list) and input_list[i]==name): return i
        return None
    
    def getUnique(self,data_list):
        """
        Get unique data
        @param data_list: list of data
        """
        # Dave Kirby
        return {}.fromkeys(data_list).keys()