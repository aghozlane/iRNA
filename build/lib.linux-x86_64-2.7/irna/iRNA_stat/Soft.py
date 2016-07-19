"""
 @brief: Handle software information
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv, bisect
from Parser import *

class Soft(Parser):
    
    
    def __init__(self,soft_file):
        """
        Instanciate soft object
        """
        Parser.__init__(self)
        self.getSoftinf(soft_file)
        
    def getSoftinf(self,soft_file):
        """
        @param soft_file: 
        """
        try:
            soft_infReader = csv.reader(open(soft_file, 'rb'), delimiter='\t')
            #pass header
            soft_infReader.next()
            softdict=[{"id":line[0].strip().lower(),"type_sol":int(line[1]),"score_type":int(line[2])} for line in soft_infReader]
            #Sort of dictionnary list 
            softdict=self.sortdict(softdict,"id")
            #Get data in separate item
            self.id= [software['id'] for software in softdict]
            self.type_sol= [software['type_sol'] for software in softdict]
            self.score_type= [software['score_type'] for software in softdict]
        except IOError: sys.exit("Error : can not open file %s"%soft_file)
        #except:
        #    sys.exit("Something went wrong with %s"%self.soft_inf)
    
    def getSoftnum(self,name):
        """
        @param name: 
        """
        #Searching the node with its name
        i=bisect.bisect_left(self.id,name.lower())
        #Object has been found
        if(i!=len(self.id) and self.id[i]==name.lower()): return i
        return None
    
    def __cmp__(self,other):
        """
        General method to compare node based on the name
        @param other: Compared value
        """
        if self.id<other: return -1
        elif self.id==other: return 0
        elif self.id>other: return 1
        else: raise ValueError