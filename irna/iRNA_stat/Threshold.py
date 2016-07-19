"""
 @brief: Handle significativity threshold information
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv, bisect
from Parser import *

class Threshold(Parser):
    
    def __init__(self,thres_file):
        """
        Instanciate soft object
        """
        Parser.__init__(self)
        self.getThresinf(thres_file)
    
    def getThresinf(self,thres_file):
        """
        Parse pValue threshold file
        """
        try:
            thes_infReader = csv.reader(open(thres_file, 'rb'), delimiter='\t')
            #pass header
            thes_infReader.next()
            thresdict=[{"soft":line[0].strip().lower(),"threshold":float(line[3])} for line in thes_infReader]
            #Sort of dictionnary list
            thresdict=self.sortdict(thresdict,"soft")
            #Get data in separate item
            self.soft= [execution['soft'] for execution in thresdict]
            self.threshold= [execution['threshold'] for execution in thresdict]
        except IOError:
            sys.exit("Error : can not open file %s"%thres_file)
            
    def getSoftnum(self,soft):
        """
        @param soft: Name of asked software 
        """
        #Searching the node with its name
        i=bisect.bisect_left(self.soft,soft.lower())
        #Object has been found
        if(i!=len(self.soft) and self.soft[i]==soft.lower()): return i
        return None