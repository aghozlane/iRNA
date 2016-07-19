"""
 @brief: Compute general data on list
 @author: Amine Ghozlane
 @version: 1.0 
"""
import numpy as np

class Analysis:
    def __init__(self):
        """
        Instanciate Analysis object
        """
        pass
    
    def commun_values(self,cumneg_nanvalues,normscore_nanvalues):
        """
        Detect commune nan values
        """
        commun=cumneg_nanvalues
        for j in xrange(len(commun)):
            if normscore_nanvalues[j]==True: commun[j]=True
        return commun
    
    def getRangePosit(self,array,col):
        """
        Get one column
        @param array: Matrix array of data
        @param col: Selected column 
        """
        return np.array([i[col] for i in array if i!=None])