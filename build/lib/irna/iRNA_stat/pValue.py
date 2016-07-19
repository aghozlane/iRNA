"""
 @brief: Compute pValue of the interaction
 @author: Amine Ghozlane
 @version: 1.0 
"""
import numpy as np
from Communication import *

class pValue:
    
    def __init__(self,data,dbmanage,rand_inf):
        """
        Instanciate pValue object
        @param data: Data object
        @param dbmanage: Access to the database
        @param rand_inf: Rand_data object
        """
        self.dbmanage=dbmanage  
        self.rand_inf=rand_inf
        self.softname=data.softname
        self.norm_score=data.norm_score
        self.interactions=data.interactions
                   
    def run(self):
        """
        Compute pValue
        """
        i=0
        pValue=np.arange(float(len(self.interactions)))
        for interactid in self.interactions:
            #Take random data
            ncomp=self.rand_inf.getRandRna(self.softname,self.dbmanage.getsRNAbyIntid(interactid[0]))
            if(ncomp!=None):
                if(self.rand_inf.slope[ncomp] and self.rand_inf.intercept[ncomp]):
                    scale=-1.0/self.rand_inf.slope[ncomp]
                    location=self.rand_inf.intercept[ncomp]*scale
                    pValue[i]=1.0-np.exp(-np.exp(-((self.norm_score[i]-location)/scale) ))             
                else: pValue[i]=np.nan
            else:
                print("No data for %s with %s"%(self.dbmanage.getsRNAbyIntid(interactid[0]),self.softname))
                pValue[i]=np.nan
            i+=1
        #Communication object
        data=Communication()
        data.pValue=pValue
        return data