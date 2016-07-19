"""
 @brief: Analyse ecdf data
 @author: Amine Ghozlane
 @version: 1.0 
"""
import numpy as np
from scipy import stats 
from Ecdf import Ecdf
from Analysis import *
from Communication import *

class Regression(Analysis):
    
    def __init__(self,data):
        """
        Instanciate Regression object
        @param data: Communication object
        """
        Analysis.__init__(self)
        self.unique_sRNAidinint=data.unique_sRNAidinint
        self.norm_score=data.norm_score
        self.sRNAid_tab=data.sRNAid_tab
        self.pValue_type=data.pValue_type
    
    def ecdf_estimate(self):
        """
        Compute empirical cumulated distribution 
        """
        self.cumneg=None
        #sRNA studied
        for sRNAid in self.unique_sRNAidinint:
            #Score of one sRNA
            numpair=np.where(self.sRNAid_tab==sRNAid)
            # Compute empirical cumulative distribution function
            #Observation of all sRNA studied
            if(self.pValue_type==1): fn=Ecdf(self.norm_score)                
            #Observation of current sRNA
            elif(self.pValue_type==2): fn=Ecdf(self.norm_score[numpair])
            else: raise ValueError            
            y=np.array([fn(score) for score in self.norm_score[numpair]]) 
            # Normalization
            y[np.where(y==0.0)]=np.nan
            y=-np.log(y)
            y[np.where(y==0.0)]=np.nan
            if(self.cumneg!=None): self.cumneg+=[np.array(np.log(y))]
            else: self.cumneg=[np.array(np.log(y))]
        
    def linear_regression(self):
        """
        Compute linear regression for each sRNA
        """
        i=0
        self.curve_param=[]
        #Vecteur de parametre        
        for srnaid in self.unique_sRNAidinint:
            slope=np.nan
            intercept=np.nan
            #Score of one sRNA
            numpair=np.where(self.sRNAid_tab==srnaid)
            cumnegtemp=self.cumneg[i]
            normscoretemp=self.norm_score[numpair]
            #Nan values
            cumneg_nanvalues=np.isnan(cumnegtemp)
            normscore_nanvalues=np.isnan(normscoretemp)
            #Filter commun nan values
            commun=self.commun_values(cumneg_nanvalues,normscore_nanvalues)
            if(False in commun):
                #only write values
                cumnegtemp=cumnegtemp[np.where(commun==False)]
                normscoretemp=normscoretemp[np.where(commun==False)]
                #Linear regression
                slope,intercept, r_value, p_value, std_err=stats.linregress(normscoretemp,cumnegtemp)
            self.curve_param+=[np.array([slope,intercept])]
            i+=1

    def run(self):
        #Ecdf estimation
        self.ecdf_estimate()        
        #Linear regression
        self.linear_regression()
        #Result
        data=Communication()
        data.curve_param=self.curve_param
        data.cumneg=self.cumneg
        return(data)