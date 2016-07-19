# -*- coding: utf-8 -*-
"""
 @brief: Normalize score
 @author: Amine Ghozlane
 @version: 1.0 
"""
import numpy as np
from Communication import *
  
class NormScore:
    
    def __init__(self,data,dbmanage):
        """
        Instanciate NormScore object
        @param data: Communication object
        @param dbmanage: Access to the database
        """
        self.normfuncmap=[self.normSeveralScore,self.normSeveralEnergy,self.normSeveralContactScore,self.normUniqueSolution]
        self.interactions=data.interactions
        self.type_sol=data.type_sol
        self.score_type=data.score_type
        self.dbmanage=dbmanage               
    
    def defnormfunction(self):
        """
        Determine the function to apply depending on the type of software
        """
        if(self.type_sol==0): return self.normfuncmap[0]
        elif(self.type_sol==1): return self.normfuncmap[1]
        elif(self.type_sol==2): return self.normfuncmap[2]
        else: return self.normfuncmap[3]       
        
    def normSeveralScore(self,score,loglen,score_type,i):
        """
        Compute normalised score in several score context
        @param score: Score
        @param loglen: Log length value
        @param score_type: Type of score
        @param i: Position in the tab 
        """
        #Plusieurs solutions - score
        temp_score=score/loglen
        self.pos[i]=np.argmax(temp_score)
        self.norm_score[i]=temp_score[self.pos[i]]
            
    def normSeveralEnergy(self,score,loglen,score_type,i):
        """
        Compute normalised energy in several energy context
        @param score: Score
        @param loglen: Log length value
        @param score_type: Type of score
        @param i: Position in the tab 
        """
        #Plusieurs solutions - energy
        #Normalisation
        temp_score=score/loglen
        self.pos[i]=np.argmin(temp_score)
        #Cas yass blastall
        if(score_type==1): self.norm_score[i]=score[self.pos[i]]
        else: self.norm_score[i]=temp_score[self.pos[i]]
                  
    def normSeveralContactScore(self,score,loglen,score_type,i):
        """
        Compute normalised energy in several energy context
        @param score: Score
        @param loglen: Log length value
        @param score_type: Type of score
        @param i: Position in the tab 
        """
        self.pos[i]=-1
        #Plusieurs contacts - score
        self.norm_score[i]=sum(score)/loglen
        
    def normUniqueSolution(self,score,loglen,score_type,i):
        """
        Compute normalised score/energy in several contact context
        @param score: Score
        @param loglen: Log length value
        @param score_type: Type of score
        @param i: Position in the tab 
        """
        self.pos[i]=-1
        #Plusieurs contacts - energy 3 ou une seule solution 4
        self.norm_score[i]=score[0]/loglen
    
    def run(self):
        """
        Compute normalized score
        """
        self.norm_score=np.arange(float(len(self.interactions)))
        self.pos=np.arange(len(self.interactions))
        #Score to normalise
        #listscore=self.dbmanage.getAllScore(self.interactions)
        #Len sRNA
        listlensRNA=self.dbmanage.getAllsRNAlenbyIntid(self.interactions)
        listlenmRNA=self.dbmanage.getAllmRNAlenbyIntid(self.interactions)
        funcspec=self.defnormfunction()
        for i in xrange(len(self.interactions)):
            if((i%500)==0): print("%d/%d"%(i+1,len(self.interactions)))
            score=self.dbmanage.getScore(self.interactions[i][0])
            #if(listscore[i]!=None):
            #    funcspec(listscore[i],np.log(float(listlensRNA[i]*listlenmRNA[i])),self.score_type,i)   
            if(score!=None):
                funcspec(score,np.log(float(listlensRNA[i]*listlenmRNA[i])),self.score_type,i)
                #Cas mfe / energy            
                if(self.score_type==2 or self.score_type==3): self.norm_score[i]=-self.norm_score[i]
            #Nothing is predicted
            else: self.norm_score[i]=np.nan 
        #Communication object
        data=Communication()
        data.norm_score=self.norm_score
        data.pos=self.pos
        return data
        