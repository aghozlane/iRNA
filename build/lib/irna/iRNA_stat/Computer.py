"""
 @brief: Handle software prediction
 @author: Amine Ghozlane
 @version: 1.0 
"""
import numpy as np

class Computer:
    
    def __init__(self, softid,name,numsoft,type_sol,score_type,dbmanage,args):
        """
        Instanciate Computer object
        @param softid: Id of the software
        @param name: Name of the software
        @param numsoft: Corresponding numsoft in soft_inf
        @param type_sol: Type of solution
        @param score_type: Type of score
        @param dbmanage: Access to the database
        @param args: Arguments
        """
        self.ppv=None
        self.sensitivity=None
        self.classtype=None
        #self.pValue=None
        self.curve_param=None        
        self.softid=softid
        self.numsoft=numsoft
        self.name=name
        self.cumneg=None
        self.curve_param=[]
        self.type_sol=type_sol
        self.score_type=score_type
        self.pValue_type=args.pValue
        #Number of interaction predict with one software
        self.nbinteract=dbmanage.getNbInteract(self.softid)
        #Interaction id of one software
        self.interactions=dbmanage.getInteract(self.softid)
        #sRNAid in interaction
        self.sRNAid_tab=dbmanage.getAllsRNAidbyIntid(self.interactions)
        self.unique_sRNAidinint=np.unique(self.sRNAid_tab)
        #NormScore
        self.norm_score=np.arange(float(self.nbinteract))
        self.pos=np.arange(self.nbinteract)
        #Interaction
        if(args.exp_inf):
            self.ppv=np.arange(float(self.nbinteract))
            self.sensitivity=np.arange(float(self.nbinteract))
            self.classtype=np.arange(self.nbinteract)
        #pValue
        #if(args.rand_inf):
        self.pValue=np.arange(float(self.nbinteract))
        self.select=[False]*self.nbinteract
        self.interact=[[None]]*self.nbinteract
        #Regression
        if(args.random):
            self.curve_param=[[None]]*len(self.unique_sRNAidinint)
            self.cumneg=[None]*len(self.unique_sRNAidinint)