"""
 @brief: Compute ppv and sensitivity of interaction prediction
 @author: Amine Ghozlane
 @version: 1.0 
"""
import numpy as np
from Communication import *
from Analysis import *

class Interaction(Analysis):
    
    def __init__(self,data,dbmanage,exp_inf):
        """
        Instanciate Interaction object
        @param data: Communication object
        @param dbmanage: Access to the database
        @param exp_inf: Exp_inf object
        """
        Analysis.__init__(self)
        self.intefuncmap=[self.intebestsol,self.inteseveralcontact,self.inteonesol]
        self.interactions=data.interactions
        self.type_sol=data.type_sol
        #Matrice partielle
        self.pos=data.pos
        self.exp_inf=exp_inf
        self.dbmanage=dbmanage
        
    def defintefunction(self):
        """
        Determine the function to apply depending of the type of solution
        """
        if(self.type_sol==0 or self.type_sol==1): return self.intefuncmap[0]
        elif(self.type_sol==2 or self.type_sol==3): return self.intefuncmap[1]
        else: return self.intefuncmap[2] 
    
    def intebestsol(self,posit,sRNA_deb,sRNA_end,mRNA_deb,mRNA_end):
        """
        Case best solution
        @param posit: Position in the list
        @param sRNA_deb: sRNA begin list
        @param sRNA_end: sRNA end list
        @param mRNA_deb: sRNA begin list
        @param mRNA_end: mRNA end list
        @return: Position selected and normalised
        """
        max_len=np.max([sRNA_end[posit]-sRNA_deb[posit],mRNA_end[posit]-mRNA_deb[posit]])
        if(max_len==0): max_len=1
        sRNA=np.arange(sRNA_deb[posit],sRNA_deb[posit]+max_len+1)
        mRNA=np.arange(mRNA_deb[posit]+max_len,mRNA_deb[posit]-1,-1)
        if(sRNA!=None and mRNA!=None): return np.array([sRNA,mRNA])
        else: return None
    
    def inteseveralcontact(self,posit,sRNA_deb,sRNA_end,mRNA_deb,mRNA_end):
        """
        Case several contacts score - several contacts energy
        @param posit: Position in the list
        @param sRNA_deb: sRNA begin list
        @param sRNA_end: sRNA end list
        @param mRNA_deb: sRNA begin list
        @param mRNA_end: mRNA end list
        @return: Position selected and normalised
        """
        inte=None
        for i in xrange(len(sRNA_deb)):
            max_len=np.max([sRNA_end[i]-sRNA_deb[i],mRNA_end[i]-mRNA_deb[i]])
            if(max_len==0): max_len=1
            sRNA=np.arange(sRNA_deb[i],sRNA_deb[i]+max_len+1)
            mRNA=np.arange(mRNA_deb[i]+max_len,mRNA_deb[i]-1,-1)
            temp=np.array([sRNA,mRNA])
            #Add to the matrix
            if(inte!=None):
                inte=np.concatenate((inte, temp), axis=1)
            else:
                inte=temp
        return inte
    
    def inteonesol(self,posit,sRNA_deb,sRNA_end,mRNA_deb,mRNA_end):
        """
        Case : one solution
        @param posit: Position in the list
        @param sRNA_deb: sRNA begin list
        @param sRNA_end: sRNA end list
        @param mRNA_deb: sRNA begin list
        @param mRNA_end: mRNA end list
        @return: Position selected and normalised
        """
        max_len=np.max([sRNA_end-sRNA_deb,mRNA_end-mRNA_deb])
        if(max_len==0): max_len=1
        sRNA=np.arange(sRNA_deb,sRNA_deb+max_len+1)
        mRNA=np.arange(mRNA_deb+max_len,mRNA_deb-1,-1)
        if(sRNA!=None and mRNA!=None): return np.array([sRNA,mRNA])
        else: return None
                       
    def getInteraction(self,posit,sRNA_deb,sRNA_end,mRNA_deb,mRNA_end,predict,intefuncpoint):
        """
        Compute interaction position
        @param posit: Position in the list
        @param sRNA_deb: sRNA begin list
        @param sRNA_end: sRNA end list
        @param mRNA_deb: sRNA begin list
        @param mRNA_end: mRNA end list
        @param predict: Flag Predict position
        @param intefuncpoint: Function pointer to analysis function
        @return: Position selected and normalised
        """
        inte=None
        #Predicted data
        if(predict): inte=intefuncpoint(posit,sRNA_deb,sRNA_end,mRNA_deb,mRNA_end)
        #Case real data    
        else:
            max_len=np.max([sRNA_end-sRNA_deb,mRNA_end-mRNA_deb])
            if(max_len==0): max_len=1
            sRNA=np.arange(sRNA_deb,sRNA_deb+max_len+1)
            mRNA=np.arange(mRNA_deb+max_len,mRNA_deb-1,-1)
            inte=np.array([sRNA,mRNA])        
        return inte
    
    #Number of correctly predicted base pairings
    def ncbp(self,pRNA,rRNA):
        """
        Compute the number of correctly predicted base pairings
        """
        temp=np.concatenate((pRNA,rRNA),axis=0)
        leninit=len(temp)
        lenfinal=len(np.unique(temp))
        return(leninit-lenfinal)
    
    def compute_ppv(self,Pinteract,Rinteract):
        """
        Compute PPV
        @param Pinteract: Predicted interaction position
        @param Rinteract: Real interaction position
        @return: PPV of one interation
        """
        sRNA_ncbp=self.ncbp(Pinteract[0],Rinteract[0])
        mRNA_ncbp=self.ncbp(Pinteract[1],Rinteract[1])
        sRNA_npbp=float(len(Pinteract[0]))
        mRNA_npbp=float(len(Pinteract[1]))
        return(np.mean(np.array([sRNA_ncbp/sRNA_npbp,mRNA_ncbp/mRNA_npbp])))
    
    def compute_sensitivity(self,Pinteract,Rinteract):
        """
        Compute sensitivity
        @param Pinteract: Predicted interaction position
        @param Rinteract: Real interaction position
        @return: sensitivity of one interation
        """
        sRNA_ncbp=self.ncbp(Pinteract[0],Rinteract[0])
        mRNA_ncbp=self.ncbp(Pinteract[1],Rinteract[1])
        sRNA_ntbp=float(len(Rinteract[0]))
        mRNA_ntbp=float(len(Rinteract[1]))
        return(np.mean(np.array([sRNA_ncbp/sRNA_ntbp,mRNA_ncbp/mRNA_ntbp])))    
    
    def compute_sens_ppv(self):
        """
        Compute sensitivity and ppv of the interaction
        """
        self.ppv=np.arange(float(len(self.interactions)))
        self.sensitivity=np.arange(float(len(self.interactions)))
        self.classtype=np.arange(len(self.interactions))
        intefunc=self.defintefunction()
        i=0
        for interactid in self.interactions:
            #Position d'interaction sRNA mRNA
            predict_RNA=self.dbmanage.getPositionsbyIntid(interactid[0])
            ncomp=self.exp_inf.getExpPair(self.dbmanage.getsRNAbyIntid(interactid[0]),self.dbmanage.getmRNAbyIntid(interactid[0]))
            #Paire qui existe
            if(ncomp!=None): self.classtype[i]=1
            else: self.classtype[i]=0
            if(ncomp and predict_RNA):                
                    #zone d'interaction predite
                    Pinteract=self.getInteraction(self.pos[i],self.getRangePosit(predict_RNA,0),self.getRangePosit(predict_RNA,1),self.getRangePosit(predict_RNA,2),self.getRangePosit(predict_RNA,3),True,intefunc)
                    #zone d'interaction reelle
                    Rinteract=self.getInteraction(self.pos[i],self.exp_inf.sRNA_deb[ncomp],self.exp_inf.sRNA_fin[ncomp],self.exp_inf.mRNA_deb[ncomp],self.exp_inf.mRNA_fin[ncomp],False,intefunc)
                    #Calcul de la precision
                    self.ppv[i]=self.compute_ppv(Pinteract, Rinteract)
                    #Calcul de la sensibilite
                    self.sensitivity[i]=self.compute_sensitivity(Pinteract,Rinteract)   
            else:                
                self.ppv[i]=np.nan
                self.sensitivity[i]=np.nan
            i+=1
        #Communication object
        data=Communication()
        data.classtype=self.classtype
        data.ppv=self.ppv
        data.sensitivity=self.sensitivity
        return data
    
    def Interaction_data(self):
        """
        Get interaction position data
        """
        self.type_sol=self.type_sol
        intefunc=self.defintefunction()
        interact=[]
        i=0
        for interactid in self.interactions:
            if((i%500)==0): print("%d/%d"%(i+1,len(self.interactions)))
            #Position d'interaction sRNA mRNA
            predict_RNA=self.dbmanage.getPositionsbyIntid(interactid[0])
            #Paire qui existe
            if(predict_RNA):
                sRNA_deb=self.getRangePosit(predict_RNA,0)
                sRNA_end=self.getRangePosit(predict_RNA,1)
                mRNA_deb=self.getRangePosit(predict_RNA,2)
                mRNA_end=self.getRangePosit(predict_RNA,3)
                #zone d'interaction predite
                Pinteract=self.getInteraction(self.pos[i],sRNA_deb,sRNA_end,mRNA_deb,mRNA_end,True,intefunc)
                #simplication de la position
                if(self.type_sol==0 or self.type_sol==1): temp=[sRNA_deb[self.pos[i]],sRNA_end[self.pos[i]],mRNA_deb[self.pos[i]],mRNA_end[self.pos[i]]]
                elif(self.type_sol==2 or self.type_sol==3): temp=[min(sRNA_deb),max(sRNA_end),min(mRNA_deb),max(mRNA_end)]                    
                else: temp=[sRNA_deb,sRNA_end,mRNA_deb,mRNA_end]
                interact+=[temp+[len(Pinteract[0]),len(Pinteract[1])]]
            else: interact+=[None]
            i+=1
        #Communication object
        data=Communication()
        data.interact=interact
        return data
    
    def run(self):
        """
        Compute interaction precision and sensitivity
        """
        if(self.exp_inf): return self.compute_sens_ppv()
        else: return self.Interaction_data()
    
  