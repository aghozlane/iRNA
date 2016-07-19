"""
 @brief: Manage Mpi run
 @author: Amine Ghozlane
 @version: 1.0 
"""
# -*- coding: utf-8 -*-
import time,math,sys, numpy as np
from mpi4py import MPI
from NormScore import *
from Interaction import *
from pValue import *
from Regression import *
from Communication import *

class Mpi:
    """
    @brief: Manage Mpi run 
    """
    def __init__(self):
        """
        The Constructor
        @note comm: Broadcast communicator
        @note nprocs: Number of thread
        @note myrank: Rank of the thread
        """
        self.comm = MPI.COMM_WORLD
        self.nprocs = self.comm.Get_size()
        self.myrank = self.comm.Get_rank()
        self.name = MPI.Get_processor_name()
        
    def waiting(self):
        """
        Test de reduction de consommation de ressource
        """
        while not self.comm.Iprobe(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG):
            time.sleep(0.1)
    
    def instanciate_analysis(self, data_obj,dbmanage,exp_inf,rand_inf):
        """
        Instanciate object of analysis
        @param data_obj: Comunication object
        @param dbmanage: Access to the database
        @param exp_inf: Experimental data
        @param rand_inf: Random data
        """
        if(data_obj.type_analysis==0):   return NormScore(data_obj,dbmanage)
        elif(data_obj.type_analysis==1): return Interaction(data_obj,dbmanage,exp_inf)
        elif(data_obj.type_analysis==2): return pValue(data_obj,dbmanage,rand_inf)
        elif(data_obj.type_analysis==3): return Interaction(data_obj,dbmanage,None)
        elif(data_obj.type_analysis==4): return Regression(data_obj)
        else: raise ValueError("type_analysis unknown = %d"%data_obj.type_analysis)
                
    def worker(self,dbmanage,exp_inf,rand_inf):
        """
        Worker process
        @param dbmanage: Access to the database
        @param exp_inf: Experimental data
        @param rand_inf: Random data
        """
        status = MPI.Status()
        while(status.tag!=0):
            data_obj=None
            #Get data
            self.waiting()
            data_obj=self.comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG,status=status)
            if(status.tag!=0):           
                #Intanciate analysis
                if(data_obj):
                    analysis_obj=self.instanciate_analysis(data_obj,dbmanage,exp_inf,rand_inf)
                    result=analysis_obj.run()
                    self.comm.send(result, dest=status.source, tag=1)
                else:
                    raise ValueError("Worker received a data_obj empty")
        #Finalisation    
        self.Mpi_finalize()
        
    def working_step(self,nbinteract):
        """
        Defines what each worker will compute
        @param nbinteract: Number of operation to do
        """
        tab=[]
        deb=0
        #Number of process
        nbproc=self.nprocs-1
        if(nbproc==0): return(tab)
        val=nbinteract/nbproc
        #More process than interaction
        if(nbproc>nbinteract):
            val=1
            nbproc=nbinteract
        #Defining repartion
        end=val
        for i in xrange(nbproc):
            tab+=[[deb,end]]
            deb+=val
            end+=val
        #Impair case
        if(tab[-1][-1]!=nbinteract):
            tab[-1][-1]+=nbinteract-tab[-1][-1]
        #Rest of proc
        #for j in xrange(i+1,(self.nprocs-1)):
        #    tab+=[0]
        return(tab) 
    
    def getCorrespondingSoft(self,computer,dbmanage,rand_inf):
        """
        Get corresponding soft between soft_param and current database
        """
        idsofts=dbmanage.getidSoftsbyname(computer.name)
        idpossible_softs=rand_inf.getallSofts(computer.name)
        softname=None
        #Search for corresponding
        if(len(idsofts)>1 and (len(idsofts)==len(idpossible_softs))):
            val=np.where(idsofts==computer.softid)
            softname=computer.name+"_"+str(idpossible_softs[val[0]])
        #Only one soft for each side
        elif(len(idsofts)==1 and (len(idsofts)==len(idpossible_softs))):
            softname=computer.name+"_"+str(idpossible_softs[0])
        #Problem nothing to do
        else:
            softname=computer.name+"_"+str(computer.softid)
            #sys.exit("There is not same number of %s soft for each database. Impossible to get corresponding soft"%dbmanage.getSoftname(computer.softid))
        return softname
       
    def getSendCom(self,computer,type_analysis,tab,dbmanage,rand_inf):
        """
        Create communication object
        """
        data=Communication()
        data.type_analysis=type_analysis
        if(type_analysis==0):
            data.interactions=computer.interactions[tab[0]:tab[1]]
            data.type_sol=computer.type_sol
            data.score_type=computer.score_type
        elif(type_analysis==1 or type_analysis==3):
            data.interactions=computer.interactions[tab[0]:tab[1]]
            data.type_sol=computer.type_sol
            #Matrice partielle
            data.pos=computer.pos[tab[0]:tab[1]]
        elif(type_analysis==2):
            data.softname=self.getCorrespondingSoft(computer,dbmanage,rand_inf)
            data.norm_score=computer.norm_score[tab[0]:tab[1]]
            data.interactions=computer.interactions[tab[0]:tab[1]]
        #regression
        elif(type_analysis==4):
            data.unique_sRNAidinint=computer.unique_sRNAidinint[tab[0]:tab[1]]
            data.norm_score=computer.norm_score
            data.sRNAid_tab=computer.sRNAid_tab
            data.pValue_type=computer.pValue_type
        else: 
            raise NotImplementedError("")
        return data
    
    def Mpi_end_processus(self):
        """
        End of all process
        """
        for i in xrange(1, self.nprocs):
            self.comm.send(i,dest=i,tag=0)
    
    def Mpi_finalize(self):
        """
        Fin
        """
        #Synchronisation des threads
        self.comm.Barrier()
        #Finalisation
        MPI.Finalize()
    
    def getNormScore(self,computer, tab, source, result):
        """
        Get results of score normalisation 
        """
        computer.norm_score[tab[source][0]:tab[source][1]]=result.norm_score
        computer.pos[tab[source][0]:tab[source][1]]=result.pos
         
    def getInteraction(self, computer, tab, source, result):
        """
        Get results of the interaction analysis
        """
        computer.classtype[tab[source][0]:tab[source][1]]=result.classtype
        computer.ppv[tab[source][0]:tab[source][1]]=result.ppv
        computer.sensitivity[tab[source][0]:tab[source][1]]=result.sensitivity
    
    def getInteraction_data(self, computer, tab, source, result):
        """
        Get results of the interaction analysis
        """
        computer.interact[tab[source][0]:tab[source][1]]=result.interact
    
    def getpValue(self, computer, tab, source, result):
        """
        Get results of pValue analysis
        """
        computer.pValue[tab[source][0]:tab[source][1]]=result.pValue
    
    def getRegression(self, computer, tab, source, result):
        """
        Get results of linear regression
        """
        computer.cumneg[tab[source][0]:tab[source][1]]=result.cumneg
        computer.curve_param[tab[source][0]:tab[source][1]]=result.curve_param
    
    def getAnalysis(self,type_analysis):
        """
        Get the analysis function to apply
        """
        if(type_analysis==0): self.setdata=self.getNormScore
        elif(type_analysis==1): self.setdata=self.getInteraction
        elif(type_analysis==2): self.setdata=self.getpValue
        elif(type_analysis==3): self.setdata=self.getInteraction_data
        elif(type_analysis==4): self.setdata=self.getRegression
        else : raise ValueError("type_analysis unknown = %d"%type_analysis)
        
    def run(self,computer,type_analysis,dbmanage=None,rand_inf=None):
        """
        Manage comunication and data sending
        @param computer: Computer object
        @param type_analysis: Type of analysis
        """
        status = MPI.Status()
        #compute table working step
        if(type_analysis!=4): tab=self.working_step(computer.nbinteract)
        else: tab=self.working_step(len(computer.unique_sRNAidinint))
        #SetAnalysis expected
        self.getAnalysis(type_analysis)
        #Sending work
        for i in xrange(len(tab)):
            if(tab[i]!=0): self.comm.send(self.getSendCom(computer,type_analysis,tab[i],dbmanage,rand_inf), dest=(i+1), tag=1)
        #Receiving result
        j=0
        while j<(len(tab)):
            self.waiting()
            result=self.comm.recv(source=MPI.ANY_SOURCE, tag=1,status=status)
            #computer.setdata(tab, status.source-1, result)
            self.setdata(computer,tab, status.source-1, result)
            j+=1
            print("proc %d/%d"%(j,len(tab)))