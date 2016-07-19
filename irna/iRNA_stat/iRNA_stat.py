#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 @brief: Create multilist file
 @author: Amine Ghozlane
 @version: 1.0 
"""
import os, sys
#, argparse
from Files import *
from Computer import *
from Soft import *
from Exp_data import *
from Rand_data import *
from Draw_data import *
from Mpi import *
from Pickling import *
from Execute import *
from Threshold import *
from pValue_selection import *
from Interaction import *
#sys.path.append(os.path.join(os.path.dirname(__file__),"..%siRNA_pred%s"%(os.sep,os.sep)))
from Sqlite_manager import *

def usage(parser,args):
    """
    Test correct usage of arguments
    @param parser: Parser object
    @param args: Arguments
    """
    #Run without arguments
    if len(sys.argv)== 1:
        parser.print_usage()
        sys.exit()
    if(not os.path.isfile(args.iRNA_db)):
        print("Error with \"%s\" : --iRNA_db required a file\n"%args.iRNA_db)
        parser.print_help()
        sys.exit()
            
##Determine les fichiers fournis en arguments
#def getArgument():
#    """
#    Determine the argument
#    @return: arguments
#    """
#    #Parsing arguments
#    parser = argparse.ArgumentParser(description='Statistical analysis of RNA-RNA predict interaction.')
#    parser.add_argument('-d', '--iRNA_db',help='Path to iRNA db',required=True)
#    parser.add_argument('-i', '--soft_inf',help='Path to soft information file',required=True)
#    parser.add_argument('-e', '--exp_inf',help='Path to experimental data file')
#    parser.add_argument('-n', '--rand_inf',help='Path to random data file') 
#    parser.add_argument('-a', '--random',help='Random analysis',action='store_true')         
#    parser.add_argument('-r', '--results',help='Path to result repertory',required=True)
#    args = parser.parse_args()
#    
#    #Verify usage
#    usage(parser,args)
#    return args

def main():
    """
    Main program function
    """
    exp_inf=None
    rand_inf=None
    #Initiation de MPI
    parallel=Mpi()
    #Get the arguments
    #args=getArgument()   
    args=Files(parallel.myrank)
    #Create db object
    dbmanage=Sqlite_manager(None,args.iRNA_db,args.fastmode)
    #Connect DB
    dbmanage.connectDB()
    #Read random data file
    if(args.rand_inf):
        rand_inf=Rand_data(args.rand_inf)
    # Main process
    if(parallel.myrank==0):
        execution_time=None
        #Data represent
        data_represent=draw_data(args.results)
        #data_writing=writer(args.results)
        savemethod=Pickling(args.save)
        #Load exec_inf
        if(args.exec_inf):
            execution_time=Execute(args.exec_inf)
        #Load threshold
        if(args.thres_inf):                
            pValueselect=pValue_selection(Threshold(args.thres_inf)) 
        #inte=Interaction(None,dbmanage,None)
        #Soft data
        soft_inf=Soft(args.soft_inf)        
        #Computer table
        mycomputer=[]            
        for softid,name in dbmanage.getAllSoft():
            obj=None
            print("soft %s %d"%(name,softid))
            #Load saved object
            if(args.save and not args.overwrite):
                obj=savemethod.loadobj(name,softid)
                if(obj): mycomputer+=[obj]
            if not obj:
                #Get soft information
                numsoft=soft_inf.getSoftnum(name.lower())
                #Computer
                mycomputer+=[Computer(softid,name,numsoft,soft_inf.type_sol[numsoft],soft_inf.score_type[numsoft],dbmanage,args)]
                #Compute normalized score
                parallel.run(mycomputer[-1],0)
                #compute interaction
                if(args.exp_inf): parallel.run(mycomputer[-1],1)
                #Compute regression
                if(args.random and args.pValue): parallel.run(mycomputer[-1],4)
                #Compute pValue          
                elif(not args.random and args.rand_inf):
                    #if(name=="ssearch35_t" or name=="yass-Linux64.bin" or name=="blastall"):
                    if(parallel.getCorrespondingSoft(mycomputer[-1],dbmanage,rand_inf).lower() not in rand_inf.soft): mycomputer[-1].pValue=mycomputer[-1].norm_score
                    else:  parallel.run(mycomputer[-1],2,dbmanage,rand_inf)
                    #compute interaction
                    parallel.run(mycomputer[-1],3) 
                else:
                    if(mycomputer[-1].score_type==2 or mycomputer[-1].score_type==3): mycomputer[-1].pValue=-mycomputer[-1].norm_score
                    else: mycomputer[-1].pValue=mycomputer[-1].norm_score
                    #compute interaction
                    parallel.run(mycomputer[-1],3)
                #Compute pValue selection
                if(args.thres_inf): pValueselect.run(mycomputer[-1],dbmanage,name,softid) 
                #Save result
                if(args.save): savemethod.saveobj(mycomputer[-1],name,softid)
        #Finalization of communication
        if((parallel.nprocs-1)>0):
            parallel.Mpi_end_processus()
            parallel.Mpi_finalize()                               
        #Plot
        data_represent.plot(args,soft_inf,mycomputer,dbmanage,execution_time)
    #Worker processes       
    else:
        #Read experimental data file  
        if(args.exp_inf):
            exp_inf=Exp_data(args.exp_inf)
        #Run work
        parallel.worker(dbmanage,exp_inf,rand_inf) 
    #Disconnect DB
    dbmanage.disconnectDB()
                 

if __name__ == "__main__":
    main()
    