#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 @brief: Predict interaction
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, os, time
from Files import *
from Fasta import *
from Mpi import *
from RanRNA import *
from Comparison import *
from Sqlite_manager import *
from Merger import *

def WriteAllExecTime(result_out, start, end):
    """
    Write execution time
    @param result_out: Time result
    @param start: Date of start
    @param end: Date of end
    """
    out=result_out+"execution_time.txt"
    try:
        exectime=open(out,"a")
        exectime.write("iRNA\t%f\n"%(end-start))
        exectime.close()
    except IOError: sys.exit("Error can not open file %s"%out)
    except : sys.exit("Something went wrong with %s"%out)

def saveData(data_files, fastmode):
    """
    Save data in a sqlite database
    """
    #Create db object
    dbmanage=Sqlite_manager(data_files.predict[1],None,fastmode)
    #Merger
    m=Merger(data_files.predict[1])
    #Delete previous database
    if(os.path.isfile(dbmanage.db_path)): os.remove(dbmanage.db_path)
    #Connect DB
    dbmanage.connectDB()
    #Create database
    dbmanage.createSQLdb()
    #set sRNA
    m.setRNA(dbmanage,data_files.fasta[2]+"seq_inf.txt",0)
    #set mRNA
    m.setRNA(dbmanage,data_files.fasta[3]+"seq_inf.txt",1)        
    #Merge trees in one database
    m.merge(dbmanage)
    #Create indexes
    dbmanage.createIndexes()
    #Disconnect database
    dbmanage.disconnectDB2()
        
def main():
    """
    @brief: Main program
    @author: Amine Ghozlane
    """    
    #Initiation de MPI
    parallel=Mpi()
    #Recuperation des fichiers
    data_files=Files(parallel.Mpi_getmyrank())
    #Creation de l'objet comparaison
    compair = Comparaison(data_files.predict[2],  data_files.fasta[2], data_files.fasta[3], data_files.predict[4], data_files.fasta[1])
   
    #Probleme de connexion
    print("I am process %d of %d on %s"%(parallel.Mpi_getmyrank(),parallel.Mpi_getnprocs(), parallel.Mpi_getname()))
    if(parallel.Mpi_getnprocs()==1 and data_files.operation[1]):  sys.exit("Rank problem, nprocs=%d"%parallel.Mpi_getnprocs())
    
    #Extraction du multifasta
    if(data_files.operation[0] and parallel.Mpi_getmyrank()==0):
        fasta_files=Fasta()
        print("Start extraction...")
        #MF sRNA
        print("sRNA:")
        fasta_files.ExtractFasta(data_files.fasta[0],data_files.fasta[2])
        #MF mRNA
        print("mRNA:")
        (GC, lenran)=fasta_files.ExtractFasta(data_files.fasta[1],data_files.fasta[3])
        print("Extraction done")
        
    #Generate random sequence
    if(data_files.operation[2] and parallel.Mpi_getmyrank()==0):
        generator=None
        #If use extract option is on
        if(data_files.operation[3]): generator=RanRNA(GC, lenran, data_files.random[2], data_files.random[3])
        #Else use given data
        else: generator=RanRNA(data_files.random[0], data_files.random[1], data_files.random[2], data_files.random[3])
        generator.GenerateFile()
    
    #Generation des couples
    if(data_files.operation[4] and parallel.Mpi_getmyrank()==0):
        #Creation des couples
        compair.prepareComparaison() 
    
    #Execution des programmes   
    if(data_files.operation[1]):
        #Debut du decompte
        a=time.gmtime()
        #Lancement des calculs
        parallel.Mpi_run(data_files.predict, compair)
        #Debut du decompte
        b=time.gmtime()
        #Calcul du temps d'execution
        #Bug en cas de changement de mois ou d'annee entre les deux
        start=a.tm_sec+a.tm_min*60+a.tm_hour*3600+a.tm_mday*24*3600
        end=b.tm_sec+b.tm_min*60+b.tm_hour*3600+b.tm_mday*24*3600
        #Ecriture du temps d'execution global
        WriteAllExecTime(data_files.predict[1], start, end)
        #Merge des fichers
        if(parallel.Mpi_getmyrank()==0): saveData(data_files,data_files.fastmode)

if __name__ == "__main__":
    main()