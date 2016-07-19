#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
 @brief: Parsing program
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),"..%sirna%siRNA_pred%s"%(os.sep,os.sep,os.sep)))
from Sqlite_manager import *
from Merger import *

#Determine les fichiers fournis en arguments
def getArgument():
    """
    Determine the argument
    @return: argument
    """
    if len(sys.argv)!= 4:
        sys.exit("Usage : testMerge.py sRNA_inf mRNA_inf repIn")
    return sys.argv[1],sys.argv[2],sys.argv[3]

def main():
    """
    Main function
    """
    sRNA_inf,mRNA_inf,repIn =getArgument()
    #Create db object
    dbmanage=Sqlite_manager(repIn,None,True)
    #Delete previous database
    if(os.path.isfile(dbmanage.db_path)):
        print("remove %s"%dbmanage.db_path)
        os.remove(dbmanage.db_path)
    #Merger
    m=Merger(repIn)
    #Connect DB
    dbmanage.connectDB()
    #Create database
    dbmanage.createSQLdb()
    #set sRNA
    m.setRNA(dbmanage,sRNA_inf,0)
    #set mRNA
    m.setRNA(dbmanage,mRNA_inf,1)        
    #Merge trees in one database
    m.merge(dbmanage)
    #Create indexes
    dbmanage.createIndexes()
    #Disconnect database
    dbmanage.disconnectDB2()
    
if __name__ == "__main__":
    main()