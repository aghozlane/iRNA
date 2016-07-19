# -*- coding: utf-8 -*-
"""
 @brief: Get sRNA information from sqlite database
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv,os
from Parser import *
#sys.path.append(os.path.join(os.path.dirname(__file__),"..%siRNA_pred%s"%(os.sep,os.sep)))
from Sqlite_manager import *

class Database_data(Parser):
    
    def __init__(self, iRNA_db,fastmode):
        """
        Instanciate pValue object
        @param pValue_file: pValue file
        """
        Parser.__init__(self)
        #Create db object
        self.dbmanage=Sqlite_manager(None,iRNA_db,fastmode)
        
    def setRNAlength(self,type_RNA,node_objects):
        """
        Set RNA length
        @param type_RNA: Type of RNA
        @param node_objects: List of node object
        """
        for name,length in self.dbmanage.getallRNAlength(type_RNA):
            node=self.getlinknode(node_objects,name)
            if(node): node.length=length           
        return node_objects
    
    def setdata(self,node_objects,edge_objects):
        """
        Add Length information to nodes
        @param node_objects: List of node objects
        @param edge_objects: List of edge objects
        """
        #Connect DB
        self.dbmanage.connectDB()
        #setsRNA
        node_objects=self.setRNAlength(0,node_objects)
        #setmRNA
        node_objects=self.setRNAlength(1,node_objects)
        #Disconnect DB
        self.dbmanage.disconnectDB()
        return(node_objects,edge_objects)