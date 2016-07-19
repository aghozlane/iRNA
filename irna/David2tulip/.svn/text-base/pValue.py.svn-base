# -*- coding: utf-8 -*-
"""
 @brief: Handle pValue information
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv
from Parser import *

class pValue(Parser):
    
    def __init__(self, pValue_file):
        """
        Instanciate pValue object
        @param pValue_file: pValue file
        """
        Parser.__init__(self)
        self.pValue_file=pValue_file
        
    def setdata(self,node_objects,edge_objects):
        """
        Add pValue information to edges
        @param node_objects: List of node objects
        @param edge_objects: List of edge objects
        """
        try:
            pValueReader = csv.reader(open(self.pValue_file, 'rb'), delimiter='\t')
            #passage de l'entete
            pValueReader.next()
            for i in pValueReader:
                if(int(i[3])==1):
                    edge_element=self.getlinkedge(edge_objects,i[0],i[1])
                    if(edge_element!=None): edge_element.pValue=float(i[2])
        except IOError: sys.exit("Error : can not open file %s"%self.pValue_file)
        except: sys.exit("Something went wrong with %s"%self.pValue_file)
        return(node_objects,edge_objects)