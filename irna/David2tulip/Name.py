"""
 @brief: Handle corresponding name information
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, re,os,csv
from Edge import *
from Parser import *

class Name(Parser):
    
    def __init__(self, name_file):
        """
        Instanciate name parser object
        @param name_file: Name file
        """
        Parser.__init__(self)
        self.name_file=name_file
    
    def setdata(self,node_objects,edge_objects):
        """
        Add sRNA - sRNA edges based on their similarity
        @param node_objects: list of node objects
        @param edge_objects: list of edge objects
        """
        try:
            nameReader = csv.reader(open(self.name_file, 'rb'), delimiter='\t')
            #passage de l'entete
            nameReader.next()
            for i in nameReader:
                node=self.getlinknode(node_objects,i[0].strip().upper())
                #Modify node informations
                if(node!=None): node.code=i[1].strip()
        except IOError:
            sys.exit("Error : can not open file %s"%self.name_file)
        #except:
        #    sys.exit("Something went wrong with %s"%self.similarity_file)
        return(node_objects,edge_objects)