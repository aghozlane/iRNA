"""
 @brief: Handle similarity data
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, re,os,csv
from Edge import *
from Parser import *

class Similarity(Parser):
    
    def __init__(self, similarity_file):
        """
        Instanciate similarity parser object
        @param similarity_file: Similarity file
        """
        Parser.__init__(self)
        self.similarity_file=similarity_file
    
    def setdata(self,node_objects,edge_objects):
        """
        Add sRNA - sRNA edges based on their similarity
        @param node_objects: list of node objects
        @param edge_objects: list of edge objects
        """
        try:
            similarityReader = csv.reader(open(self.similarity_file, 'rb'), delimiter='\t')
            #passage de l'entete
            similarityReader.next()
            for i in similarityReader:
                node1=self.getlinknode(node_objects,i[0])
                node2=self.getlinknode(node_objects,i[1])
                #Create edge
                if(node1 !=None and node2 !=None):
                    edge_element=Edge(node1,node2)
                    #set similarity value and scategory
                    edge_element.similarity=float(i[2])
                    edge_element.category=2
                    #Add edge element in the list
                    self.addelements(edge_objects,edge_element)
                    #edge_objects+=[edge_element]
        except IOError:
            sys.exit("Error : can not open file %s"%self.similarity_file)
        #except:
        #    sys.exit("Something went wrong with %s"%self.similarity_file)
        return(node_objects,edge_objects)