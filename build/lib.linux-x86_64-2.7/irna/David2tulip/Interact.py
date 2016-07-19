"""
 @brief: Handle interaction data
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv
from Parser import *

class Interact(Parser):
    
    def __init__(self, interact_file):
        """
        Instanciate similarity parser object
        @param interact_file: interact file
        """
        Parser.__init__(self)
        self.interact_file=interact_file
    
    def setdata(self,node_objects,edge_objects):
        """
        Add sRNA - sRNA edges based on their similarity
        @param node_objects: list of node objects
        @param edge_objects: list of edge objects
        """
        try:
            interactReader = csv.reader(open(self.interact_file, 'rb'), delimiter='\t')
            #passage de l'entete
            interactReader.next()
            for i in interactReader:
                if(int(i[10])==1):
                    edge_element=self.getlinkedge(edge_objects,i[0],i[1])
                    if(edge_element!=None):
                        edge_element.sRNA_deb=int(i[2])
                        edge_element.sRNA_end=int(i[3])
                        edge_element.mRNA_deb=int(i[4])
                        edge_element.mRNA_end=int(i[5])
        except IOError:
            sys.exit("Error : can not open file %s"%self.interact_file)
        except:
            sys.exit("Something went wrong with %s"%self.interact_file)
        return(node_objects,edge_objects)