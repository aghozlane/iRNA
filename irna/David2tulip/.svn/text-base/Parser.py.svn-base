# -*- coding: utf-8 -*-
"""
 @brief: General operation on graph
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, bisect, cPickle
from Node import *

class Parser:
    
    def __init__(self):
        """
        Instanciate Parser object
        """
        pass
    
    def createnodes(self,list_element,gender,number):
        """
        Create the list of nodes
        @param list_element: 
        @param gender: Gender of the list
        @return: list of nodes 
        """
        a=0
        listnodes=[]
        #Creation of nodes
        for i in list_element:
            listnodes+=[Node(i,gender,a+number)]
            a+=1
        #Sorted list
        listnodes.sort()
        return(listnodes)
    
    def addelements(self,list_objects,element):
        """
        Add an element in list of objects
        @param list_objects: list of objects
        @param element: object
        """
        #Add element with bisect method
        bisect.insort(list_objects,element)
   
    def getlinknode(self,node_objects, name):
        """
        Get a node based on its name
        @param node_objects: List of node objects
        @param name: Name of searched object
        @return: node
        """
        #Searching the node with its name
        i=bisect.bisect_left(node_objects,name)
        #Object has been found
        if(i!=len(node_objects) and node_objects[i].name==name): return node_objects[i]
        return None
    
    def getlinkedge(self,edge_objects,sRNA,mRNA):
        """
        @param edge_objects: List of edge objects
        @param sRNA: sRNA
        @param mRNA: mRNA
        """
        #Searching the node with its id
        i=bisect.bisect_left(edge_objects,sRNA+mRNA)
        #Object has been found
        if(i!=len(edge_objects) and edge_objects[i].node1.name==sRNA and edge_objects[i].node2.name==mRNA):
            return edge_objects[i]
        return None
    
    def getsRNA(self, node_objects):
        """
        Get list of sRNA
        @param node_objects: List of node objects
        @return: List of sRNA
        """
        return [i.name for i in node_objects if(i.gender==1)]
    
    def getTarget(self,sRNA,edge_objects):
        """
        Get targets of one sRNA if category
        @param sRNA: sRNA name
        @param edge_objects: List of edge objects
        """
        return [i.node2.name for i in edge_objects if(i.node1.name==sRNA and i.category==1)]
    
    def loadObject(self,pickledump):
        """
        Load object
        @param pickledump: Filename
        @return: Object loaded
        """
        try:
            mypickle=open(pickledump,"r")
            obj=cPickle.load(mypickle)
            mypickle.close()
        except IOError:
            sys.exit("Error : can not open file %s"%pickledump)
        #except:
        #    sys.exit("Something went wrong with %s"%results)
        return obj
    
    def printObject(self,obj,pickledump):
        """
        Write object
        @param obj: Object to dump
        @param pickledump: Filename
        """
        try:
            mypickle=open(pickledump,"w")
            cPickle.dump(obj,mypickle)
            mypickle.close()
        except IOError:
            sys.exit("Error : can not open file %s"%pickledump)
        #except:
        #    sys.exit("Something went wrong with %s"%results)  
        