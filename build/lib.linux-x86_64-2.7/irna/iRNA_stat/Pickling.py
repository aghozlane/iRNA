# -*- coding: utf-8 -*-
"""
 @brief: Handle DAVID information pickling
 @author: Amine Ghozlane
 @version: 1.0 
"""
import cPickle,sys,os

class Pickling:
    
    def __init__(self,path):
        """
        Instanciate Pickling object
        @param path: Path of pickling object
        """
        self.path=path
        
    def loadobj(self,name,num):
        """
        Load the pickling object of the software result 
        @param name: Name of the software
        @param num: Software id
        """
        obj=None
        pickle=self.path+name+"_"+str(num)+".pkl"
        #file does not exist
        if not os.path.isfile(pickle): return None
        #Open file
        try:
            mypickle=open(pickle,"r" )
            obj=cPickle.load(mypickle)
            mypickle.close()
        except IOError: sys.exit("Error : can not open file %s"%pickle)
        return obj
    
    def saveobj(self,obj,name,num):
        """
        Pickle the object
        @param obj: Computed object
        @param name: Name of the software
        @param num: Software id
        """
        pickle=self.path+name+"_"+str(num)+".pkl"
        try:
            mypickle=open(pickle,"w" )
            cPickle.dump(obj,mypickle)
            mypickle.close()
        except IOError: sys.exit("Error : can not open file %s"%pickle)