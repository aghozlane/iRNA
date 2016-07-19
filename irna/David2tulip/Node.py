"""
 @brief: Handle RNA information
 @author: Amine Ghozlane
 @version: 1.0 
"""
class Node:
    
    def __init__(self,name,gender,num):
        """
        Instanciate Node object
        @param name: Name of the node
        @param gender: Gender of the node
        @param num: Num of the node
        """
        self.name=name
        self.code=name
        self.gender=gender
        self.num=num
        self.group=[]
        self.length=0
        #Set the group of the object
        if(gender==1): self.group+=[num]
    
    def addgroupelement(self,numgroup):
        """
        Add a group element
        @param numgroup: Number of a group
        """
        self.group+=[numgroup]
        
    def setuniquegroup(self):
        """
        Set unique element in the group
        """
        self.group={}.fromkeys(self.group).keys()
        
    def printgroupdata(self):
        """
        Print group data
        """
        #Get unique value of group
        self.setuniquegroup()
        #Print group data
        return (','.join(str(n) for n in self.group))
        
    def println(self):
        """
        Print object values
        """
        return ("%d\t%d\t%s\t%s\t%d\t(%s)\n"%(self.num,self.gender,self.code,self.name,self.length,self.printgroupdata()))
    
    def __cmp__(self,other):
        """
        General method to compare node based on the name
        @param other: Compared value
        """
        if self.name<other: return -1
        elif self.name==other: return 0
        elif self.name>other: return 1
        else: raise ValueError