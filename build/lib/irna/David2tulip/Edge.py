"""
 @brief: Handle interaction
 @author: Amine Ghozlane
 @version: 1.0 
"""
class Edge:

    def __init__(self,node1,node2):
        """
        Instanciate Edge object
        @param node1: sRNA node
        @param node2: mRNA node
        """
        self.node1=node1
        self.node2=node2
        self.pValue=0.0
        self.similarity=0.0
        self.category=1
        self.categoryName=[]
        self.termName=[]
        self.sRNA_deb=0
        self.sRNA_end=0
        self.mRNA_deb=0
        self.mRNA_end=0
        
    def println(self):
        """
        Print edge value
        @return: String converted values 
        """
        self.dbrec=len({}.fromkeys(self.categoryName).keys())
        result="%d\t%d\t%d\t%f\t%f\t%f\t(%d,%d)\t(%d,%d)\t%d\n"%(self.node1.num,self.node2.num,self.category,self.similarity,1.0-self.similarity,self.pValue,self.sRNA_deb,self.sRNA_end,self.mRNA_deb,self.mRNA_end,self.dbrec)
        return(result)
    
    def __cmp__(self,other):
        """
        General method to compare edge based on the name of nodes
        @param other: Compared value
        """
        temp1=self.node1.name+self.node2.name
        if(temp1<other): return -1
        elif(temp1==other): return 0
        elif(temp1>other): return 1
        else: raise ValueError