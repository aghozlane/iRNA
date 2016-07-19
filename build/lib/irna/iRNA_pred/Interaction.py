"""
 @brief: Save interactions into an xml tree
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, re
from lxml import etree

class Interaction:
    """
    @brief: Store information on interaction.
    """
    matrix=None
    
    def __init__(self):
        """
        Init Interaction object
        """
        self.root= etree.Element("result")
        
    def setSoft(self, soft_name, ref):
        """
        Set software node
        @param soft_name: Software name
        @param ref: Software num        
        """
        self.schema_root=etree.SubElement(self.root,"soft", attrib={"reference":str(ref),"id":soft_name})
            
    def getlastmRNA(self):
        """
        Get mRNA attribute
        @return: Value of mRNA attribute
        """
        if(list(self.schema_root)): return self.getLastchildnode().attrib.get("mRNA")
        return ""
    
    def getlastsRNA(self):
        """
        Get sRNA attribute
        @return: Value of sRNA attribute
        """
        if(list(self.schema_root)): return self.getLastchildnode().attrib.get("sRNA")
        return ""
        
    def addListOfComparison(self, mRNA, sRNA):
        """
        Add a list of comparison node
        @param mRNA: mRNA id
        @param sRNA: sRNA id
        """
        self.schema_root.append(etree.Element("listOfComparison", attrib={"mRNA":mRNA, "sRNA":sRNA}))
        
    def addComparison(self, data):
        """
        Add a comparison node
        @param data: Comparison result 
        """
        #debut mRNA fin mRNA debut sRNA fin sRNA local_score
        self.getLastchildnode().append(etree.Element("comparison", attrib={"mRNA_begin":str(data[0]), "mRNA_end":str(data[1]), "sRNA_begin":str(data[2]), "sRNA_end":str(data[3]), "score":str(data[4])}))
        
    def getLastchildnode(self):
        """
        @return: The last child from root
        """
        return list(self.schema_root)[-1]
    
    def setScoreList(self, score):
        """
        Set ListOfComparison score
        @param score: Score the whole comparison
        """
        self.getLastchildnode().set("score",str(score))        
        
    def getResult(self):
        """
        @return: The tree in text format
        """
        #return etree.tostring(self.schema_root)
        return(etree.tostring(self.root, method="xml", pretty_print=True))
              
    def println(self):
        """
        Print object data
        """
        print(etree.tostring(self.schema_root, pretty_print=True))        
        
    #def getMatrix(self, matrix):
    def setMatrix(self, matrix):
        """
        Manage matrix load and format
        @param matrix: Matrix file
        """
        #Get matrix
        mat = self.loadMatrix(matrix)
        #Format matrix
        self.matrix = self.formatMatrix(mat)
        #return mat
    
    def formatMatrix(self, mat):
        """
        Format into integer the matrix
        @param mat: Matrix data in string format
        @return: Matrix in integer format
        """
        try:
            #Convert into integers data
            for i in range(1,(len(mat[0])+1)):
                for j in range(1,(len(mat[0])+1)):
                    mat[i][j]=int(mat[i][j])
        except: sys.exit("Something went wrong in formatMatrix with %s"%mat[i][j])
        return mat
    
    def loadMatrix(self, matrix):
        """
        Load the matrix of score
        @param matrix: Matrix file
        @return: Matrix of score
        """
        regex_matrix=re.compile("^(?!#)\S*\s+\S*")
        mat=[]
        try:
            #Open file
            mat_file=open(matrix,"r")
            for i in mat_file:
                #Detect matrix data
                a=regex_matrix.match(i)
                #Load data into matrix
                if a: mat.append(i.rsplit())
            #Close file
            mat_file.close()
        #Read error
        except IOError: sys.exit("Error : Can not open %s"%matrix)
        #Other error
        except : sys.exit("Something went wrong with %s"%matrix)
        return mat
    
    def getPosition(self, letter):
        """
        Get the position on matrix
        @param letter: One letter in the sequence
        @return: The position on matrix 
        """
        count=1
        #Get the position on the matrix
        for i in self.matrix[0]:
            #Letter found
            if(i==letter): return count
            else: count+=1
        #Letter not found
        sys.exit("Letter not found : %s"%letter)
        return
    
    def scoreSequence(self, seqs):
        """
        Score two sequences
        @param seqs: Table of two sequences
        @return: The score of the alignement
        """
        score=0
        for i in range(len(seqs[0])):
            #Sum score                
            score+=self.matrix[self.getPosition(seqs[0][i])][self.getPosition(seqs[1][i])]
        return score