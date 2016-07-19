"""
 @brief: Select interaction depending on threshold
 @author: Amine Ghozlane
 @version: 1.0 
"""
import bisect,sys
from Parser import *

class pValue_selection(Parser):
    
    def __init__(self,pValue_thres):
        """
        Instanciate pValue_selection object
        @param pValue_thres: pValue_thres object
        """
        self.pValue_thres=pValue_thres
    
    def getnumsRNA(self,allsRNA,sRNA):
        """
        Get an sRNA in the list
        @param allsRNA: List of sRNAs 
        @param sRNA: An sRNA
        """
        #Searching the node with its name
        i=bisect.bisect_left(allsRNA,sRNA)
        #Object has been found
        if(i!=len(allsRNA) and allsRNA[i]==sRNA): return i
        else: sys.exit("sRNA %s not present in the database"%sRNA)
    
    def compute_Davidmatrix(self,computer,dbmanage,name,softid):
        """
        Compute multilist target
        @param computer: Computer object
        @param dbmanage: Access to the database
        @param name: Software name
        @param softid: Software id
        """
        thres=self.pValue_thres.threshold[self.pValue_thres.getSoftnum(name+"_"+str(softid))]
        #Get All sRNA
        computer.allsRNA=dbmanage.getallsRNAname()
        computer.allsRNA.sort()
        #Select mRNA under the pValue threshold
        selection=[]
        tabsRNA=[0]*len(computer.allsRNA)
        tabmRNA=[]
        for i in xrange(computer.nbinteract):
            if(computer.pValue[i]<=thres):
                computer.select[i]=True
                sRNA=dbmanage.getsRNAbyIntid(computer.interactions[i][0])
                mRNA=dbmanage.getmRNAbyIntid(computer.interactions[i][0])
                selection+=[{'sRNA': sRNA, 'mRNA': mRNA}]
                tabsRNA[self.getnumsRNA(computer.allsRNA,sRNA)]+=1
                tabmRNA+=[mRNA]
        #Unique mRNA
        unique_tabmRNA=self.getUnique(tabmRNA)
        #create david file
        matrix=[["" for j in range(0,len(computer.allsRNA))] for i in range(0,max(tabsRNA))]
        matrix_sRNA=[[None]]*len(computer.allsRNA)
        i=0
        for sRNA in computer.allsRNA:
            j=0
            tab=[]
            for select in selection:
                if(select['sRNA']==sRNA):
                    matrix[j][i]=select['mRNA']
                    tab+=[select['mRNA']]
                    j+=1
            matrix_sRNA[i]=tab
            i+=1
        return(matrix,matrix_sRNA,tabmRNA,unique_tabmRNA)
    
    def compute_frequency(self,computer,tabmRNA,unique_tabmRNA):
        """
        Compute frequency of mRNA
        @param computer: Computer object
        @param tabmRNA: Redondant list of mRNA
        @param unique_tabmRNA: Unique mRNA
        """
        frequencydict=[{"mRNA":unique_tabmRNA[i],"frequency":(float(tabmRNA.count(unique_tabmRNA[i]))/float(len(computer.allsRNA))*100.0)} for i in xrange(len(unique_tabmRNA))]
        #Sort of dictionnary list
        return self.sortdict(frequencydict,"frequency")
    
    def getCommune(self,list1,list2):
        """
        Count Commune elements between lists
        @param list1: First list
        @param list2: Second list
        @return: Count the number of elements
        """
        commune=0
        for element1 in list1:
            if(element1 in list2): commune+=1
        return commune
    
    def compute_similarity(self,computer,matrix_sRNA):
        """
        Compute similarity between target list of each sRNA
        @param computer: Computer object
        @param matrix_sRNA: List of sRNA targets
        """
        simdict=[]
        for i in xrange(len(matrix_sRNA)):
            for j in xrange(i+1,len(matrix_sRNA[i+1:])):
                commune=self.getCommune(matrix_sRNA[i],matrix_sRNA[j])
                #Jaccard distance
                denom=(float(len(matrix_sRNA[i]))+float(len(matrix_sRNA[j]))-float(commune))
                if denom!=0:
                    similarity=float(commune)/denom
                else: similarity=0.0
                simdict+=[{"sRNA_1":computer.allsRNA[i],"sRNA_2":computer.allsRNA[j],"similarity":similarity}]
        return self.sortdict(simdict,"similarity")
    
    def run(self,computer,dbmanage,name,softid):
        """
        Select targets depending on score and compute frequency and similarity between target groups
        @param computer: Computer object
        @param dbmanage: Access to the database
        @param name: Software name
        @param softid: Software id
        """
        computer.matrix,matrix_sRNA,tabmRNA,unique_tabmRNA=self.compute_Davidmatrix(computer,dbmanage,name,softid)
        computer.frequency=self.compute_frequency(computer,tabmRNA,unique_tabmRNA)
        computer.similarity=self.compute_similarity(computer,matrix_sRNA)

