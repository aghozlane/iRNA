"""
 @brief: Handle random interaction information 
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv,re
from Parser import *

class Rand_data(Parser):
    
    def __init__(self,rand_file):
        """
        Instanciate rand_data object
        """
        Parser.__init__(self)
        self.getRandinf(rand_file)
        
    def getRandinf(self,rand_file):
        """
        Load soft parameters
        @param rand_file: Random file
        """
        randict=[]
        try:
            randReader = csv.reader(open(rand_file, 'rb'), delimiter='\t')
            #pass header
            randReader.next()
            for line in randReader:
                if(len(line)==4): randict+=[{"combined":line[0].strip().lower()+line[1].strip().lower(),"soft":line[0].strip().lower(),"sRNA":line[1].strip().lower(),"slope":float(line[2]),"intercept":float(line[3])}]
            #Sort of dictionnary list
            randict=self.sortdict(randict,"combined")
            #Get data in separate item
            self.combined= [rand['combined'] for rand in randict]
            self.soft= [rand['soft'] for rand in randict]
            self.sRNA= [rand['sRNA'] for rand in randict]
            self.intercept= [rand['intercept'] for rand in randict]
            self.slope= [rand['slope'] for rand in randict]
            self.unique_soft=self.getUnique(self.soft)
        except IOError:
            sys.exit("Error : can not open file %s"%rand_file)
        #except:
        #    sys.exit("Something went wrong with %s"%rand_file)
    
    def getallSofts(self,soft):
        """
        Get soft corresponding to one name
        @param soft: Soft name
        """
        softid=[]
        regex=re.compile("^%s_([0-9]+)"%soft.lower())
        for i in self.unique_soft:
            a=regex.match(i)
            if(a): softid+=[int(a.group(1))]
        return softid
    
    def getRandRna(self,soft,sRNA):
        """
        Get analysis of one sRNA by one soft
        @param soft: software name
        @param sRNA: sRNA name
        """
        name=soft.lower()+sRNA.lower()
        return self.getData(self.combined,name)