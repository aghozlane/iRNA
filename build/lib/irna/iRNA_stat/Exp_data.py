"""
 @brief: Handle experimental data
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv
from Parser import *

class Exp_data(Parser):
    
    
    def __init__(self,exp_file):
        """
        Instanciate Exp_data object
        """
        Parser.__init__(self)
        self.getExpinf(exp_file)
        
    def getExpinf(self,exp_file):
        """
        Parse experiment file
        @param exp_file: Experiment file
        """
        expdict=[]
        try:
            expReader = csv.reader(open(exp_file, 'rb'), delimiter='\t')
            #pass header
            expReader.next()
            for line in expReader:
                expdict+=[{"combined":line[0].strip().lower()+line[1].strip().lower(),"sRNA":line[0].strip().lower(),"mRNA":line[1].strip().lower(),"sRNA_deb":int(line[2]),"sRNA_fin":int(line[3]),"mRNA_deb":int(line[4]),"mRNA_fin":int(line[5])}]
            #Sort of dictionnary list
            expdict=self.sortdict(expdict,"combined")
            #Get data in separate item
            self.combined= [exp['combined'] for exp in expdict]
            self.sRNA= [exp['sRNA'] for exp in expdict]
            self.mRNA= [exp['mRNA'] for exp in expdict]
            self.sRNA_deb= [exp['sRNA_deb'] for exp in expdict]
            self.sRNA_fin= [exp['sRNA_fin'] for exp in expdict]
            self.mRNA_deb= [exp['mRNA_deb'] for exp in expdict]
            self.mRNA_fin= [exp['mRNA_fin'] for exp in expdict]
            self.lenexp_data=len(self.sRNA)
        except IOError:
            sys.exit("Error : can not open file %s"%exp_file)
        except:
            sys.exit("Something went wrong with %s"%exp_file)
        
    def getExpPair(self,sRNA,mRNA):
        """
        Get the experimental experience for one sRNA-mRNA interaction
        @param sRNA: name of sRNA
        @param mRNA: name of mRNA
        """
        name=sRNA.lower()+mRNA.lower()
        return self.getData(self.combined,name)
