"""
 @brief: Handle RNA information
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv

class RNA:
    def __init__(self):
        """
        Instanciate name parser object
        """
        pass
    
    def readRNA(self,RNA_file):
        """
        Read RNA file
        @param RNA_file: 
        """
        self.RNA_obj=[]
        try:
            rnaReader = csv.reader(open(RNA_file, 'rb'), delimiter='\t')
            #passage de l'entete
            self.srnareg=rnaReader.next()[0]
            for line in rnaReader:
                if(line[3]=="-"): self.RNA_obj+=[{'GeneID':None,"GeneName":line[0].strip().lower(),"begin":int(line[1]),"end":int(line[2]),"complement":True}]
                else: self.RNA_obj+=[{'GeneID':None,"GeneName":line[0].strip().lower(),"begin":int(line[1]),"end":int(line[2]),"complement":False}]
            #Sort srna tab
            self.RNA_obj.sort()
        except IOError:
            sys.exit("Error : can not open file %s"%RNA_file)
        #except:
        #    sys.exit("Something went wrong with %s"%sRNA_file)