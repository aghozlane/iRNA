"""
 @brief: Handle sRNA information
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv

class sRNA:
    def __init__(self):
        """
        Instanciate name parser object
        """
        pass
    
    def readsRNA(self,sRNA_file):
        """
        Add sRNA - sRNA edges based on their similarity
        @param node_objects: list of node objects
        @param edge_objects: list of edge objects
        """
        self.sRNA_obj=[]
        try:
            srnaReader = csv.reader(open(sRNA_file, 'rb'), delimiter='\t')
            #passage de l'entete
            self.srnareg=srnaReader.next()[0]
            for line in srnaReader:
                if(line[3]=="-"): self.sRNA_obj+=[{"GeneName":line[0].strip().lower(),"begin":int(line[2]),"end":int(line[1]),"complement":True}]
                else: self.sRNA_obj+=[{"GeneName":line[0].strip().lower(),"begin":int(line[1]),"end":int(line[2]),"complement":False}]
            #Sort srna tab
            self.sRNA_obj.sort()
        except IOError:
            sys.exit("Error : can not open file %s"%sRNA_file)
        #except:
        #    sys.exit("Something went wrong with %s"%sRNA_file)