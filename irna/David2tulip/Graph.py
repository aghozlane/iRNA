# -*- coding: utf-8 -*-
"""
 @brief: Handle Graph
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, bisect,re
from Node import *
from Edge import *

class Graph:

    def __init__(self):
        """
        Instanciate graph object
        """
        self.node_objects=[]
        self.edge_objects=[]
        self.fcategoryid=["INTERPRO","KEGG_PATHWAY","PIR_SUPERFAMILY","SMART"]
        self.gocategoryid=["GOTERM_BP_FAT","GOTERM_CC_FAT","GOTERM_MF_FAT"]
        self.fregex=re.compile("(.+):(.+)")
        self.goregex=re.compile("(.+)~(.+)")
    
    def parse(self,obj):
        """
        Parse and information to node and edge objects
        @param obj: Parser-linked object
        """
        self.node_objects,self.edge_objects=obj.setdata(self.node_objects,self.edge_objects)
    
    def printnodes(self,results):
        """
        Write nodes csv
        @param results: Path to result repertory
        """
        results=results+"nodes_iRNA.csv"
        try:
            #Open file
            nodes_file=open(results,"wt")
            #Write heading
            nodes_file.write("Node\tGender\tname\tcode\tlength\tgroup\n")
            for i in self.node_objects:
                nodes_file.write(i.println())
            nodes_file.close()
        except IOError: sys.exit("Error : can not open file %s"%results)
        except: sys.exit("Something went wrong with %s"%results) 
            
    
    def getCategories(self):
        """
        Get active categories
        @return: Unique categories
        """
        listCategoryName=[]
        for edge in self.edge_objects:
            listCategoryName+=edge.categoryName
        if(len(listCategoryName)>0): return({}.fromkeys(listCategoryName).keys())
        else: return None
    
    def getNumCategory(self,element,liste):
        """
        Get the position of on element in a list
        @param element: an element of the list
        @param liste: a list 
        """
        return bisect.bisect_left(liste,element)
    
    def verif(self):
        """
        Verif value for sRNA_100
        """
        listTarget=[]
        for edge in self.edge_objects:
            if(edge.node1.name=="sRNA_100" and edge.category==1):
                listTarget+=[edge.node2.name]
        print(listTarget)
    
    def filterTermName(self,category,term):
        """
        Detect the method to use for filtering
        """
        test=None
        #Detect category and parse
        if(category in self.fcategoryid): test=self.fregex.match(term)
        elif(category in self.gocategoryid): test=self.goregex.match(term)
        else:
            filterm=term
            filterid=None
        #Set filtering
        if(test):
            filterid=test.group(1)
            filterm=test.group(2)
        return(filterm,filterid)
    
    def printedges(self,results):
        """
        Write edges csv
        @param results:  Path to result repertory
        """
        results=results+"edges_iRNA.csv"
        #self.verif()        
        header="source_id\tdestination_id\tcategory\tsimilarity\tdistance\tpValue\tsRNA_positions\tmRNA_positions\tDatabase_recurrency"
        #Get categories
        categories=self.getCategories()
        if(categories):
            db=""
            #Add id categories
            for i in categories:
                header+="\tDBID_%s"%i
                db+="\tDB_%s"%i
            #Db categories
            header+=db
        try:
            #Open file
            edges_file=open(results,"wt")
            #Write heading
            edges_file.write(header+"\n")
            #Print David result
            if(categories):
                for edge in self.edge_objects:
                    text=edge.println().translate(None,"\n")
                    dbdata=""
                    #print categories
                    for i in categories:
                        lencategory=len(edge.categoryName)
                        termName=[]
                        idtab=[]
                        # Case category available
                        if(lencategory>0):
                            for pos in range(len(edge.categoryName)): 
                                if(i==edge.categoryName[pos]):
                                    filtermName,filteridtab=self.filterTermName(edge.categoryName[pos],edge.termName[pos])
                                    termName+=[filtermName]
                                    if(filteridtab): idtab+=[filteridtab]
                            if(len(idtab)>0): text+="\t(\""+"\",\"".join(idtab)+"\")"
                            else: text+="\t()"
                            if(len(termName)>0): dbdata+="\t(\""+"\",\"".join(termName)+"\")"
                            else: dbdata+="\t()"
                        # Case no category
                        else:
                            text+="\t()"
                            dbdata+="\t()"
                    #Write data
                    text+=dbdata
                    edges_file.write(text+"\n")                   
            else:
                #Write linked sRNA-mRNA
                for i in self.edge_objects:
                    edges_file.write(i.println())            
            #Close file
            edges_file.close()
        except IOError:
            sys.exit("Error : can not open file %s"%results)
        #except:
        #    sys.exit("Something went wrong with %s"%results) 
    
    def writefile(self, results):
        """
        Write node and edge csv
        @param results: Path to result repertory
        """
        self.printnodes(results)
        self.printedges(results)