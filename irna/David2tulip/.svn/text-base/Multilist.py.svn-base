# -*- coding: utf-8 -*-
"""
 @brief: Handle multilist information
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, re, csv
from Edge import *
from Parser import *

class Multilist(Parser):
    
    def __init__(self,multilist_file):
        """
        Instanciate multilist object
        @param multilist_file: Multilist file
        """
        Parser.__init__(self)
        self.multilist_file=multilist_file

#    def treatSequence(self,sequence):
#        """
#        Filter multilist data
#        @param sequence: line in multilist file
#        @return: cleaned line
#        """
#        test=sequence.translate(None,'"').replace('NA\t',"").replace('\tNA',"").replace('NA\n',"\n").translate(None,'\n').translate(None,'\r')
#        if(test=="NA"):
#            return None
#        return(test.split("\t"))
   
#    def getnodelists(self):
#        """
#        Extract unique nodes from multilist file
#        @return: lists of unique nodes for header and  tail
#        """
#        tail_list=[]
#        regex=re.compile('\"[1-9]+\"\t(.+)\n')
#        try:
#            multilist=open(self.multilist_file,"rt")
#            #Filter header        
#            header_list=self.treatSequence(multilist.readline())
#            for i in multilist:
#                #Regex detection
#                a=regex.match(i)
#                if(a):
#                    #Filter target
#                    treated=self.treatSequence(a.group(1))                    
#                    if(treated):                    
#                        tail_list+=treated
#            #Selection of unique id
#            tail_list={}.fromkeys(tail_list).keys()
#            #Close file
#            multilist.close()
#        except IOError:
#            sys.exit("Error : can not open file %s"%self.multilist_file)
#        except:
#            sys.exit("Something went wrong with %s"%self.multilist_file)   
#        return(header_list,tail_list)
    def getnodelists(self):
        """
        Read nodes
        """
        tail_list=[]
        header_list=[]
        try:
            multilistReader = csv.reader(open(self.multilist_file, 'rb'), delimiter='\t')
            #Write head text
            header_list=multilistReader.next()
            for row in multilistReader:
                tail_list+=row
            tail_list={}.fromkeys(tail_list).keys()
            tail_list.remove('')
        except IOError: sys.exit("Error : can not open file %s"%self.multilist_file)
        return(header_list,tail_list)
    
    def getPosTable(self,header_nodes,header):
        """
        @param header_nod: 
        """
        indextable=[]
        for sRNA in header:
            i=bisect.bisect_left(header_nodes,sRNA)
            if(i!=len(header_nodes) and header_nodes[i]==sRNA): indextable+=[i]
        return indextable

    def createedges(self, header_nodes, tail_nodes):
        """
        Create edges
        @param header_nodes: list of header nodes
        @param tail_nodes: list of tail nodes
        @return: list of edges
        """
        listedges=[]
        try:
            multilistReader = csv.reader(open(self.multilist_file, 'rb'), delimiter='\t',  quotechar='"')
            #passage de l'entete
            header=multilistReader.next()
            #Get the index
            indextable=self.getPosTable(header_nodes,header)
            for line in multilistReader:
                #for j in range(1,len(line)):
                for j in range(0,len(line)):    
                    #If it's a true element
                    #if line[j]!="NA":
                    if line[j]!='':
                        node_element=self.getlinknode(tail_nodes,line[j])
                        #print(line[j]+" correspond a "+node_element.name+" et est rattaché à "+header_nodes[indextable[j-1]].name)                        
                        if(node_element!=None):
                            #Add edge element
                            listedges+=[Edge(header_nodes[indextable[j]],node_element)]
                            #Set group element
                            #node_element.addgroupelement(j-1)
                            node_element.addgroupelement(header_nodes[indextable[j]].num)
            #Sort edges
            listedges.sort()
        except IOError: sys.exit("Error : can not open file %s"%self.multilist_file)
        #except:
        #    sys.exit("Something went wrong with %s"%self.multilist_file) 
        return listedges

    def setdata(self,node_objects,edge_objects):
        """
        Build the graph
        @param node_objects: list of node objects
        @param edge_objects: list of edge objects
        @return: List of node and edge objects
        """
        #Read node information on multilist file
        header_list,tail_list=self.getnodelists()
        #Build node objects
        sRNA_nodes=self.createnodes(header_list,1,0)
        mRNA_nodes=self.createnodes(tail_list,2,len(header_list)+1)
        node_objects=sRNA_nodes+mRNA_nodes
        #sort the new list
        node_objects.sort()
        #Read edge information on multilist file
        edge_objects=self.createedges(sRNA_nodes, mRNA_nodes)
        return(node_objects,edge_objects)
    
    def getFilter(self,filter_file):
        """
        Parse gene known 
        @param filter_file: GO file
        @return: List of mRNA known from the GO
        """
        go_list=None
        try:
            filterReader=csv.reader(open(filter_file, 'rb'), delimiter='\t',  quotechar='"')
            #Read accession number
            go_list=[row[0] for row in filterReader]
            #Selection of unique id
            go_list={}.fromkeys(go_list).keys()  
        except IOError: sys.exit("Error : can not open file %s"%filter_file)
        except: sys.exit("Something went wrong with %s"%filter_file)     
        return go_list
    
    def writefilterlist(self,go_list,results):
        """
        Write the filtered multilist
        @param go_list: List of mRNA known
        @param results: Path to result repertory
        """
        self.filtered=results+'filtered_multilist.txt'
        try:
            multilistReader = csv.reader(open(self.multilist_file, 'rb'), delimiter='\t')
            writer= csv.writer(open(self.filtered,"wb"), delimiter='\t',quoting=csv.QUOTE_ALL)
            #Write head text
            writer.writerow(multilistReader.next())
            for row in multilistReader:
                for j in range(1,len(row)):
                    #If it's a true element
                    if(not row[j] in go_list): row[j]="NA"
                writer.writerow(row)                        
        except IOError: sys.exit("Error : can not open file %s or %s"%(self.multilist_file,self.filtered))
        #except:
        #    sys.exit("Something went wrong with %s or %s"%(self.multilist_file,self.filtered))
    
    def GOfilter(self,filter_file,results):
        """
        Filter mRNA based on GO-known gene
        @param filter_file: GO file
        @param results: Path to result repertory
        """
        #List obtained by the geneontology
        go_list=self.getFilter(filter_file)
        #If Go list obtained
        if(go_list!=None): self.writefilterlist(go_list,results)
        