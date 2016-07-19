"""
 @brief: Configure and connect to DAVID database
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys,logging
import suds.metrics as metrics
from suds.client import Client
from suds import *
from Parser import *
from Gene_list import *

class David(Parser):

    def __init__(self, david_file=None):
        """
        Instanciate david object
        @param david_file: Path to david file
        """
        Parser.__init__(self)
        self.listType=0
        self.david_chart=[]
        #Load Category Names object
        if(david_file): self.david_chart=self.loadObject(david_file)

    def davidConnection(self):
        """
        Connection to david
        """
        logging.getLogger('suds.client').setLevel(logging.DEBUG)
        client = Client(self.url,faults=False)
        #authentification
        if(client.service.authenticate(self.login)): return client
        else: return False
        
    def analysis(self,graph):
        """
        Proceed to DAVID enrichment and store useful data
        @param graph: Graph object
        """
        counter=1
        client=self.davidConnection()
        if(client!=None):
            listsRNA=self.getsRNA(graph.node_objects)
            lensRNA=len(listsRNA)
            for sRNA in listsRNA:
                print("sRNA submit = %d/%d"%(counter,lensRNA))
                #Get list of Target as a string
                listTarget=','.join(self.getTarget(sRNA, graph.edge_objects))
                #Submit the list to DAVID
                client.service.addList(listTarget, self.idType, sRNA, self.listType)
                #Parse david chart
                self.parseChart(client.service.getChartReport(self.thd, self.count),sRNA)
                #client.service.getDefaultCategoryNames()
                counter+=1
        else: sys.exit("Connection failed to DAVID")
    
    def parseChart(self,david_chart,sRNA):
        """
        Parse DAVID chart data
        @param david_chart: DAVID result
        @param sRNA: sRNA 
        """
        self.david_chart+=[Gene_list(i,sRNA) for i in david_chart[1]]
    
    def writefile(self,results):
        """
        Pickle DAVID information
        @param results: Path to result file
        """
        self.printObject(self.david_chart,results+"david_chart.pkl")    
            
    def setdata(self,node_objects,edge_objects):
        """
        Set data from DAVID into the graph
        @param node_objects: List of node objects
        @param edge_objects: List of edge objects
        @return: Node and edge objects
        """
        for chart in self.david_chart:
            chart.geneIds.encode('utf-8')
            for mRNA in chart.geneIds.split(", "):
                edge=self.getlinkedge(edge_objects,chart.sRNA,mRNA)
                if(edge!=None):
                    #Get the chart category
                    edge.categoryName+=[chart.categoryName]
                    #Get functionnal category
                    edge.termName+=[chart.termName]
                    
        return(node_objects,edge_objects)
