"""
 @brief: Handle execution data
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, csv, bisect
from Parser import *


class Execute(Parser):
    
    
    def __init__(self,exec_file):
        """
        Instanciate Execute object
        """
        Parser.__init__(self)
        self.getExecinf(exec_file)
    
    def getExecinf(self,exec_file):
        """
        Parse exection information file
        @param exec_file: Execution file
        """
        execdict=[]
        try:
            exec_infReader = csv.reader(open(exec_file, 'rb'), delimiter='\t')
            #pass header
            exec_infReader.next()
            for line in exec_infReader:
                if(line[0].strip()!="iRNA"):
                    execdict+=[{"soft":line[0].strip(),"duration":float(line[1])}]
            #Sort of dictionnary list
            execdict=self.sortdict(execdict,"soft")
            #Get data in separate item
            self.soft= [execution['soft'] for execution in execdict]
            self.duration= [execution['duration'] for execution in execdict]
        except IOError:
            sys.exit("Error : can not open file %s"%exec_file)