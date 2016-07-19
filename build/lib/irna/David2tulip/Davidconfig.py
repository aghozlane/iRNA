"""
 @brief: Handle DAVID configuration
 @author: Amine Ghozlane
 @version: 1.0 
"""
import ConfigParser,sys,os

class Davidconfig:
    
    def __init__(self,config):
        """
        Instanciate Davidconfig object
        @param config: config path
        """
        self.config = ConfigParser.RawConfigParser()
        if(config!=None): self.davidconfig_file=config
        else: self.davidconfig_file='.%sdavidconfig.cfg'%os.sep
    
    def readconfig(self,davidobject):
        """
        Read david config
        @param davidobject: David pkl
        """
        #If config parser empty
        if(not os.path.isfile(self.davidconfig_file)): self.writeconfig()
        #Read config file
        self.config.read(self.davidconfig_file)
        #Get parameter value
        davidobject.url = self.config.get('DAVID_config', 'url')
        davidobject.idType = self.config.get('DAVID_config', 'idType')
        davidobject.login = self.config.get('DAVID_config', 'login')
        davidobject.thd = self.config.getfloat('DAVID_config', 'thd')
        davidobject.count = self.config.getint('DAVID_config', 'count')
        return davidobject
        
    def writeconfig(self):
        """
        Write david config
        """
        self.config.add_section('DAVID_config')
        self.config.set('DAVID_config', 'url', 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl')
        login = raw_input("Enter your DAVID login :")
        self.config.set('DAVID_config','login',login)
        idtype = raw_input("Enter the idtype :")
        self.config.set('DAVID_config', 'idType',idtype)
        #"UNIPROT_ACCESSION"
        thres = raw_input("Enter the threshold :")
        self.config.set('DAVID_config','thd',thres)
        count = raw_input("Enter the minimum number of genes :")
        self.config.set('DAVID_config','count',count)
        #Write data
        try:
            # Writing our configuration file to 'example.cfg'
            with open(self.davidconfig_file, 'wb') as configfile:
                self.config.write(configfile)
        except IOError: sys.exit("Error : can not open file %s"%self.davidconfig_file)
        except: sys.exit("Something went wrong with %s"%self.davidconfig_file)   