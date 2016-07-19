# -*- coding: utf-8 -*-
"""
 @brief: Get files
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, getopt, os

class Files:
        
    def __init__(self):
        """
        Instanciate Files object
        """
        self.multilist = None
        self.filter = None
        self.similarity = None
        self.pValue = None 
        self.results = None
        self.david = False
        self.interact = None
        self.enrichment = None
        self.fastmode= False       
        self.iRNA_db = None
        self.config = None
        self.name = None
        self.getfiles()

    def usage(self,info):
        """
         Give information to use iRNA
         @param info: Error texte
         @return: Use of iRNA
        """
        text=None
        text = "Create nodes and edge csv for Tulip from multiple datasets.\n\n"
        if(info): text += info
        temp = "Option\t\t\tfile\t\tDescription\n"
        text += temp
        text += '-'*(len(temp) + 60)
        text += '\n'
        text += "-m, --multilist\t\tmultilist.txt\tMultilist file from iRNA_stat\n"
        text += "-r, --results\t\tno\t\tPath to result repertory\n"
        text += "-s, --similarity\tsimilarity.txt\tSet similarity of sRNA groups\n"
        text += "-p, --pValue\t\tpValue.txt\tpValue of selected files\n"
        text += "-d, --david\t\tno\t\tSubmit multilist to DAVID\n"
        text += "-c, --config\t\tDavid.cfg\tConfiguration file for David\n"
        text += "-i, --interact\t\tinteraction.txt\tPath to interaction file\n"
        text += "-e, --enrichment\tDavid.pkl\tPath to DAVID enrichment file\n"
        text += "-y, --iRNA_db\t\tiRNA.db\tPath to iRNA db file\n"
        text += "-f, --fastmode\t\tno\t\tUse sqlitebck package to get faster sqlite implementation (optional)\n"
        text += "-l, --filter\t\tfilter.txt\tFilter multilist for only some metabolites and use it for analysis\n"
        text += "-n, --name\t\tname.txt\tRename uniprot genes using datafile\n"
        return text
 
    def case(self):
        """
          Test if necessary document are available
          @param operation: list of options called
          @param fasta: list fasta related information
          @param predict: list comparison related information
        """
        #Test des fichiers et repertoires
        if(not self.multilist):
            sys.exit(self.usage("multilist (-m,--multilist) : \"%s\" must be indicated\n" % (self.multilist)))
        if(not self.results):
            sys.exit(self.usage("results (-r,--results) : \"%s\" must be indicated\n" % (self.results)))
                               
    def data_format(self):
        """
        Check if information are correct
        """
        #Run without arguments
        if len(sys.argv)== 1:
            sys.exit(self.usage(None)) 
        #Test multilist file argument
        if self.multilist:
            if(not os.path.isfile(self.multilist)): sys.exit(self.usage("Error with \"%s\" : -m required a multilist file\n"%self.multilist))            
        #Test result repertory argument
        if self.results:
            if(not os.path.isdir(self.results)): sys.exit(self.usage("Error with \"%s\" : -r required a repertory\n"%self.results)) 
        #Test pValue file argument
        if self.pValue:
            if(not os.path.isfile(self.pValue)): sys.exit(self.usage("Error with \"%s\" : -p required a file\n"%self.pValue)) 
        #Test enrichment repertory argument
        if self.enrichment:
            if(not os.path.isfile(self.enrichment)): sys.exit(self.usage("Error with \"%s\" : -e required a file\n"%self.enrichment)) 
        #Test filter file
        if self.filter:
            if(not os.path.isfile(self.filter)): sys.exit(self.usage("Error with \"%s\" : -f required a file\n"%self.filter)) 
        #Test name file
        if self.name:
            if(not os.path.isfile(self.name)): sys.exit(self.usage("Error with \"%s\" : -n required a file\n"%self.name)) 
        #Test similarity file
        if self.similarity:
            if(not os.path.isfile(self.similarity)): sys.exit(self.usage("Error with \"%s\" : -s required a file\n"%self.similarity)) 
        #Test interact file
        if self.interact:
            if(not os.path.isfile(self.interact)): sys.exit(self.usage("Error with \"%s\" : -i required a file\n"%self.interact)) 
        #Test david arguments
        if self.enrichment and self.david: print("Warning : enrichment file \"%s\" will be ignored\n"%self.enrichment)
        #Test iRNA db
        if self.iRNA_db:
            if(not os.path.isfile(self.iRNA_db)): sys.exit(self.usage("Error with \"%s\" : -y required a file\n"%self.iRNA_db))
        #Test config
        if self.config:
            if(not os.path.isfile(self.config)): sys.exit(self.usage("Error with \"%s\" : -c required a file\n"%self.config))
                
    #Determine les fichiers fournis en arguments
    def getfiles(self):
        """
        Determine the files provided as arguments
        @return: Choosen options
        """
        #Sans argument
        if len(sys.argv) <= 1: sys.exit("Do 'David2tulip.py -h' for a usage summary")        
        #test des option
        try: (opts, args) = getopt.getopt(sys.argv[1:], "m:r:l:s:p:i:e:dy:c:n:fh", ["multilist=","results=","filter=","similarity=","pValue=","interact=","enrichment=","iRNA_db=","name=","config","david","fastmode","help"])
        except getopt.GetoptError, err:
            # print help information and exit:
            print str(err) # will print something like "option -a not recognized"
            sys.exit(self.usage(None))
        #Identification of options
        for (o, a) in opts:
            if o in ("-m","--multilist"): self.multilist = a
            elif o in ("-l", "--filter"): self.filter = a
            elif o in ("-n", "--name"): self.name = a
            elif o in ("-s","--similarity"): self.similarity = a
            elif o in ("-p","--pValue"): self.pValue = a 
            elif o in ("-r", "--results"): self.results = a
            elif o in ("-d","--david"): self.david = True
            elif o in ("-i","--interact"): self.interact = a
            elif o in ("-e","--enrichment"): self.enrichment = a
            elif o in ("-f","--fastmode"): self.fastmode= True
            elif o in ("-y","--iRNA_db"): self.iRNA_db = a
            elif o in ("-c","--config"): self.config = a
            elif o in ("-h", "--help"): sys.exit(self.usage(None))
            else: assert False, "unhandled option"            
        #Verification of cases    
        self.case()
        self.data_format()