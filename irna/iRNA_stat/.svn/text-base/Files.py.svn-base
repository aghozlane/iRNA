# -*- coding: utf-8 -*-
"""
 @brief: Get files
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, getopt, os

class Files:
        
    def __init__(self,myrank):
        """
        Instanciate Files object
        """
        self.iRNA_db = None
        self.soft_inf = None
        self.rand_inf = None
        self.exp_inf = None 
        self.random = False
        self.save = None
        self.results = None
        self.exec_inf = None
        self.thres_inf=None
        self.myrank=myrank
        self.overwrite=False
        self.fastmode=False
        self.pValue=None
        self.getfiles()

    def usage(self,info):
        """
         Give information to use iRNA_stat
         @param info: Error texte
         @return: Use of iRNA
        """
        text=None
        if(self.myrank==0):
            text = "Statistical analysis of sRNA-mRNA predicted interaction.\n\n"
            if(info):
                text += info
            temp = "Option\t\tcase\t\tfile\t\tDescription\n"
            text += temp
            text += '-'*(len(temp) + 60)
            text += '\n'
            text += "-d, --iRNA_db\tall\t\tiRNA.db\tiRNA db file\n"
            text += "-r, --results\tall\t\tno\t\tPath to output repertory\n"
            text += "-a, --random\t-d -r\t\tno\t\tRandom analysis\n"
            text += "-i, --soft_inf\t-d -r\t\tsoft_inf.txt\tsoft_inf\n"
            text += "-e, --exp_inf\t-d -r\t\texp_inf.txt\tPredict position analysis\n"
            text += "-n, --rand_inf\t-d -r\t\tsoft_param.txt\tCompute pValue\n"
            text += "-p, --pValue\t-d -r -n\tno\t\tType of pValue calculation: 1: global 2: self \n"
            text += "-t, --thres_inf\t-d -r\t\tpValuethres.txt\tSelect pair with their pValue threshold\n"
            text += "-x, --exec_inf\t-d -r\t\texectime.txt\tPlot execution time\n"
            text += "-s, --save\t-d -r\t\tno\t\tSave temporary object to this path\n"
            text += "-o, --overwrite\t-d -r -s\tno\t\tOverwrite previous save files\n"
            text += "-f, --fastmode\tOptional\tno\t\tUse sqlitebck package to get faster sqlite implementation\n"
        return text

    #Test if necessary document are available
    def case(self):
        """
          Test if necessary document are available
          @param operation: list of options called
          @param fasta: list fasta related information
          @param predict: list comparison related information
        """
        #Test des fichiers et repertoires
        if(not self.iRNA_db):
            sys.exit(self.usage("iRNA_db (-d,--iRNA_db) : \"%s\" must be indicated\n" % (self.iRNA_db)))
        if(not self.soft_inf):
            sys.exit(self.usage("soft_inf (-i,--soft_inf) : \"%s\" must be indicated\n" % (self.soft_inf)))
        if(not self.results):
            sys.exit(self.usage("results (-r,--results) : \"%s\" must be indicated\n" % (self.results)))
        if(self.random and not self.pValue):
            sys.exit(self.usage("pValue type (-p,--pValue) : \"%s\" must be indicated\n" % (self.pValue)))
    
    def data_format(self):
        """
        Check if information are correct
        """
        #Test des fichiers et repertoires
        if(self.iRNA_db):
            if(not os.path.isfile(self.iRNA_db)):
                sys.exit(self.usage("iRNA_db (-d,--iRNA_db) : \"%s\" requires a file\n" % (self.iRNA_db)))
        if(self.results):
            if(not os.path.isdir(self.results)):
                sys.exit(self.usage("results (-r,--results) : \"%s\" requires a path to a repertory\n" % (self.results)))
        if(self.rand_inf):
            if(not os.path.isfile(self.rand_inf)):
                sys.exit(self.usage("rand_inf (-n,--rand_inf) : \"%s\" requires a file\n" % (self.rand_inf)))
        if(self.soft_inf):
            if(not os.path.isfile(self.soft_inf)):
                sys.exit(self.usage("soft_inf (-i,--soft_inf) : \"%s\" requires a file\n" % (self.soft_inf)))
        if(self.exec_inf):
            if(not os.path.isfile(self.exec_inf)):
                sys.exit(self.usage("soft_inf (-x,--exec_inf) : \"%s\" requires a file\n" % (self.exec_inf)))
        if(self.save):
            if(not os.path.isdir(self.save)):
                sys.exit(self.usage("save (-s,--save) : \"%s\" requires a path to a repertory\n" % (self.save)))
        if(self.exp_inf):
            if(not os.path.isfile(self.exp_inf)):
                sys.exit(self.usage("exp_inf (-e,--exp_inf) : \"%s\" requires a file\n" % (self.exp_inf)))
        if(self.thres_inf):
            if(not os.path.isfile(self.thres_inf)):
                sys.exit(self.usage("thres_inf (-t,--thres_inf) : \"%s\" requires a file\n" % (self.thres_inf)))
        if(self.pValue):
            if(self.pValue!=1 and self.pValue!=2):
                sys.exit(self.usage("pValue type (-p,--pValue) : \"%s\" must be equal to 1 or 2\n" % (self.pValue)))    
    
    #Determine les fichiers fournis en arguments
    def getfiles(self):
        """
        Determine the files provided as arguments
        @return: Choosen options
        """        
        #Sans argument
        if len(sys.argv) <= 1:
            sys.exit("Do 'iRNA_stat.py -h' for a usage summary")
        #test des option
        try:
            (opts, args) = getopt.getopt(sys.argv[1:], "d:i:e:n:r:s:t:x:p:afoh", ["random","overwrite","fastmode","pValue=","iRNA_db=","soft_inf=","exp_inf=","rand_inf=", "results=","save=","exec_inf=","thres_inf=","help"])
        except getopt.GetoptError, err:
            # print help information and exit:
            print str(err) # will print something like "option -a not recognized"
            sys.exit(self.usage(None))
        #Identification of options
        for (o, a) in opts:
            if o in ("-d","--iRNA_db"): self.iRNA_db = a
            elif o in ("-i", "--soft_inf"): self.soft_inf = a
            elif o in ("-n","--rand_inf"): self.rand_inf = a
            elif o in ("-e","--exp_inf"): self.exp_inf = a 
            elif o in ("-r", "--results"): self.results = a
            elif o in ("-s","--save"): self.save = a
            elif o in ("-p","--pValue"):
                try: self.pValue=int(a)
                except: sys.exit("Impossible to get pValue : %s, integer value expected"%a)
            elif o in ("-x","--exec_inf"): self.exec_inf = a
            elif o in ("-t","--thres_inf"): self.thres_inf = a
            elif o in ("-a","--random"): self.random = True
            elif o in ("-f","--fastmode"): self.fastmode= True
            elif o in ("-o","--overwrite"): self.overwrite = True
            elif o in ("-h", "--help"): sys.exit(self.usage(None))
            else: assert False, "unhandled option"            
        #Verification of cases    
        self.case()
        self.data_format()