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
        self.genbank = None
        self.results = None
        self.sRNA_file=None
        self.mRNA_file = None
        self.begin = -150
        self.end = 50 
        self.nucleo = None
        self.complete = False
        self.getfiles()

    def usage(self,info):
        """
         Give information to use iRNA_seq
         @param info: Error texte
         @return: Use of iRNA
        """
        text=None
        text = "Create mRNA multifasta file for iRNA.\n\n"
        if(info): text += info
        temp = "Option\t\t\tfile\t\tDescription\n"
        text += temp
        text += '-'*(len(temp) + 60)
        text += '\n'
        text += "-g, --genbank\t\tgenbank.gbk\tGenbank file\n"
        text += "-m, --mRNA_file\t\tmrna.csv\tIndicate extra mRNA gene\n"
        text += "-s, --sRNA_file\t\tsrna.csv\tIndicate which gene correspond to sRNA\n"
        text += "-b, --begin\t\tno\t\tBegin position (default: -150 before start codon)\n"
        text += "-e, --end\t\tno\t\tEnd position (default: +50 after start codon)\n"
        text += "-r, --results\t\tno\t\tPath to result repertory\n"
        text += "-c, --complete\t\tno\t\tComplete genes (optional)\n"
        return text
 
    def case(self):
        """
          Test if necessary document are available
          @param operation: list of options called
          @param fasta: list fasta related information
          @param predict: list comparison related information
        """
        #Test des fichiers et repertoires
        if(not self.genbank):
            sys.exit(self.usage("genbank (-g,--genbank) : \"%s\" must be indicated\n" % (self.genbank)))
        if(not self.results):
            sys.exit(self.usage("results (-r,--results) : \"%s\" must be indicated\n" % (self.results)))
        if(self.begin>self.end):
            sys.exit(self.usage("Begin position of the mRNA (%d) must be inferior to end position of the mRNA (%d)\n" % (self.begin,self.end)))
                               
    def data_format(self):
        """
        Check if information are correct
        """
        #Run without arguments
        if len(sys.argv)== 1:
            sys.exit(self.usage(None)) 
        #Test genbank file argument
        if self.genbank:
            if(not os.path.isfile(self.genbank)):
                sys.exit(self.usage("Error with \"%s\" : -g required a genbank file\n"%self.genbank))
        #Test sRNA file argument
        if self.mRNA_file:
            if(not os.path.isfile(self.mRNA_file)):
                sys.exit(self.usage("Error with \"%s\" : -m required a mRNA_file file\n"%self.mRNA_file))       
        #Test sRNA file argument
        if self.sRNA_file:
            if(not os.path.isfile(self.sRNA_file)):
                sys.exit(self.usage("Error with \"%s\" : -s required a sRNA_file file\n"%self.sRNA_file))              
        #Test result repertory argument
        if self.results:
            if(not os.path.isdir(self.results)):
                sys.exit(self.usage("Error with \"%s\" : -r required a repertory\n"%self.results))
                
    #Determine les fichiers fournis en arguments
    def getfiles(self):
        """
        Determine the files provided as arguments
        @return: Choosen options
        """
        #Sans argument
        if len(sys.argv) <= 1: sys.exit("Do 'iRNA_seq.py -h' for a usage summary")        
        #test des option
        try: (opts, args) = getopt.getopt(sys.argv[1:], "g:r:m:s:b:e:ch", ["genbank=","mRNA_file=","sRNA_file=","results=","begin=","end=","complete","help"])
        except getopt.GetoptError, err:
            # print help information and exit:
            print str(err) # will print something like "option -a not recognized"
            sys.exit(self.usage(None))
        #Identification of options
        for (o, a) in opts:
            if o in ("-g","--genbank"): self.genbank = a
            elif o in ("-b", "--begin"):
                try: self.begin = int(a)
                except: sys.exit(self.usage("Option \"b,begin\" must be followed by a number: %s\n"%(a)))            
            elif o in ("-e","--end"):
                try: self.end = int(a)
                except: sys.exit(self.usage("Option \"e,end\" must be followed by a number: %s\n"%(a)))
            elif o in ("-c","--complete"): self.complete = True
            elif o in ("-s", "--sRNA_file"): self.sRNA_file = a
            elif o in ("-m", "--mRNA_file"): self.mRNA_file = a
            elif o in ("-r", "--results"): self.results = a
            elif o in ("-h", "--help"): sys.exit(self.usage(None))
            else: assert False, "unhandled option"            
        #Verification of cases    
        self.case()
        self.data_format()