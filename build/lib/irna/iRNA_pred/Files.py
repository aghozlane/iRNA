"""
 @brief: Handle arguments
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, getopt, os

class Files:
    """
    @brief: 
    """
    
    def __init__(self,myrank):
        """
        Instanciate Files object
        """
        self.myrank=myrank
        self.operation=[False]*5
        self.fasta = [None]*4
        self.predict = [None]*6
        self.random = [None]*4
        self.predict[3]=10
        self.fastmode=False
        #Load files
        self.getfiles()

    def usage(self,info):
        """
         Give information to use iRNA
         @param info: Error texte
         @return: Use of iRNA
        """
        text=None
        if(self.myrank==0):
            text = "Compute prediction of sRNA-mRNA interaction.\n\n"
            if(info): text += info
            temp = "Option\t\tCase\t\tFilename\t\tType\tDescription\n"
            text += temp
            text += '-'*(len(temp) + 60)
            text += '\n'
            text += "-a\t\tno\t\tno\t\t\tno\tCreate Random Sequences\n"
            text += "-m\t\tno\t\tno\t\t\tno\tExtract multifasta files to fasta\n"
            text += "-c\t\t-m\t\tno\t\t\tno\tUse to generate couple files\n"
            text += "-r\t\tno\t\tno\t\t\tno\tRun comparisons\n"
            text += "-f\t\t-r(Optional)\tno\t\t\tno\tUse sqlitebck package to get faster sqlite implementation\n"
            text += "--use_extract\t-a -m\t\tno\t\t\tno\tUse GC and mean sequence length for generated sequence\n"
            
            text += "--GC\t\t-a\t\tNumber of GC%\t\tInput\tGC percentage of the generated sequence\n"
            text += "--seqlen\t-a\t\tNumber of letter\tInput\tLength of the generated sequence\n"
            text += "--nbseq\t\t-a\t\tNumber of sequence\tInput\tNumber of sequence to generate\n"
            text += "--outrand\t-a\t\tRepertory/\t\tOutput\tOutput Repertory\n"
            
            text += "--mf_sRNA\t-m\t\tfile.fasta\t\tInput\tsRNA Multifasta\n"
            text += "--mf_mRNA\t-m -r\t\tfile.fasta\t\tInput\tmRNA Multifasta\n"
            text += "--fasta_sRNA\t-m -r\t\tRepertory/\t\tIn/Out\tsRNA fasta repertory\n"
            text += "--fasta_mRNA\t-m -r\t\tRepertory/\t\tIn/Out\tmRNA fasta repertory\n"
            text += "--soft_list\t-r\t\tsoft_list.txt\t\tInput\tmList of soft to test\n"
            text += "--soft_path\t-r\t\tRepertory/\t\tInput\tmPath to the soft root\n"
            text += "--matrix\t-r\t\tmatrix.txt\t\tInput\tMatrix of score\n"
            text += "--result_out\t-r\t\tRepertory/\t\tOutput\tResult repertory\n"
            text += "--comp_list\t-c -r(Optional) comp_list.txt\t\tInput\tComparison list\n"
            
            text += "--buffer_size\t-r (Optional)\tNumber of element\tInput\tBuffer size (default value 10)\n"
        return text

    #Test if necessary document are available
    def case(self):
        """
          Test if necessary document are available
          @param operation: list of options called
          @param fasta: list fasta related information
          @param predict: list comparison related information
        """
        #All cases
        if(self.operation[0]==False and self.operation[1]==False and self.operation[2]==False):
            sys.exit(self.usage("You must indicate an activity : Multifasta extract (option -m),  Generate random sequence (option -a), compute interaction (option -r)\n"))
    
        #Test des fichiers et repertoires
        if(self.operation[0]):
            if (self.fasta[0] == None):
                sys.exit(self.usage("multifasta file (--mf_sRNA) : \"%s\" must be indicated\n" % (self.fasta[0])))
        #Les deux
        if(self.operation[0] or self.operation[1]):
            if (self.fasta[1] == None):
                sys.exit(self.usage("multifasta file (--mf_mRNA) : \"%s\" must be indicated\n" % (self.fasta[1])))
            if (self.fasta[2] == None):
                sys.exit(self.usage("Fasta repertory (--fasta_sRNA) : \"%s\" must be indicated\n" % (self.fasta[2])))
            if (self.fasta[3] == None):
                sys.exit(self.usage("Fasta repertory (--fasta_mRNA) : \"%s\" must be indicated\n" % (self.fasta[3])))    
        #Comparaison
        if(self.operation[1]):
            if (self.predict[0] == None):
                sys.exit(self.usage("Soft list file (--soft_list) : \"%s\" must be indicated\n" % (self.predict[0])))
            if (self.predict[1] == None):
                sys.exit(self.usage("Result repertory (--result_out) : \"%s\" must be indicated\n" % (self.predict[1])))
            if (self.predict[4] == None):
                sys.exit(self.usage("Soft path (--soft_path) : \"%s\" must be indicated\n" % (self.predict[4])))
            if (self.predict[5] == None):
                sys.exit(self.usage("Matrix of score (--matrix) : \"%s\" must be indicated\n" % (self.predict[5])))
        #Random sequence
        if(self.operation[2]):
            if(self.random[2]==None):
                sys.exit(self.usage("Number of sequence to generate (--nbseq) : \"%s\" must be indicated\n"%self.random[2]))
            if(self.random[3]==None):
                sys.exit(self.usage("Random sequence output repertory (--outrand) : \"%s\" must be indicated\n"%self.random[3]))
            #case use_extract true and m not indicated
            if(self.operation[3]):
                if(not self.operation[0]):
                    sys.exit(self.usage("Option (-m) must be indicated with (--use_extract)\n"))
            #case m and use_extract not indicated
            if(not self.operation[1] and not self.operation[3]):
                if(not self.random[0]):
                    sys.exit(self.usage("GC percentage (--GC) : \"%s\" must be indicated\n"%self.random[0]))
                if(not self.random[1]):
                    sys.exit(self.usage("Sequence length (--seqlen) : \"%s\" must be indicated\n"%self.random[1]))
        #Couple
        if(self.operation[4]):
            if(not self.operation[0]):
                sys.exit(self.usage("Option (-m) must be indicated with (-c)\n"))
    
    def data_format(self):
        """
        Check if information are correct
        """
        #Test des fichiers et repertoires
        if(self.operation[0]):
            #test si c'est bien un fichier
            if(not os.path.isfile(self.fasta[0])):
                sys.exit(self.usage("Error with \"%s\" : --mf_sRNA requires a Multifasta file\n"%self.fasta[0]))
        #Les deux
        if(self.operation[0] or self.operation[1]):
            if(not os.path.isfile(self.fasta[1])):
                sys.exit(self.usage("Error with \"%s\" : --mf_mRNA requires a Multifasta file\n"%self.fasta[1]))
            if(not os.path.isdir(self.fasta[2])):
                sys.exit(self.usage("Error with \"%s\" : --fasta_sRNA requires a path to a repertory\n"%self.fasta[2]))
            if(not os.path.isdir(self.fasta[3])):
                sys.exit(self.usage("Error with \"%s\" : --fasta_mRNA requires a path to a repertory\n"%self.fasta[3]))    
        #Comparaison
        if(self.operation[1]):
            if(not os.path.isfile(self.predict[0])):
                sys.exit(self.usage("Error with \"%s\", --soft_list requires a file\n"%self.predict[0]))
            if(not os.path.isdir(self.predict[1])):
                sys.exit(self.usage("Error with \"%s\" : --result_out requires a path to a repertory\n"%self.predict[1]))
            if(not os.path.isdir(self.predict[4])):
                sys.exit(self.usage("Error with \"%s\" : --soft_path requires a path to a repertory\n" % (self.predict[4])))
            if(not os.path.isfile(self.predict[5])):
                sys.exit(self.usage("Error with \"%s\", --matrix requires a file\n"%self.predict[5]))        
        #Random sequence
        if(self.operation[2]):
            if(not os.path.isdir(self.random[3])):
                sys.exit(self.usage("Error with \"%s\", --outrand requires a path to a repertory\n"%self.random[3]))
             
    
    #Determine les fichiers fournis en arguments
    def getfiles(self):
        """
        Determine the files provided as arguments
        @return: Choosen options
        """
        #Sans argument
        if len(sys.argv) <= 1:
            sys.exit("Do 'iRNA_pred.py -h' for a usage summary")
        #test des option
        try:
            (opts, args) = getopt.getopt(sys.argv[1:], "afcmrh", ["GC=","seqlen=","nbseq=","outrand=","use_extract","mf_sRNA=", "mf_mRNA=", "fasta_sRNA=", "fasta_mRNA=", "comp_list=", "result_out=", "soft_list=", "buffer_size=","soft_path=", "matrix=","help"])
        except getopt.GetoptError, err:
            # print help information and exit:
            print str(err) # will print something like "option -a not recognized"
            sys.exit(self.usage(None))
            
        temp = "iRNA_pred chosen activity : "
        #Identification of options
        for (o, a) in opts:
            if o == "-m":
                self.operation[0] = True
                temp += 'Multifasta extract\t'
            elif o == "-r":
                self.operation[1] = True
                temp += 'Compute interaction\t'
            elif o == "-a":
                self.operation[2] = True
                temp += 'Generate random sequence\t'
            elif o == "-c":
                self.operation[4] = True
                temp += 'Generate couple files\t'
            elif o == "--use_extract":
                self.operation[3] = True
                temp += 'Use sequences data\t'
            elif o == "-f":
                self.fastmode = True
                temp += 'Use sequences data\t'
            elif o == "--GC":
                try: self.random[0] = float(a)
                except: sys.exit(self.usage("Option \"GC\" must be followed by a number: %s\n"%(a)))
            elif o == "--seqlen":
                try: self.random[1] = int(a)
                except: sys.exit(self.usage("Option \"seqlen\" must be followed by a number: %s\n"%(a)))
            elif o == "--nbseq":
                try: self.random[2] = int(a)
                except: sys.exit(self.usage("Option \"nbseq\" must be followed by a number: %s\n"%(a)))
            elif o == "--outrand": self.random[3] = a
            elif o == "--mf_sRNA": self.fasta[0] = a
            elif o == "--mf_mRNA": self.fasta[1] = a
            elif o == "--fasta_sRNA": self.fasta[2] = a
            elif o == "--fasta_mRNA": self.fasta[3] = a
            elif o == "--soft_list": self.predict[0] = a
            elif o == "--result_out": self.predict[1] = a
            elif o == "--comp_list":
                self.predict[2] = a
                if(not os.path.isfile(self.predict[2])):
                    sys.exit(self.usage("Error with %s, --soft_list required a file\n"%self.predict[2]))
            elif o == "--buffer_size":
                try:
                    self.predict[3] = int(a)
                except:
                    sys.exit(self.usage("Option \"buffer_size\" must be followed by a number: %s\n"%(a)))
            elif o == "--soft_path": self.predict[4] = a
            elif o == "--matrix": self.predict[5] = a
            elif o in ("-h", "--help"): sys.exit(self.usage(None))
            else: assert False, "unhandled option"
        #Print selected options
        if( (self.operation[0] or self.operation[1])and self.myrank==0): print temp            
        #Verification of cases    
        self.case()
        self.data_format()