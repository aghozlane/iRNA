"""
 @brief: Compute different statistic on sequences : GC percent and mean length
 @author: Amine Ghozlane
 @version: 1.0
"""
import sys, re, os, string

class Stat:
           
    def __init__(self):
        """
        Stat constructor
        """
        pass

    def ComputeGC(self, name, seq):    
    #def ComputeGC(self, name, seq, dbmanage): 
        """
        Compute the length and the number of GC
        @param name: Name of the sequence
        @param seq: RNA sequence
        """
        #Sequence length
        seqlen=len(seq)
        #Save seq info
        self.SeqInf(name, seqlen)
        #dbmanage.setRNA(name,seqlen)
        self.LenRan+=[seqlen]
        #Add sequence length        
        self.LenSeq+=seqlen
        #Number of GC
        self.GC+=len(seq.translate(None,"ATUatu"))
        #Sequence length
        self.NbSeq+=1

        
    def GetGC(self):
        """
        Get the GC percentage
        @return: The GC percentage
        """
        return(self.GC/self.LenSeq*100.0)
    
    def GetLenSeq(self):
        """
        Get the mean length of sequence
        @return: Mean length
        """
        return(self.LenSeq/self.NbSeq)
    
    def SeqInf(self, name, seqlen):
        """
        Write Seq Inf data
        """
        try: self.file.write("%s\t%d\n"%(name,seqlen))
        except IOError: sys.exit("The program cannot write in %s"%self.out)
    
    def openSeqInf(self):
        """
        Open seq_inf file and write header
        """ 
        try:
            #Open seq inf file
            self.file=open(self.out,"w")
            self.file.write("SeqName\tLength\n")
        except IOError: sys.exit("The program cannot open %s"%self.out)
        except: sys.exit("Something went wrong with %s"%self.out)
        
    def CloseSeqInf(self):
        """
        Close the Seq Inf file
        """
        try: self.file.close()
        except: sys.exit("Something went wrong with %s"%self.out)