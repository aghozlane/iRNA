"""
 @brief: Generate Random Sequences
 @author: Amine Ghozlane
 @version: 1.0
"""
#!/usr/bin/python
# -*- coding: utf-8 -*-
import  sys, numpy as np
#pygsl.rng,

class RanRNA:

    #UniValues=pygsl.rng.uni()
    def __init__(self, GC, SeqLength, NbSeq, OutRep):
        """
        @param GC: GC percentage
        @param SeqLength: Length of sequences
        @param NbSeq: Number of sequences
        @param OutRep: Output repertory
        """
        self.GC=GC/100.0
        self.SeqLength=SeqLength
        self.NbSeq=NbSeq
        self.OutRep=OutRep
    
    def GetValues(self,seqlength):
        """
        Get double random values
        @return: Table of double random Values
        """
        #return self.UniValues(self.SeqLength)
        return np.random.uniform(size=seqlength)
    
    def GetIntValues(self,val,nbseq):
        """
        Get integer random values
        @return: Table of integer random values
        """
        return np.random.randint(val,size=nbseq)
    
    def SetLetters(self, L1, L2):
        """
        Set the nucleotide
        @param L1: Nucleotide 1
        @param L2: Nucleotide 1
        @return: Nucleotide
        """
        #val=self.UniValues()
        val=np.random.uniform()
        #Si la valeur est superieure a 0.5
        if(val>=0.5): return L1
        return L2
    
    def GenerateSeq(self, seqlength):
        """
        Generate a random sequence
        @return: a random sequence
        """
        #Get random Values
        values=self.GetValues(seqlength)
        #Create a list
        seq=list("Z"*seqlength)
        #Generate a sequence of SeqLength length
        for i in xrange(seqlength):
            #Generate AT
            if(values[i]>=self.GC): seq[i]=self.SetLetters("A","T")
            #Generate GC
            else: seq[i]=self.SetLetters("G","C")
        return "".join(seq)
    
    #Ecriture du fichier fasta
    def WriteFasta(self, fastafile, seq):
        """
          Write the fasta.
          @param fastafile: Name of the fasta file
          @param seq: Generated sequence
        """
        try:
            #Ouverture du fichier de compilation
            fasta = open(self.OutRep+fastafile+".fasta", "w")
            #Ecriture des donnees
            a=0
            for i in seq:
                fasta.write(">NC_%d_ranD\n"%(a))
                fasta.write(i[0]+"\n")
                a+=1
            #Fermeture du fichier
            fasta.close()
        except IOError:
            sys.exit("Error: We can not open the fasta file : %s"%(self.fastafile))
    
    def GenerateFile(self):
        """
        Generate a multifasta with random sequences
        """
        seq=[]
        values=self.GetIntValues(len(self.SeqLength),self.NbSeq)
        #Generate NbSeq sequences
        for i in xrange(self.NbSeq):
            seq+=[[self.GenerateSeq(self.SeqLength[values[i]])]]
        #Write the fasta
        self.WriteFasta("NC_ranD",seq)