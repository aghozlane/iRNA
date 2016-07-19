"""
 @brief: Extraction of multifasta
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, re, os, string
from Stat import *


class Fasta(Stat):
    """
    @brief: Write fasta needed for the different software from a multifasta
    """
        
    def __init__(self):
        """
        Build fasta object
        """
        Stat.__init__(self)
        "@note cacheFasta: Fasta sequence" 
        self.cacheFasta = ""
        "@note: fastafile: Name of the fastafile"
        self.fastafile = None
    
    #Ecriture du fichier fasta
    def WriteFasta(self):
        """
        Write the fasta
        """
        try:
            #Ouverture du fichier de compilation
            fasta = open(self.fastafile, "w")
            #Ecriture des donnees
            fasta.write(self.cacheFasta)
            #Fermeture du fichier
            fasta.close()
        except IOError:
            sys.exit("Error: We can not open the fasta file : %s"%(self.fastafile))
    
    def EmptyRep(self, fasta_out):
        """
        Empty the folder of fasta files
        @param fasta_out: Fasta repertory
        """
        #Regex
        regex_fasta=re.compile(".*.fasta")
        #Get the list of file the repertory
        files = os.listdir(fasta_out)
        for file in files:
            a=regex_fasta.match(file)
            #matched file
            if a:
                #Path to the file
                path=fasta_out+file
                #Remove the file
                os.remove(path)
    
    #Extraction des donnees du fichier multifasta
    def ExtractFasta(self, mf, fasta_out):
    #def ExtractFasta(self, mf, fasta_out, dbmanage):
        """
        Extract multifasta to fasta.
        @param mf: Multifasta file
        @param fasta_out: Fasta repertory
        """
        self.GC=0.0
        self.LenSeq=0.0
        self.NbSeq=0.0
        self.LenRan=[]
        self.out=fasta_out+"seq_inf.txt"
        #Open seqinf file
        self.openSeqInf()
        #Vider le fasta out de ses fasta
        self.EmptyRep(fasta_out)
        #Definition de la regex
        regex = re.compile("^>(\S+)\s*\S*")
        regex_seq=re.compile("^[ATGCatgc]+\S*")
        try:
            #Ouverture du fichier de compilation
            multifasta = open(mf, "r")
            #Lecture de tous les fichiers
            cacheData = False 
            for i in multifasta:
                a = regex.match(i)
                if a:
                    seqname=a.group(1)
                    #Ecriture du fichier fasta
                    if(self.cacheFasta):
                        self.WriteFasta()
                    #Definition du nom du fichier 
                    self.fastafile=fasta_out+seqname+".fasta"
                    cacheData=True
                    self.cacheFasta=">%s\n"%seqname
                #Mise en cache de la sequence du fasta
                elif(cacheData and regex_seq.match(i)):
                    self.cacheFasta+=i
                    #self.ComputeGC(seqname,i.translate(None,'\n'),dbmanage)                    
                    self.ComputeGC(seqname,i.translate(None,'\n'))
            #Ecriture du dernier fichier fasta
            if(self.cacheFasta):
                self.WriteFasta()               
            #Fermeture du fichier
            multifasta.close()
            print("GC=%f Mean_len=%f"%(self.GetGC(),self.GetLenSeq()))
        #Erreur d'ouverture de fichier
        except IOError: sys.exit("Error: We can not open the multifasta file : %s"%mf)
        #Erreur inattendue
        except : sys.exit("Something went wrong with multifasta file : %s"%mf)
        #Close seq inf file
        self.CloseSeqInf()
        return (self.GetGC(),self.LenRan)