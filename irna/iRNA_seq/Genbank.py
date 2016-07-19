"""
 @brief: Analyse and extract useful information from a genbank
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, re

class Genbank:
    
    def __init__(self,genbank_file):
        """
        Instanciate Genbank object
        """
        #Read Genbank file
        self.genbankString=self.readGenbank(genbank_file)
    
    def readGenbank(self,genbank_file):
        """
        Read Genbank file and save into a set of readlines.
        @param genbank_file: Path to Genbank file
        @return: List contening Genbank text
        """
        genbankString=None
        try: 
            genbankReader=open(genbank_file,'r')
            genbankString=genbankReader.readlines()
            genbankReader.close()
        except IOError: sys.exit("Error : can not open file %s"%genbank_file)
        except: sys.exit("Something went wrong with %s"%genbank_file)
        return genbankString
    
    def getGenes(self):
        """
        Get genes data
        @return: List of dictionnary with all genes data
        """
        posit=None
        GI=None
        locus_tag=None
        gene_name=None
        list_genes=[]
        #Regex of searched data
        regex_geneid=re.compile("^\s*/db_xref=\"GeneID:([0-9]+)\"")
        regex_position=re.compile("^\s*gene\s*(complement\()?([0-9]+)..([0-9]+)")
        regex_locus=re.compile("^\s*/locus_tag=\"(\S+)\"")
        regex_gi=re.compile("^\s*/db_xref=\"GI:([0-9]+)\"")
        regex_gene=re.compile("^\s*/gene=\"(\S+)\"")         
        for i in self.genbankString :
            #Cette expression reguliere permet de recuperer les bornes des genes et si ils sont sur le brin anti-sens ou le brin sens
            position_line=regex_position.search(i)
            #Ajout du dictionnaire a la liste de genes
            if position_line : posit=position_line 
            elif posit!=None:
                locus_line=regex_locus.search(i)
                GI_line=regex_gi.search(i)
                geneid_line=regex_geneid.search(i)
                gene_name_line=regex_gene.search(i)
                #Get the locus
                if locus_line: locus_tag=locus_line.group(1)
                #Get the GI
                elif GI_line: GI= int(GI_line.group(1))
                elif gene_name_line: gene_name=gene_name_line.group(1)
                #Get the Geneid 
                elif geneid_line:
                    list_genes.append({'begin':int(posit.group(2)),'end':int(posit.group(3)),'complement':posit.group(1),'GeneID':int(geneid_line.group(1)),'GeneName':gene_name,"GI":GI,"locus_tag":locus_tag.strip().lower()})
                    posit=None
                    GI=None
                    locus_tag=None
                    gene_name=None
        return list_genes
    
    def getDNA(self):
        """
        Read DNA sequence and transform into string
        @param nom_fichier: File
        @param ADN: DNA string
        """
        DNA=""
        regex=re.compile("^\s+[0-9]+\s+([atcg ACTG]*)$")
        for i in self.genbankString :
            # Expression reguliere qui permet de recuperer seulement les lignes contenant la sequence d'ADN
            ligne=regex.search(i) 
            if ligne :
                if ligne.group(1)!=None:
                    Tmp=str.split(ligne.group(1))
                    for j in Tmp:
                        DNA+=j
        DNA=DNA.lower()
        return DNA
    
    def getOrganism(self):
        """
        Get organism name
        """
        regex=re.compile("ORGANISM\s+(.+)")
        for i in self.genbankString :
            #Recuperation du nom d'organisme
            a=regex.search(i)
            if(a): return "_".join(a.group(1).split())
                    
    def setdata(self,organism,list_genes,DNA):
        """
        Add genbank information
        @param organism: 
        @param list_genes: 
        @param DNA: 
        @return: 
        """
        if(self.genbankString!=None):
            #Get DNA sequence
            DNA=self.getDNA()
            #Get Genes information
            list_genes=self.getGenes()
            #Get Organism
            organism=self.getOrganism()
        return organism,list_genes,DNA