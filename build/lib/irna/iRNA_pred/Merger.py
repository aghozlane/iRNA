"""
 @brief: Merge the xml files
 @author: Amine Ghozlane
 @version: 1.0 
"""
import os,re,sys, bisect,csv
from lxml import etree

class Merger:
    """
    Merge the xml files in a sqlite database
    """
    #Resulting tree
    result_data=etree.Element("result")

    def __init__(self, repIn):
        """
        Initiate the constructor
        @param repIn: Result repertory
        """
        #Get repertory
        self.repIn=repIn
        #Get xml files
        self.files=os.listdir("%s"%self.repIn)
        self.files.sort()
        self.softprevious=[None]*3     
    
                     
    def buildTree(self, tree, dbmanage):
        """
        Read the xml tree and copy the data into the sqlite database
        @param tree: xml tree parsed
        @param dbmanage: Database manager
        """      
        for soft in tree:
            soft_name=soft.get("id")
            soft_ref=soft.get("reference")
            if(self.softprevious[0]==soft_name and self.softprevious[1]==soft_ref):
                softid=self.softprevious[2]
            else:
                softid=dbmanage.setSoft(soft_name)
                self.softprevious[0]=soft_name
                self.softprevious[1]=soft_ref
                self.softprevious[2]=softid
            for interaction in list(soft):
                contact=[]
                #Add interaction
                sRNAid=dbmanage.getRNA(interaction.get("sRNA"),0)
                mRNAid=dbmanage.getRNA(interaction.get("mRNA"),1)
                #interactid=None
                if(sRNAid and mRNAid): interactid=dbmanage.setInteract(sRNAid,mRNAid,softid)
                else: sys.exit("RNA unknown in the results")
                #Add contact
                for i in list(interaction):
                    contact+=[[int(i.get("sRNA_begin")),int(i.get("sRNA_end")),int(i.get("mRNA_begin")),int(i.get("mRNA_end")),float(i.get("score"))]]
                dbmanage.setContact(interactid,contact)
                del(contact)
                    
    def getData(self, file):
        """
        Parse xml tree
        @param file: Result part in Xml format
        @return: resulting tree
        """
        tree=None
        #Parse tree
        try:
            parser= etree.XMLParser(remove_blank_text=True)
            tree = etree.parse(file,parser).getroot()
        except: pass
        return(tree)
    
    def setRNA(self,dbmanage,RNA_inf,type):
        """
        Parse RNA_inf file
        @param RNA_inf: Path to RNA_inf file
        @param dbmanage: Access to the database
        @param type: Type of RNA
        """
        try:
            RNA_infReader=csv.reader(open(RNA_inf, 'rt'), delimiter='\t')
            #passage de l'entete
            RNA_infReader.next()
            RNAtab=[[i[0],int(i[1]),type] for i in RNA_infReader]
            dbmanage.setRNA(RNAtab)
        except IOError: sys.exit("Error : can not open file %s"%RNA_inf)
            
    def merge(self,dbmanage):
        """
        Merge iRNA xml files into a sqlite database
        @param dbmanage: Database manager
        """
        tree=None
        #Trouver les fichiers xml
        #^(?!iRNA_result\.xml)
        regex=re.compile(".*\.xml")
        #Lecture de la liste de fichier
        for i in self.files:
            a=regex.match(i)
            if(a):
                #print(i)
                #Path to the file
                filepath=self.repIn+i
                #Get the tree
                tree=self.getData(filepath)
                #Copy tree
                if(tree is not None): self.buildTree(tree,dbmanage)     