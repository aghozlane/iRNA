"""
 @brief: Handle sequence from Genbank file 
 @author: Amine Ghozlane
 @version: 1.0 
"""
import string,sys,csv

class Sequences:
    
    def __init__(self):
        """
        Instanciate Sequence object
        """
        self.organism=""
        self.list_genes=[]
        self.DNA=""
    
    def parse(self,obj):
        """
        Get parsed data
        @param obj: Data object
        """
        self.organism,self.list_genes,self.DNA=obj.setdata(self.organism,self.list_genes,self.DNA)
    
    def invers_compl(self,gene):
        """
        Inverse and complemente gene sequence
        @param gene: Gene sequence
        """
        gene=string.replace(gene,'a','*')
        gene=string.replace(gene,'t','a')
        gene=string.replace(gene,'*','t')
        gene=string.replace(gene,'c','*')
        gene=string.replace(gene,'g','c')
        gene=string.replace(gene,'*','g')
        # A cette etape, la sequence est complemente
    
        gene=list(gene)
        gene.reverse()
        gene="".join(gene)
        # A cette etape, la sequence est inverse egalement
        return gene

    def fractionalGene(self,genedict,begin,end):
        """
        Get Gene between indicated set of position
        @param genedict: list of Gene dictionnary
        @return: Fractionnal gene sequence
        """
        return self.DNA[(genedict['begin']-1+begin):genedict['begin']-1+end]
    
    def completeGene(self,genedict,begin=None,end=None):
        """
        Get the complete Gene
        @param genedict: list of Gene dictionnary
        @return: Complete gene sequence
        """
        return self.DNA[(genedict['begin']-1):genedict['end']]
    
    def setGeneFunction(self,complete):
        """
        Define the getGene function to use
        @param complete: Complete flag
        @return: Function pointer
        """
        if(not complete): return self.fractionalGene
        return self.completeGene
    
    def setdefaultmRNAlen(self,begin,end,complete):
        """
        Compute default mRNA length
        @param begin: sRNA begin 
        @param end: sRNA end
        @param complete: Complete flag
        @return: RNA length
        """
        #Case: fractionnal gene
        if(not complete):
            #Compute length of mRNA
            if(begin<0 and end>0): return(abs(begin)+end)
            else: return(begin + end)
        #Case complete gene
        else: return None
    
    def writemRNA(self,begin,end,results,complete):
        """
        Write mRNA multifasta file
        @param begin: sRNA begin 
        @param end: sRNA end
        @param results: Path to result repertory
        @param complete: Complete flag
        """
        mRNA_multifastaName=results+"%s_mRNA.fasta"%self.organism
        #Set parameters depending on complete boolean
        getGene=self.setGeneFunction(complete)
        mRNAlen=self.setdefaultmRNAlen(begin,end,complete)  
        #Open multifasta
        try:
            multifastaWriter=open(mRNA_multifastaName,"w")
            for genedict in self.list_genes:
                #print(genedict)
                #Get the gene sequence
                gene=getGene(genedict,begin,end)
                # si ce gene est sur le brin complementraire, il est inverse pour etre dans le sens 5'-3' 
                if genedict['complement']!=None: gene=self.invers_compl(gene)
                #Check mRNA length if not complete gene
                genedict['genefound']=False
                if((len(gene)==mRNAlen or (complete and len(gene)>0))):
                    genedict['genefound']=True
                    #Write multifasta
                    if genedict['GeneID']!=None: multifastaWriter.write(">%s\n"%genedict['GeneID'])
                    else: multifastaWriter.write(">%s\n"%genedict['GeneName'])
                    multifastaWriter.write("%s\n"%gene)
                elif(complete and len(gene)<=0): 
                    if genedict['GeneID']!=None: print("mRNA \"%s\" is not found"%(genedict['GeneID']))
                    else: print("mRNA \"%s\" is not found"%(genedict['GeneName']))
                else: 
                    if genedict['GeneID']!=None: print("mRNA \"%s\" has not the proper length : %s %d!=%d (default)"%(genedict['GeneID'],gene,len(gene),mRNAlen))
                    else: print("mRNA \"%s\" has not the proper length : %s %d!=%d (default)"%(genedict['GeneName'],gene,len(gene),mRNAlen))
            multifastaWriter.close()
        except IOError: sys.exit("Error : can not open file %s"%mRNA_multifastaName)
        #except: sys.exit("Something went wrong with %s"%mRNA_multifastaName)
        
    def writesRNA(self,results,srna_data):
        """
        Write sRNA multifasta file
        @param results: Path to result repertory
        @param srna_data: sRNA list
        """
        sRNA_multifastaName=results+"%s_sRNA.fasta"%self.organism
        i=0         
        #Open multifasta
        try:
            multifastaWriter=open(sRNA_multifastaName,"w")
            for srnadict in srna_data.RNA_obj:
                gene=self.completeGene(srnadict)
                # si ce gene est sur le brin complementraire, il est inverse pour etre dans le sens 5'-3' 
                if srnadict['complement']: gene=self.invers_compl(gene)
                if(len(gene)>0):
                    #Write multifasta
                    if(srnadict['GeneName']!=None): 
                        multifastaWriter.write(">%s\n"%srnadict['GeneName'])
                        multifastaWriter.write("%s\n"%gene)
                    else: print("sRNA name missing line %d"%i)
                else:
                    if srnadict['complement']: print("sRNA \"%s\" at position %d - %d on antisense strand is not found"%(srnadict['GeneName'],srnadict['end'],(srnadict['begin']-1)))
                    else: print("sRNA \"%s\" at position %d - %d on sense strand is not found"%(srnadict['GeneName'],(srnadict['begin']-1),srnadict['end']))
                    
                i+=1
            multifastaWriter.close()
        except IOError: sys.exit("Error : can not open file %s"%sRNA_multifastaName)
        #except: sys.exit("Something went wrong with %s"%sRNA_multifastaName)
        
    def writeCorresponding(self,results):
        """
        Write corresponding geneID to genename for mRNA
        """
        correspondingName=results+"%s_corresponding_mRNAcode.txt"%self.organism
        try:
            with open(correspondingName, 'wt') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(["GeneID","gene"])
                for genedict in self.list_genes:
                    if(genedict['genefound'] and genedict['GeneID']!=None):
                        if (genedict['GeneName']): writer.writerow([genedict['GeneID'],genedict['GeneName']])
                        elif(genedict['locus_tag']): writer.writerow([genedict['GeneID'],genedict['locus_tag']])
                        else: writer.writerow([genedict['GeneID'],genedict['GeneID']])
                    elif(genedict['genefound'] and genedict['GeneID']==None): writer.writerow([genedict['GeneName'],genedict['GeneName']])
        except IOError: sys.exit("Error : can not open file %s"%correspondingName)
        #except: sys.exit("Something went wrong with %s"%correspondingName)    