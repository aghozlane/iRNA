"""
 @brief: Parse soft output
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sys, re
from lxml import etree

class Parse:
    """
    @brief: Parse data from the different software.
    """
    def parseGuugle(self, data, matrix, root):
        """
        Parse Guugle result
        @param data: result
        @param matrix: score matrix
        @param root: xml root node
        """
        #result=Interaction([sRNA, mRNA],matrix)
        #result=Interaction([mRNA,sRNA],matrix)
        regex=re.compile("MatchLength:\s+(\S+)\s+\S+\s+at\s+(\S+)\s+vs.\s+\S+\s+at\s+(\S+)")
        #Read data
        seqs=[]
        flag=False
        j=0
        nb=0
        for i in data.split("\n"):
            a=regex.match(i)
            if a:
                #debut mRNA fin mRNA debut sRNA fin sRNA score
                longueurRNA= int(a.group(1))
                #Add data
                #tab=[int(a.group(3)),int(a.group(3))+longueurRNA,int(a.group(2)),int(a.group(2))+longueurRNA]
                tab=[int(a.group(2)),int(a.group(2))+(longueurRNA-1),int(a.group(3)),int(a.group(3))+(longueurRNA-1)]
                #result.addTable("%d\t%d\t%d\t%d\t"%(int(a.group(3)),int(a.group(3))+longueurRNA,int(a.group(2)),int(a.group(2))+longueurRNA))
                #result.addTable([int(a.group(3)),int(a.group(3))+longueurRNA,int(a.group(2)),int(a.group(2))+longueurRNA])
                #It will take the next two lines
                j=2
                flag=True
            #Lecture des deux morceaux de sequence
            elif j>0:
                #seq = [mRNA, sRNA]
                seqs+=[(i.strip("35\n")).upper()]
                j-=1
            #Ajout de la sequence
            if j==0 and flag and seqs:
                #result.addTable("%d\n"%result.scoreSequence(seqs))
                tab+=[root.scoreSequence(seqs)]
                if(len(tab)==5): root.addComparison(tab)
                else : sys.exit("Probleme with the regex Guugle\n%s"%tab)
                nb+=1
                #Reinit
                seqs=[]
                flag=False
    
#    def parseIntaRNA(self, data, matrix, root):
#        """
#        Parse IntaRNA result
#        @param data: Soft result
#        @param matrix: score matrix
#        @param root: xml root node
#        """
#        tab=[]
#        #Definition des regex
#        position_mRNA=re.compile("^positions\s+with\s+dangle\(target\):\s+(\S+)\s+--\s+(\S+)\s*")
#        position_sRNA=re.compile("^positions.+with.+dangle\(miRNA.+\):\s+(\S+)\s+--\s+(\S+)\s*")
#        energy=re.compile("^energy:\s*(\S+)\s*kcal/mol\s*")
#        for i in data.split("\n"):
#            a=position_mRNA.match(i)
#            b=position_sRNA.match(i)
#            c=energy.match(i)
#            #debut mRNA fin mRNA debut sRNA fin sRNA score
#            #mRNA
#            if a:
#                tab=[int(a.group(1)),int(a.group(2))]
#                flag_a=True
#            #sRNA
#            elif b:
#                tab+=[int(b.group(1)),int(b.group(2))]
#                flag_b=True
#            #energie
#            elif c:
#                tab+=[float(c.group(1))]
#                #Ajout de la comparison
#                try:
#                    if(len(tab)==5 and flag_a and flag_b): 
#                        root.addComparison(tab)
#                        flag_a=False
#                        flag_b=False
#                        #Reinitialisation
#                        del(tab)
#                        tab=[]
#                    else :
#                        sys.exit("Probleme with the regex IntaRNA\n%s"%tab)
#                except:
#                    sys.exit("Something went wrong with parse IntaRNA")

    def parseIntaRNA(self, data, matrix, root):
        """
        Parse IntaRNA result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        tab=[]
        #Definition des regex
        position_mRNA=re.compile("^positions\s+with\s+dangle\(target\):\s+(\S+)\s+--\s+(\S+)\s*")
        position_sRNA=re.compile("^positions.+with.+dangle\(ncRNA\):\s+(\S+)\s+--\s+(\S+)\s*")
        energy=re.compile("^energy:\s*(\S+)\s*kcal/mol\s*")
        for i in data.split("\n"):
            a=position_mRNA.match(i)
            b=position_sRNA.match(i)
            c=energy.match(i)
            #debut mRNA fin mRNA debut sRNA fin sRNA score
            #mRNA
            if a:
                tab=[int(a.group(1)),int(a.group(2))]
                flag_a=True
            #sRNA
            elif b:
                tab+=[int(b.group(1)),int(b.group(2))]
                flag_b=True
            #energie
            elif c:
                tab+=[float(c.group(1))]
                #Ajout de la comparison
                try:
                    if(len(tab)==5 and flag_a and flag_b): 
                        root.addComparison(tab)
                        flag_a=False
                        flag_b=False
                        #Reinitialisation
                        del(tab)
                        tab=[]
                    else :
                        sys.exit("Probleme with the regex IntaRNA\n%s"%tab)
                except:
                    sys.exit("Something went wrong with parse IntaRNA")                   
    
    def parseRNAup(self, data, matrix, root):
        """
        Parse RNAup result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        tab=[]
        #Definition de la regex
        regex=re.compile("\S+&\S+\s+(\S+),(\S+)\s+:\s+(\S+),(\S+)\s+\((\S+)\s+=.+\)\s*")
        #Parsing
        for i in data.split("\n"):
            a=regex.match(i)
            if a:
                mRNA_begin=int(a.group(1))
                mRNA_end=int(a.group(2))
                sRNA_begin=int(a.group(3))
                sRNA_end=int(a.group(4))
                energy=float(a.group(5))
                if(mRNA_begin!=0 or mRNA_end!=0 or sRNA_begin!=0 or sRNA_end!=0 or energy!=100.0):
                    #debut mRNA fin mRNA debut sRNA fin sRNA score
                    tab+=[mRNA_begin,mRNA_end,sRNA_begin,sRNA_end,float(a.group(5))]
                    #Ajout du resultat
                    if(len(tab)==5): root.addComparison(tab)
                    else :
                        sys.exit("Probleme with the regex RNAup\n%s"%tab)                     
    
    def parseRNAplex(self, data, matrix, root):
        """
        Parse RNAplex and RNAduplex result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        #Definition de la regex
        regex=re.compile("\S+&\S+\s+(\S+),(\S+)\s+:\s+(\S+),(\S+)\s+\(\s*(\S+)\)\s*")
        for i in data.split("\n"):
            a=regex.match(i)
            if a:
                #debut mRNA fin mRNA debut sRNA fin sRNA score
                tab=[int(a.group(1)),int(a.group(2)),int(a.group(3)),int(a.group(4)),float(a.group(5))]
                #Ajout du resultat
                if(len(tab)==5): root.addComparison(tab)
                else :
                    sys.exit("Probleme with the regex RNAplex\n%s"%tab)
    
    def parseRNAduplex(self, data, matrix, root):
        """
        Parse RNAplex and RNAduplex result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        #Definition de la regex
        regex=re.compile("\S+&\S+\s+(\S+),(\S+)\s+:\s+(\S+),(\S+)\s+\(\s*(\S+)\)\s*")
        for i in data.split("\n"):
            a=regex.match(i)
            if a:
                #debut mRNA fin mRNA debut sRNA fin sRNA score
                tab=[int(a.group(3)),int(a.group(4)),int(a.group(1)),int(a.group(2)),float(a.group(5))]
                #Ajout du resultat
                if(len(tab)==5): root.addComparison(tab)
                else :
                    sys.exit("Probleme with the regex RNAplex\n%s"%tab)
    
    def getposition_RNAcofold0(self, result):
        """
        Get the position of interaction for RNAcofold
        @param result: sRNA result
        @return: tabmRNA, tabsRNA
        """
        tab=[]
        for i in xrange(len(result)):
            if(result[i]=="("):
                tab.append(i+1)
            elif(result[i]==")"):
                tab.pop()
            i-=1
        return tab
    
    def getposition_RNAcofold1(self, result):
        """
        Get the position of interaction for RNAcofold
        @param result: mRNA result 
        @return: tabmRNA, tabsRNA
        """
        tab=[]
        i=len(result)-1
        while(i>-1):
            #Ajout dans la pile
            if(result[i]==")"):
                #Ajout de la position
                tab.append(i+1)
            #On retire de la pile
            elif(result[i]=="("):
                tab.pop()
            i-=1
        return tab  
    
    def reverseElongate(self,tab):
        """
        Get a double reverse for a tab of tab
        @param tab: Table of table
        @return: Reversed table of table
        """
        #Reverse du premier tableau
        tab.reverse()
        for i in tab:
            #reverse les sous tableau
            i.reverse()
        return tab
    
    def forwardTable(self, tab_RNA):
        """
        Reverse a table
        @param tab_RNA: Table of RNA
        @return: Reversed table
        """
        maxValue=tab_RNA[0][0]
        for i in tab_RNA:
            i[0]=maxValue-i[0]+1
            i[1]=maxValue-i[1]+1
        return tab_RNA,maxValue
    
    def parseRNAcofold(self, data, matrix, root):
        """
        Parse RNAcofold result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        tag=True
        #tag=False
        result=[]
        energy=0
        #Regex
        regex_result=re.compile("(^[.\[\]\(\)\&\{\},\|]+)")
        #regex_energy=re.compile(".*delta G binding=\s*(\S+)")
        regex_energy=re.compile("(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s*")
        #Lecture du fichier
        for i in data.split("\n"):
            a=regex_result.match(i)
            b=regex_energy.match(i)
            if(a):
                #On prend la premiere sequence de resultat
                if(tag):
                    #sRNA mRNA
                    result=re.split("&",a.group(1))
                tag=False
                #tag=True
            elif(b):
                energy=b.group(1)
        #Result
        if(result):
            #Detection des positions
            tab_sRNA=self.getposition_RNAcofold0(result[0])
            tab_mRNA=self.getposition_RNAcofold1(result[1])
            #Si des positions sont detectees
            if(tab_sRNA and tab_mRNA):
                #On reverse le resultat de mRNA pour le traiter facilement
                tab_mRNA.reverse()
                #Elongation des positions 
                tab_sRNA=self.elongatePositions_rip(tab_sRNA)
                tab_mRNA=self.elongatePositions_rip(tab_mRNA)
                #Reversion pour l'elongate
                tab_mRNA=self.reverseElongate(tab_mRNA)
                #Reverse the table
                tab_mRNA,maxValue=self.forwardTable(tab_mRNA)
                #Find common position
                (tab_sRNA,tab_mRNA)=self.consensus([],[],tab_sRNA, tab_mRNA, tab_sRNA[0][0], tab_mRNA[0][0], 0, 0)
                #Reverse again the table
                tab_mRNA=self.reverseTablemax(tab_mRNA,maxValue)
                #Add comparison
                for i in range(len(tab_mRNA)):
                    #debut mRNA fin mRNA debut sRNA fin sRNA score
                    root.addComparison([tab_mRNA[i][1],tab_mRNA[i][0],tab_sRNA[i][0],tab_sRNA[i][1],energy])
            
    def getResultSequence_RNAhybrid(self, out_RNA, in_RNA,  RNA_begin, length_RNA, int_char):
        """
        Get contact sequence.
        @param out_RNA: Non contact nucleotid
        @param in_RNA: Contact nucleotid
        @param RNA_begin: Begining position of the RNA
        @param length_RNA: Length of the RNA sequence
        @param int_char: Interaction caracter
        """
        RNA_R=list("."*length_RNA)
        a=0
        i=RNA_begin
        while( a<len(in_RNA)):
            #Point de contact
            if(in_RNA[a]=="A" or in_RNA[a]=="U" or in_RNA[a]=="G" or in_RNA[a]=="C"):
                RNA_R[i]=int_char
                a+=1
                i+=1
            #Point de non contact
            elif(out_RNA[a]=="A" or out_RNA[a]=="U" or out_RNA[a]=="G" or out_RNA[a]=="C"):
                a+=1
                i+=1
            #Cas d'un saut
            else:
                a+=1
        return "".join(RNA_R)
  
    def parseRNAhybrid(self, data, matrix, root):
        """
        Parse RNAhybrid result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        mRNA_flag=False
        sRNA_flag=False
        length_sRNA_flag=True
        energy=0
        #Definition des regex
        length_mRNA_regex=re.compile("length:\s+(\S+)")
        mRNA_begin_regex=re.compile("position\s+(\S+)")
        mRNA_regex=re.compile("target\s+5' (.+) 3'")
        sRNA_regex=re.compile("miRNA\s+3' (.+) 5'")
        energy_regex=re.compile("mfe:\s+(\S+)\s+kcal/mol")
        #Lecture du fichier
        for i in data.split("\n"):
            a = mRNA_regex.match(i)
            b = sRNA_regex.match(i)
            c = mRNA_begin_regex.match(i)
            d = length_mRNA_regex.match(i)
            e = energy_regex.match(i)
            #Outside mRNA
            if(a):
                out_mRNA=a.group(1)
                mRNA_flag=True
            #outside sRNA
            elif(b):
                out_sRNA=b.group(1)
            #mRNA begin
            elif(c):
                mRNA_begin=int(c.group(1))
            #Length mRNA
            elif (d and length_sRNA_flag):
                length_mRNA=int(d.group(1))
                length_sRNA_flag=False
            #Energy
            elif(e):
                energy=e.group(1)
            #Length sRNA                
            elif(d and not length_sRNA_flag):
                length_sRNA=int(d.group(1))
            #Contact mRNA
            elif(mRNA_flag):
                in_mRNA=i[10:-3]
                sRNA_flag=True
                mRNA_flag=False
            #Contact sRNA
            elif (sRNA_flag):
                in_sRNA=i[10:-3]
                sRNA_flag=False
        #creation de la sequence mRNA_R 5'-3'
        mRNA_r=self.getResultSequence_RNAhybrid(out_mRNA, in_mRNA, mRNA_begin-1, length_mRNA, "[")
        sRNA_r=self.getResultSequence_RNAhybrid(out_sRNA, in_sRNA, 0, length_sRNA, "]")
        #Get position
        tab_mRNA=self.getPositions_rip(mRNA_r,"[")
        tab_sRNA=self.getPositions_rip(sRNA_r,"]")
        #Si des positions sont identifies
        if(tab_mRNA and tab_sRNA):
            #Elonguate
            tab_mRNA=self.elongatePositions_rip(tab_mRNA)
            tab_sRNA=self.elongatePositions_rip(tab_sRNA)
            #Find common position
            (tab_sRNA,tab_mRNA)=self.consensus([],[],tab_sRNA, tab_mRNA, tab_sRNA[0][0], tab_mRNA[0][0], 0, 0)
            #Reverse the table
            tab_sRNA=self.reverseTablemax(tab_sRNA,length_sRNA)
            #Add comparison
            for i in range(len(tab_mRNA)):
                #debut mRNA fin mRNA debut sRNA fin sRNA score
                root.addComparison([tab_mRNA[i][0],tab_mRNA[i][1],tab_sRNA[i][1],tab_sRNA[i][0],energy])
    
    def getPositions_rip(self, RNA_r, seq):
        """
        Get the positions for rip problem output
        @param RNA_r: RNA result output 
        @return: Position of contact for RNA   
        """
        tab_RNA=[]
        #Find contact for RNA
        for i in range(len(RNA_r)):
            if(RNA_r[i]==seq):
                tab_RNA.append((i+1))
        return tab_RNA
    
    def elongatePositions_rip(self, tab):
        """
        Elongate the contact position
        @param tab: Table of contact position
        @return: Table elongate
        """
        begin=tab[0]
        end=tab[0]
        tab_r=[]
        #Elegation du tableau
        for i in tab[1:]:
            #Ecart superieur a 1
            if((end+1)!=i):
                tab_r+=[[begin,end]]
                begin=i
                end=i
            #Elongation
            else:
                end+=1
        #Last position
        tab_r+=[[begin,end]]
        return tab_r
    
    def getSequences_53(self, RNA, tab_RNA):
        """
        Get the sequence in contact in 5'3' sense
        @param RNA: RNA sequence
        @param tab_RNA: Position of contact for RNA
        @return: RNA tab of contact sequences 
        """
        seq_RNA=[]
        #Get RNA sequence
        for i in tab_RNA:
            #On compense le debut a 0 du tableau
            seq_RNA+=[RNA[(i[0]-1):i[1]]]
        return seq_RNA
    
    def getSequences_local_reversed(self, RNA):
        """
        Get the sequence in contact in 5'3' sense
        @param RNA: RNA sequence
        @return: RNA tab of contact sequences 
        """
        seq_RNA=[]
        #Get RNA sequence
        for i in RNA:
            #On compense le debut a 0 du tableau
            seq_RNA+=[self.reverseSequence(i)]
        return seq_RNA
    
    def consensus(self,resultsRNA,resultmRNA, tab_sRNA, tab_mRNA, debutsRNA, debutmRNA, i, j):
        """
        Find the contact position between sRNA and mRNA
        @param resultsRNA: New table of contact for sRNA
        @param resultmRNA: New table of contact for mRNA
        @param tab_sRNA: Old table of contact for sRNA
        @param tab_mRNA: Old table of contact for mRNA
        @param debutsRNA: sRNA begin
        @param debutmRNA: mRNA begin
        @param i: Position in the old table of sRNA
        @param j: Position in the old table of mRNA
        """
        #Calcul de la lon
        lensRNA=(tab_sRNA[i][1]-debutsRNA)
        lenmRNA=(tab_mRNA[j][1]-debutmRNA)
        #sRNA est le plus long
        if(lensRNA>lenmRNA):
            reste=lenmRNA
        #sRNA est le plus court
        elif(lensRNA<lenmRNA):
            reste=-lensRNA
        #Ils ont la meme taille
        else:
            reste=0   
        #Il y a un reste positif, il faut couper le mRNA
        if(reste<0):
            #casser tabmRNA
            finmRNA=debutmRNA+abs(reste)
            resultmRNA+=[[debutmRNA,finmRNA]]
            resultsRNA+=[[debutsRNA,tab_sRNA[i][1]]]
            #On avance au sRNA suivant
            i+=1
            #On avance le debut du mRNA
            (resultsRNA,resultmRNA)=self.consensus(resultsRNA, resultmRNA, tab_sRNA, tab_mRNA, tab_sRNA[i][0], finmRNA+1,i,j)
        #Il y a un reste negatif, il faut couper le sRNA 
        elif(reste>0):
            #casser tabsRNA
            finsRNA=debutsRNA+abs(reste)
            resultsRNA+=[[debutsRNA,finsRNA]]
            resultmRNA+=[[debutmRNA,tab_mRNA[j][1]]]
            #On avance au mRNA suivant
            j+=1
            #On avance le debut du sRNA
            (resultsRNA,resultmRNA)=self.consensus(resultsRNA, resultmRNA, tab_sRNA, tab_mRNA,  finsRNA+1, tab_mRNA[j][0],i,j)
        #Les deux ont la meme longueur ou l'un a une taille 0
        else:
            #Le sRNA a une taille 0
            if(lensRNA==0 and lenmRNA>0):
                resultsRNA+=[[debutsRNA,tab_sRNA[i][1]]]
                #N'avance que d'un le mRNA
                resultmRNA+=[[debutmRNA,debutmRNA]]
                #on avance au sRNA suivant
                i+=1
                (resultsRNA,resultmRNA)=self.consensus(resultsRNA,resultmRNA,tab_sRNA, tab_mRNA,  tab_sRNA[i][0],debutmRNA+1,i,j)
            #Le mRNA a une taille 0
            elif(lenmRNA==0 and lensRNA>0):
                #N'avance que d'un le sRNA
                resultsRNA+=[[debutsRNA,debutsRNA]]
                resultmRNA+=[[debutmRNA,tab_mRNA[j][1]]]
                #On avance au mRNA suivant
                j+=1
                #On avance le sRNA
                (resultsRNA,resultmRNA)=self.consensus(resultsRNA,resultmRNA,tab_sRNA, tab_mRNA,  debutsRNA+1,tab_mRNA[j][0],i,j)
            #Les deux ont une taille 0
            else:
                resultsRNA+=[[debutsRNA,tab_sRNA[i][1]]]
                resultmRNA+=[[debutmRNA,tab_mRNA[j][1]]]
                #On avance au mRNA et sRNA suivant
                i+=1
                j+=1
                #Seulement si on est pas arrive au bout du tablea pour les deux
                if(i<len(tab_sRNA) or j<len(tab_mRNA)):
                    (resultsRNA,resultmRNA)=self.consensus(resultsRNA,resultmRNA,tab_sRNA, tab_mRNA,tab_sRNA[i][0],tab_mRNA[j][0],i,j)
        return (resultsRNA, resultmRNA)
            
    
    def reverseSequence(self, seq):
        """
        Reverse a sequence
        @param seq: Sequence of text
        @return: Reversed sequence
        """
        listseq=list(seq)
        listseq.reverse()
        return "".join(listseq)
    
    def reverseTablemax(self, tab_RNA, maxValue):
        """
        Reverse a table
        @param tab_RNA: Table of RNA
        @return: Reversed table
        """
        for i in tab_RNA:
            i[0]=maxValue-i[0]+1
            i[1]=maxValue-i[1]+1
        return tab_RNA
    
    def reverseMinTable(self, tab_RNA):
        """
        Reverse a table into 5'3'
        @param tab_RNA: Table of RNA
        @return: Reversed table
        """
        maxValue=tab_RNA[-1][1]
        for i in tab_RNA:
            i[0]=maxValue-i[0]+1
            i[1]=maxValue-i[1]+1
        return tab_RNA
    
    def parseRactip(self, data, matrix, root):
        """
        Parse Ractip result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        flag_seq=True
        flag_res=True
        sRNA_r=""
        mRNA_r=""
        mRNA=""
        sRNA=""
        #Definition des regex
        regex_sequence=re.compile("(^[ATGCatgc]+(?!LPK))")
        regex_result=re.compile("(^[.\[\]\(\)]+)")
        #Lecture du fichier
        for i in data.split("\n"):
            a=regex_sequence.match(i)
            b=regex_result.match(i)
            #Sequence
            if a:
                if(flag_seq):
                    mRNA=a.group(1).upper()
                    flag_seq=False
                else:
                    sRNA=a.group(1).upper()
            #Score
            elif b:
                if(flag_res):
                    mRNA_r=b.group(1)
                    flag_res=False
                else:
                    sRNA_r=b.group(1)
        #Compute position
        maxValue=len(sRNA)
        sRNA=self.reverseSequence(sRNA)
        sRNA_r=self.reverseSequence(sRNA_r)
        #Get position
        tab_mRNA=self.getPositions_rip(mRNA_r,"[")
        tab_sRNA=self.getPositions_rip(sRNA_r,"]")
        #Si des positions sont identifies
        if(tab_mRNA and tab_sRNA):
            #Elonguate
            tab_mRNA=self.elongatePositions_rip(tab_mRNA)
            tab_sRNA=self.elongatePositions_rip(tab_sRNA)
            #Find common position
            (tab_sRNA,tab_mRNA)=self.consensus([],[],tab_sRNA, tab_mRNA, tab_sRNA[0][0], tab_mRNA[0][0], 0, 0)
            #Find sequence
            seq_mRNA=self.getSequences_53(mRNA, tab_mRNA)
            seq_sRNA=self.getSequences_53(sRNA, tab_sRNA)
            #Reverse sequence into 3' 5'
            #tab_sRNA=self.reverseTablemax(tab_sRNA)
            tab_sRNA=self.reverseTablemax(tab_sRNA, maxValue)
            #Add comparison
            for i in range(len(tab_mRNA)):
                #debut mRNA fin mRNA debut sRNA fin sRNA score
                root.addComparison([tab_mRNA[i][0],tab_mRNA[i][1],tab_sRNA[i][1],tab_sRNA[i][0],root.scoreSequence([seq_mRNA[i], seq_sRNA[i]])])
    
    def readCompfile(self, complist_file):
        """
        Read the comp_list file
        @param complist_file: Comparison list file
        @return: Comparison data
        """
        comp_data=[]
        regex=re.compile("(\S+\.fasta)\s+(\S+\.fasta)")
        try:
            comp_file=open(complist_file,"r")
            for i in comp_file:
                a=regex.match(i)
                if a:
                    comp_data+=[[a.group(1),a.group(2)]]
            comp_file.close()
        #Read error
        except IOError:
            sys.exit("Error : Can not open %s"%complist_file)
        #Other error
        except :
            sys.exit("Something went wrong with %s"%complist_file)
            
        return comp_data
    
    def comparaison(self, mRNA, sRNA, comp_data):
        """
        Determine if comparison is needed
        @param mRNA: mRNA name
        @param sRNA: sRNA name
        @param comp_data: Comparison data
        @return: True if comparison is expected, False otherwise
        """
        for i in comp_data:
            #La comparaison est demandee
            if(sRNA==i[0][:-6] and mRNA==i[1][:-6]):
                return True
        #La comparaison n'est pas attendue
        return False
    
    def parseBlastall(self, data, matrix, root, complist_file):
        """
        Parse Blastall result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        @param complist_file: Comparison list file
        """
        regex=re.compile("(\S+)\s+"*11+"(\S+)\s*")
        if(complist_file):
            comp_data=self.readCompfile(complist_file)
            #Trie des donnees pour eviter les redondances    
            data=data.split("\n")
            data.sort()
            #Lecture du fichier
            for i in data:
                a=regex.match(i)
                if a:
                    #mRNA sRNA
                    if(self.comparaison(a.group(2), a.group(1), comp_data)):
                        #si la comparaison en question est differente de la derniere
                        if(root.getlastmRNA()!= a.group(2) or root.getlastsRNA()!=a.group(1)):
                            #Ajout du nouveau noeud
                            #mRNA sRNA
                            root.addListOfComparison(a.group(2), a.group(1))
                        #debut mRNA fin mRNA debut sRNA fin sRNA score
                        root.addComparison([int(a.group(10)),int(a.group(9)),int(a.group(7)),int(a.group(8)),float(a.group(11))])
        else:
            #Trie des donnees pour eviter les redondances    
            data=data.split("\n")
            data.sort()
            #Lecture du fichier
            for i in data:
                a=regex.match(i)
                if a:
                   
                    #si la comparaison en question est differente de la derniere
                    if(root.getlastmRNA()!= a.group(2) or root.getlastsRNA()!=a.group(1)):
                        #Ajout du nouveau noeud
                        #mRNA sRNA
                        root.addListOfComparison(a.group(2), a.group(1))                    
                    else:
                        #debut mRNA fin mRNA debut sRNA fin sRNA score
                        root.addComparison([int(a.group(10)),int(a.group(9)),int(a.group(7)),int(a.group(8)),float(a.group(12))])
    
    def parseBistaRNA(self, data, matrix, root):
        """
        Parse BistaRNA result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        flag_seq=True
        flag_res=True
        sRNA_r=""
        mRNA_r=""
        mRNA=""
        sRNA=""
        #Definition des regex
        #(^[ATGCatgc]+)
        regex_sequence=re.compile("^([ATGCatgc]+)")
        bad_antisense=re.compile("Antisense")
        bad_target=re.compile("Target")
        regex_result=re.compile("(^[\.\,\[\]\(\)]+)")
        #Lecture du fichier
        for i in data.split("\n"):
            a=regex_sequence.match(i)
            b=regex_result.match(i)
            c=bad_antisense.match(i)
            d=bad_target.match(i)
            #Sequence
            if a and not c and not d:
                if(flag_seq):
                    sRNA=a.group(1).upper()
                    flag_seq=False
                else:
                    mRNA=a.group(1)
            #Score
            elif b:
                if(flag_res):
                    sRNA_r=b.group(1)
                    flag_res=False
                else:
                    mRNA_r=b.group(1)
        #Compute position
        tab_mRNA=self.getPositions_rip(mRNA_r,"]")
        tab_sRNA=self.getPositions_rip(sRNA_r,"[")
        #Si des positions sont identifies
        if(tab_mRNA and tab_sRNA):
            #Elongate
            tab_mRNA=self.elongatePositions_rip(tab_mRNA)
            tab_sRNA=self.elongatePositions_rip(tab_sRNA)
            #Find common position
            (tab_sRNA,tab_mRNA)=self.consensus([],[],tab_sRNA, tab_mRNA, tab_sRNA[0][0], tab_mRNA[0][0], 0, 0)
            #Find sequence
            seq_mRNA=self.getSequences_53(mRNA, tab_mRNA)
            seq_sRNA=self.getSequences_53(sRNA, tab_sRNA)
            #Add comparison
            for i in range(len(seq_mRNA)):
                #debut mRNA fin mRNA debut sRNA fin sRNA score
                #root.addComparison([tab_mRNA[i][0],tab_mRNA[i][1],tab_sRNA[i][1],tab_sRNA[i][0],root.scoreSequence([seq_mRNA[i], seq_sRNA[i]])])            
                root.addComparison([tab_mRNA[i][0],tab_mRNA[i][1],tab_sRNA[i][0],tab_sRNA[i][1],root.scoreSequence([seq_mRNA[i], seq_sRNA[i]])])
    
    def parseYass(self, data, matrix, root, complist_file):
        """
        Parse Yass output
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        @param complist_file: Comparison list file
        """
        #seq="(?!\*\.)"
        #seq+="(\S+)\s+"
        seq="^(NC_[0-9]+_[A-Za-z0-9]+)\S+\s+"
        seq+="(\S+)\s+"*10
        seq+="(\S+)\s*"
        #"(?!\*)^(\S+)\s+"*11+"(\S+)\s*"
        regex=re.compile(seq)
        #Lecture du fichier
        for i in data.split("\n"):
            a=regex.match(i)
            if a:
                #debut mRNA fin mRNA debut sRNA fin sRNA score
                root.addComparison([int(a.group(10)),int(a.group(9)),int(a.group(7)),int(a.group(8)),float(a.group(11))])

    
    def parseSsearch(self, data, sRNA, matrix, root, complist_file):
        """
        Parse ssearch output
        @param data: Soft result
        @param sRNA: sRNA name
        @param matrix: score matrix
        @param root: xml root node
        @param complist_file: Comparison list file
        """
        
        #seq="^(NC_[0-9]+_[A-Za-z0-9]+).*\[r\]\s+"
        seq="^(\S+).*\[r\]\s+"
        seq+="(\S+)\s+"*17
        #seq+="(\S+)\s*"
        previous_sRNA=[]
        regex=re.compile(seq)
        #Presence d'une liste de comparaison
        if(complist_file):
            comp_data=self.readCompfile(complist_file)
            #Lecture du fichier
            for i in data.split("\n"):
                a=regex.match(i)
                if a:
                    #mRNA sRNA
                    if(self.comparaison(a.group(1), sRNA, comp_data) and (a.group(1) not in previous_sRNA)):
                        #eviter les doublons
                        previous_sRNA+=[a.group(1)]
                        #mRNA sRNA
                        root.addListOfComparison(a.group(1), sRNA)
                        #debut mRNA fin mRNA debut sRNA fin sRNA score
                        root.addComparison([int(a.group(13)),int(a.group(14)),int(a.group(9)),int(a.group(10)),float(a.group(4))])
        else:
            for i in data.split("\n"):
                a=regex.match(i)
                if (a  and (a.group(1) not in previous_sRNA)):
                    #eviter les doublons
                    previous_sRNA+=[a.group(1)]
                    #mRNA sRNA
                    root.addListOfComparison(a.group(1), sRNA)
                    #debut mRNA fin mRNA debut sRNA fin sRNA score
                    root.addComparison([int(a.group(13)),int(a.group(14)),int(a.group(9)),int(a.group(10)),float(a.group(4))])
    
    def parsePairfold(self, data, matrix, root):
        """
        Parse Pairfold result
        @param data: Soft result
        @param matrix: score matrix
        @param root: xml root node
        """
        sRNA_r=""
        mRNA_r=""
        energy=None
        #Definition des regex
        regex_result=re.compile("^MFE:\s+([\.\(\)]+)\s+([\.\(\)]+)\s+(\S+)")        
        #Lecture du fichier
        for i in data.split("\n"):
            a=regex_result.match(i)
            if a:
                sRNA_r=a.group(1)
                mRNA_r=a.group(2)
                energy=a.group(3)
        #Analyse des positions
        tab_sRNA=self.getposition_RNAcofold0(sRNA_r)
        tab_mRNA=self.getposition_RNAcofold1(mRNA_r)
        #Si des positions sont detectees
        if(tab_sRNA and tab_mRNA):
            #On reverse le resultat de mRNA pour le traiter facilement
            tab_mRNA.reverse()
            #Elongation des positions 
            tab_sRNA=self.elongatePositions_rip(tab_sRNA)
            tab_mRNA=self.elongatePositions_rip(tab_mRNA)
            #Reversion pour l'elongate
            tab_mRNA=self.reverseElongate(tab_mRNA)
            #Reverse the table
            tab_mRNA,maxValue=self.forwardTable(tab_mRNA)
            #Find common position
            (tab_sRNA,tab_mRNA)=self.consensus([],[],tab_sRNA, tab_mRNA, tab_sRNA[0][0], tab_mRNA[0][0], 0, 0)
            #Reverse again the table
            tab_mRNA=self.reverseTablemax(tab_mRNA,maxValue)
            #Add comparison
            for i in range(len(tab_mRNA)):
                #debut mRNA fin mRNA debut sRNA fin sRNA score
                root.addComparison([tab_mRNA[i][1],tab_mRNA[i][0],tab_sRNA[i][0],tab_sRNA[i][1],energy])
    
    def runParsing(self, soft_name, matrix, complist_file, sRNA, mRNA, data, root):
        """
        Run the parsing
        @param cmd: Soft command
        @param complist_file: Comparison list file
        @param data: Soft result
        """
        #Set matrix
        if(not root.matrix): 
            root.setMatrix(matrix)        
        #Parse result
        if(soft_name=="guugle"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseGuugle(data, matrix, root)
        elif(soft_name=="IntaRNA"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseIntaRNA(data, matrix, root)
        elif(soft_name=="RNAup"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseRNAup(data, matrix, root)
        elif(soft_name=="RNAplex"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseRNAplex(data, matrix, root)
        elif(soft_name=="RNAduplex"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseRNAduplex(data, matrix, root)
        elif(soft_name=="RNAcofold"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseRNAcofold(data, matrix, root)
        elif(soft_name=="RNAhybrid"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseRNAhybrid(data, matrix, root) 
        elif(soft_name=="ractip"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseRactip(data, matrix, root)
        elif(soft_name=="bistarna"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseBistaRNA(data, matrix, root)
        elif(soft_name=="blastall"):
            self.parseBlastall(data, matrix, root, complist_file)
        elif(soft_name=="yass-Linux64.bin"):
            #Set setlement for comparison
            root.addListOfComparison(mRNA, sRNA)
            self.parseYass(data, matrix, root, complist_file)
        elif(soft_name=="ssearch35_t"):
            self.parseSsearch(data, sRNA, matrix, root, complist_file)
        elif(soft_name=="pairfold"):
            root.addListOfComparison(mRNA, sRNA)
            self.parsePairfold(data, matrix, root)
        else:
            sys.exit("Soft name is not recognize = %s"%soft_name)