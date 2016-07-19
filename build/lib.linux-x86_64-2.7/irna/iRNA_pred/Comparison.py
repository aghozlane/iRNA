"""
 @brief: Define the comparison
 @author: Amine Ghozlane
 @version: 1.0 
"""
#from mpi4py import MPI
import re, os, string, subprocess, sys, time, shutil


class Comparaison:
    """
    @brief: Create coupled files and the list of comparison
    """
    def __init__(self, comp_list, sRNA_out, mRNA_out, soft_path, mf_mRNA):
        """
        The constructor
        """
        self.comp_list=comp_list  
        self.sRNA_out=sRNA_out 
        self.mRNA_out=mRNA_out
        self.soft_path=soft_path
        self.mf_mRNA=mf_mRNA
        try:
            #Indique le path du soft
            if(self.soft_path[-1]=="/"):
            #gestion des fins de path
                self.soft_path=self.soft_path[:-1]
        except:
            pass
    
    def runCommand(self, cmd):
        """
        Run command
        @param cmd: Command to run
        """
        try:
            #Execution de la ligne de commande
            retcode = subprocess.call(cmd,shell=True)
            #Cas aucun retour du soft
            if retcode == None:
                sys.exit("Child was terminated")
        except OSError, e:
            sys.exit("Execution failed: %s"%(e))
        except:
            sys.exit("There is something wrong with the command: %s"%(cmd))
    
    def createCouple(self, sRNA, mRNA):
        """
        Create the couple
        @param sRNA: name of the sRNA file
        @param mRNA: name of the mRNA file
        """
        #Definition du fichier
        #out=self.sRNA_out+"couple_"+sRNA+"_"+mRNA+".fasta"
        out_name='%scouple_%s_%s.fasta'% (self.sRNA_out, sRNA, mRNA)
        try:
            out = open(out_name, 'w')
            src1 = open(self.sRNA_out+sRNA+".fasta", 'r')
            out.write(src1.read())
            src1.close()
            src2 = open(self.mRNA_out+mRNA+".fasta", 'r')
            out.write(src2.read())
            src2.close()
        except IOError:
            sys.exit("Execution failed with %s"%out_name)
        except:
            sys.exit("There is something wrong with the copy of %s"%out_name)
        
 
    
    def createRNAcofoldCouple(self, sRNA, mRNA):
        """
        Create the couple for RNAcofold
        @param sRNA: name of the sRNA file
        @param mRNA: name of the mRNA file
        """
        #Definition du fichier
        out=self.sRNA_out+"couple_RNAcofold_"+sRNA+"_"+mRNA+".fasta"
        #Ecriture du fichier
        try:
            cofold_couple=open(out,"w")
            src1 = open(self.sRNA_out+sRNA+".fasta", 'r')
            src2 = open(self.mRNA_out+mRNA+".fasta", 'r')
            cofold_couple.write(src1.next())
            cofold_couple.write(src2.next())
            cofold_couple.write("%s&"%(src1.next().translate(None,'\n')))
            cofold_couple.write(src2.next())
            src1.close()
            src2.close()
            cofold_couple.close()
        except IOError:
            sys.exit("Error : Can not open file %s"%out)
        except:
            sys.exit("Something went wrong with %s"%out)   
    
    def getCouple(self, sRNA, mRNA):
        """
        Get a couple in sRNA repertory
        @param sRNA: name of the sRNA file
        @param mRNA: name of the mRNA file
        @return: Return the file couple 
        """
        #Definition du fichier
        out=self.sRNA_out+"couple_"+sRNA+"_"+mRNA+".fasta"
        #Test si le fichier existe
        if(not os.path.isfile(out)):
            sys.exit("File is not found %s"%out)
        #couple file
        return out
    
    def getRNAcofoldCouple(self, sRNA, mRNA):
        """
        Get a couple in sRNA repertory
        @param sRNA: name of the sRNA file
        @param mRNA: name of the mRNA file
        @return: Return the file couple 
        """
        #Definition du fichier
        out=self.sRNA_out+"couple_RNAcofold_"+sRNA+"_"+mRNA+".fasta"
        #Test si le fichier existe
        if(not os.path.isfile(out)):
            sys.exit("File is not found %s"%out)
        #couple file
        return out
    
    def getName(self, regex, fasta_file):
        """
        Get the name of the fasta file
        @param regex: Regex used to detect the name 
        @param fasta_file: A standard fasta file with her path
        @return: The name of the fasta file
        """
        #Recupere le nom du fichier fasta
        a=regex.match(fasta_file)
        if not a:
            sys.exit("Something went wrong with the fasta file: %s"%fasta_file)
        #Retourne le nom du fichier
        return a.group(1)
    
    def getCommand(self, build_command, sRNA, mRNA,RNAcofold_flag):
        """
        Build the command 
        @param build_command: Soft command
        @param sRNA: name of the sRNA file
        @param mRNA: name of the mRNA file
        @param RNAcofold_flag: Flag for RNAcofold special case
        @return: Comparison command
        """
        regex_couple=re.compile(".*couple.*")
        #Indique le sRNA
        build_command = string.replace(build_command,'%sRNA',self.sRNA_out+sRNA+".fasta")
        #Indique le mRNA
        build_command = string.replace(build_command,'%mRNA',self.mRNA_out+mRNA+".fasta")
        #Replacement path du soft
        build_command=string.replace(build_command,'%soft',self.soft_path)
        #Indique le multifasta mRNA
        build_command=string.replace(build_command,'%mf_mRNA',self.mf_mRNA)
        #Indique les couples : cas RNAcofold 
        if(regex_couple.match(build_command) and RNAcofold_flag):
            build_command=string.replace(build_command,'%couples-2',self.getRNAcofoldCouple(sRNA,mRNA))
        #Indique les couples cas normal
        elif(regex_couple.match(build_command) and not RNAcofold_flag):
            build_command=string.replace(build_command,'%couples',self.getCouple(sRNA,mRNA))
        #Retourne la commande completee
        return build_command.translate(None,'\n')
    
    def getFastaFiles(self, repertory):
        """
        Get the list of data
        @param repertory: Interest repertory
        @return: Return the list of file of interest
        """
        #return subprocess.Popen("ls %s*.fasta|xargs"%repertory, shell=True, stdout=subprocess.PIPE).stdout.readlines()
        return os.listdir("%s"%repertory)
    
    def getUnique(self,data_list):
        """
        Get unique data
        @param data_list: list of data
        """
        # Dave Kirby
        # Order preserving
        seen = set()
        return [x for x in data_list if x not in seen and not seen.add(x)]
        #return {}.fromkeys(data_list).keys()
    
    def getComparison(self, command, RNAcofold_flag):
        """
        Define the comparison
        @param command: Soft command
        @param RNAcofold_flag: Flag for RNAcofold special case
        @return: Comparison list and the number of comparison
        """
        result_list=[]
        name_list=[]
        ncomparison=0
        blast_mult=0
        #Liste de comparaison
        if(self.comp_list):
            regex=re.compile("(\S+\.fasta)\s+(\S+\.fasta)")
            #Case mf fasta
            try:
                #Ouverture du fichier de comparaison
                comp=open(self.comp_list,"r")
                #Constitution de la liste des commandes a partir de liste
                for i in comp:
                    a = regex.match(i)
                    if a:
                        #Test si les fichiers existent
                        if(os.path.isfile(self.sRNA_out+a.group(1)) and os.path.isfile(self.mRNA_out+a.group(2))):
                            #Transforme la commande
                            result_list.append(self.getCommand(command, a.group(1)[:-6], a.group(2)[:-6],  RNAcofold_flag))
                            name_list.append([a.group(1)[:-6], a.group(2)[:-6]])
                            blast_mult+=1
                        else:
                            sys.stderr.write("One file is missing : %s  or %s"%((self.sRNA_out+a.group(1)), (self.mRNA_out+a.group(2))))
                #Creation de commande unique
                result_list=self.getUnique(result_list)
                ncomparison=len(result_list)
                #cas mf_mRNA
                if(ncomparison!=blast_mult):
                    list_temp=[]
                    #Get data
                    for x in name_list:
                        list_temp.append(x[0])
                    #Recuperation des id uniques
                    list_temp=self.getUnique(list_temp)
                    #reconstitution de la liste
                    name_list=[]
                    for x in list_temp: 
                        name_list.append([x, ""])
                #Fermeture du fichier
                comp.close()
            except IOError:
                sys.exit("Comparison list file is not present : %s"%(self.comp_list))
        #Sans liste de comparaison
        else:
            regex_mf=re.compile(".*\%mf_mRNA*")
            regex_couple=re.compile(".*couple.+.fasta")
            regex_fasta=re.compile(".*.fasta")
            #Listes des fichiers a traiter
            sRNA_files = self.getFastaFiles(self.sRNA_out) 
            mRNA_files = self.getFastaFiles(self.mRNA_out)        
            #Constitution de la liste des commandes
            #Cas d'une comparaison contre un multifasta 
            if(regex_mf.match(command)):
                for i in sRNA_files:
                    #Enleve les fichiers couples
                    if(not regex_couple.match(i) and regex_fasta.match(i)):
                        sRNA=os.path.basename(i)[:-6]
                        result_list.append(self.getCommand(command, sRNA,"", RNAcofold_flag))
                        name_list.append([sRNA, ""])
                        ncomparison+=1
            #Cas d'une comparaison deux a deux
            else:
                for i in sRNA_files:
                    for j in mRNA_files:
                        #Enleve les fichiers couples
                        if(not regex_couple.match(i) and regex_fasta.match(i) and regex_fasta.match(j)):
                            sRNA=os.path.basename(i)[:-6]
                            mRNA=os.path.basename(j)[:-6]
                            result_list.append(self.getCommand(command, sRNA, mRNA, RNAcofold_flag))
                            name_list.append([sRNA, mRNA])
                            ncomparison+=1            

        #return result_list, ncomparison, name_list, comp_mult, blast_mult
        return result_list, ncomparison, name_list
    
    def runComparison(self, cmd):
        """
        Run a comparison
        @param cmd: Soft command
        @return: Return soft result
        """
        output, errors=None,None
        try:
            pipe =  subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT)
            if pipe == None:
                sys.exit("Child was terminated")
            output, errors = pipe.communicate(input=input)
        #Erreur d'execution
        except OSError, e:
            sys.exit("Execution failed: %s"%(e))
        #Erreur inattendue
        except:
            sys.exit("There is something wrong with the command: %s"%(cmd))
        #assert not errors
        return output
    
    def manageCouple(self):
        """
        Manage the creation of couples
        """
        tabcomp=[]
        #Cas d'une liste de comparaison
        if(self.comp_list):
            regex=re.compile("(\S+\.fasta)\s+(\S+\.fasta)")
            try:
                #Ouverture du fichier de comparaison
                comp=open(self.comp_list,"r")
                #Constitution de la liste des commandes a partir de liste
                for i in comp:
                    a = regex.match(i)
                    if a:
                        #Test si les fichiers existent
                        if(os.path.isfile(self.sRNA_out+a.group(1)) and os.path.isfile(self.mRNA_out+a.group(2))):
                            tabcomp+=[[a.group(1)[:-6], a.group(2)[:-6]]]
                        else:
                            sys.stderr.write("One file is missing : %s  or %s"%((self.sRNA_out+a.group(1)), (self.mRNA_out+a.group(2))))                    
            except IOError:
                sys.exit("Comparison list file is not present : %s"%(self.comp_list))
        else:
            #Listes des fichiers a traiter
            sRNA_files = self.getFastaFiles(self.sRNA_out) 
            mRNA_files = self.getFastaFiles(self.mRNA_out)
            regex_couple=re.compile(".*couple.+.fasta")
            regex=re.compile(".*.fasta")
            for i in sRNA_files:
                for j in mRNA_files:
                    #Enleve les fichiers couples
                    if(not regex_couple.match(i) and regex.match(i) and regex.match(j)):
                        tabcomp+=[[os.path.basename(i)[:-6], os.path.basename(j)[:-6]]]
        return tabcomp
  
    def prepareComparaison(self):
        """
        Control Mpi procedure
        """
        #Cas du Maitre
        j=0
        tab=self.manageCouple()
        max_files=len(tab)
        print("Create coupled files :")
        while(j<max_files):
            if(j%1000==0):
                print("(%d/%d)"%(j+1,max_files))
            #Creation du couple
            self.createCouple(tab[j][0], tab[j][1])
            #Creation du couple pour RNAcofold
            self.createRNAcofoldCouple(tab[j][0], tab[j][1])           
            j+=1
        print("(%d/%d)"%(j,max_files))