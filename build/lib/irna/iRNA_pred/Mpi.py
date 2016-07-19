"""
 @brief: Run mpi
 @author: Amine Ghozlane
 @version: 1.0
"""
import sys, os, time
#import zlib
from mpi4py import MPI
from lxml import etree
from copy import deepcopy
from Comparison import *
from Parse import *
from Interaction import *

class Mpi:
    """
    @brief: Manage Mpi run 
    """
    def __init__(self):
        """
        The Constructor
        @note comm: Broadcast communicator
        @note nprocs: Number of thread
        @note myrank: Rank of the thread
        """
        self.comm = MPI.COMM_WORLD
        self.nprocs = self.comm.Get_size()
        self.myrank = self.comm.Get_rank()
        self.name = MPI.Get_processor_name()

    def Mpi_getComm(self):
        """
        @return: Communicator
        """
        return self.comm
    
    def Mpi_getnprocs(self):
        """
        @return: The number of processor
        """
        return self.nprocs
    
    def Mpi_getmyrank(self):
        """
        @return: The rank value
        """
        return self.myrank
    
    def Mpi_getname(self):
        """
        @return: The name
        """
        return self.name

    def Mpi_write_exectime(self, exec_time, result_out):
        """
        Write exectime  of command
        @param exec_time: Tab of execution time for each command
        @param result_out: Result repertory
        """
        out=result_out+"execution_time.txt"
        #Recuperation de la liste des softs
        soft_list=[]
           
        try:
            #Ouverture du fichier de resultat
            exec_file=open(out,"w")
            exec_file.write("Soft_command\tDuration(s)\n")
            #Ecriture des temps d'execution par commande
            for i in exec_time:
                #Calcul du temps d'execution
                #Bug en cas de changement de mois ou d'annee entre les deux
                start=i[1].tm_sec+i[1].tm_min*60+i[1].tm_hour*3600+i[1].tm_mday*24*3600
                end=i[2].tm_sec+i[2].tm_min*60+i[2].tm_hour*3600+i[2].tm_mday*24*3600
                #Ecriture des donnees
                soft_list.append(self.Mpi_getSoft(soft_list, i[0]))
                exec_file.write("%s\t%f\n"%(soft_list[-1],(end-start)))
            #Fermeture du fichier
            exec_file.close()
        except IOError:
            sys.stderr.write("Error : can not open file %s"%out)
        except :
            sys.stderr.write("Something went wrong with %s"%out)
    
    def Mpi_master(self, soft_list, comp_list, result_out, matrix, compair):
        """
        Master process
        @param soft_list: List of soft
        @param comp_list: List of comparison
        @param result_out: Result repertory
        @param matrix: Matrix file
        @param compair: Comparison object
        """
        status = MPI.Status()
        regex=re.compile("^(?!#)\S+\/(\S+)\s+\S+")
        Lsoft=""
        work_failed=False
        ref=0
        #Tableau de mesure du temps d'execution
        exec_time=[]
    
        try:
            #Ouverture du fichier
            Lsoft=open(soft_list,"r")
            #Pour chaque ligne de commande
            for i in Lsoft:
                #Nom du soft
                soft_name=""
                a=regex.match(i)
                RNAcofold_flag=False                
                #Soft a lancer
                if(a):
                    #Recuperation du nom du soft
                    soft_name=a.group(1)
                    #RNAcofold
                    if(soft_name=="RNAcofold"):
                        RNAcofold_flag=True     
                    #Recuperation du tableau des comparaisons
                    #comparison_list, max_comparison, interaction_list, comp_mult, blast_mult=compair.getComparison(i,RNAcofold_flag)
                    comparison_list, max_comparison, interaction_list=compair.getComparison(i,RNAcofold_flag)
                    #Test du nombre comparaison
                    if(max_comparison==0):
                        work_failed=True
                        #Fin de processus
                        self.Mpi_exit_processus()
                        sys.stderr.write("There is no comparison to do\n")
                        break
                    else:
                        #Envoie le nom du soft au processus ecrivain 
                        #self.comm.send(soft_name, dest=1, tag=1)
                        #Envoie du nombre de comparaison au processus ecrivain
#                        if(comp_mult!=0):
#                            self.comm.send((max_comparison*comp_mult), dest=1, tag=10)
#                        elif(blast_mult!=0):
#                            self.comm.send(blast_mult, dest=1, tag=10)
#                        else:
#                            self.comm.send(max_comparison, dest=1, tag=10)
                        #Affichage du nom du soft et du nombre de comparaison
                        print("soft : %s - Comparison : %d"%(soft_name, max_comparison))
                        #Debut du decompte de temps
                        a=time.gmtime()
                        #Envoie du soft etudie et flag liste de comparaison
                        tab=[soft_name, matrix, comp_list,ref]
                        #Lancement des calculs
                        j=0
                        while(j<max_comparison):
                            #Test de reduction de consommation de ressource
                            while not self.comm.Iprobe(source=MPI.ANY_SOURCE, tag=1):
                                time.sleep(0.1)
                            #Processus libre
                            self.comm.recv(source=MPI.ANY_SOURCE, tag=1,status=status)
                            #Envoie d'une tache
                            #self.comm.send(comparison_list[j], dest=status.source, tag=1)
                            #Envoie du soft et du couple
                            self.comm.send(list([comparison_list[j]])+tab+interaction_list[j], dest=status.source, tag=1)                        
                            j+=1
                            print("(%d/%d)"%(j,max_comparison))
                        #Signal d'ecriture des donnees
                        for i in range(1,self.nprocs):
                            self.comm.send([soft_name,ref], dest=i, tag=2)
                        #Test de reduction de consommation de ressource
                        #while not self.comm.Iprobe(source=1, tag=3):
                        #    time.sleep(1.0)
                        #Attendre que le processus 1 soit ok
                        #self.comm.recv(source=1, tag=3)
                        #Fin du decompte de temps
                        b=time.gmtime()
                    exec_time.append((soft_name,a,b))
                    ref+=1
            #Fermeture du fichier
            Lsoft.close()
        except IOError:
            sys.exit("Error: We can not open the soft list file : %s"%(soft_list))        
        #Fin de processus
        if(not work_failed):
            self.Mpi_sup_process(status,j,max_comparison)
            self.Mpi_end_processus()
        
        #Ecriture des temps de calcul
        self.Mpi_write_exectime(exec_time,result_out)
    
    def Mpi_sup_process(self,status,j,max_comparison):
        """
        Termine les processus supplementaire
        @param status: 
        @param j: 
        @param max_files: 
        """
        #terminaison des threads suppelementaire
        if((self.nprocs-1)>max_comparison):
            while(j<(self.nprocs-1)):
                #Test de reduction de consommation de ressource
                while not self.comm.Iprobe(source=MPI.ANY_SOURCE, tag=1):
                    time.sleep(0.1)
                #Processus libre
                self.comm.recv(source=MPI.ANY_SOURCE, tag=1,status=status)
                #decompte reception des autres signals
                j+=1
            
#    def Mpi_flush_tampon(self, soft_result, result_data):
#        """
#        Flush the tampon
#        @param soft_result: Result file
#        @param result_data: List of data from soft
#        """
#        if len(result_data)==1:
#            test=zlib.decompress(result_data[0])
#            soft_result.write(test)
#        else:
#            #Ecriture des donnees dans le fichier
#            for data in result_data:
#                test=zlib.decompress(data)
#                soft_result.write(test)
#
#        #Effacement de la liste
#        del(result_data)
        
    def Mpi_getSoft(self, soft_list, soft_name):
        """
        Test if that soft has already been used
        @param soft_list: List of soft already used
        @param soft_name: Soft name
        @return: soft file name
        """
        soft=soft_name
        num=1
        #Compare par rapport a liste
        for i in soft_list:
            if(soft==i):
                soft=soft_name+"_"+str(num)
                num+=1       
        #Nom du soft modifie
        return soft
    
    def Mpi_write_data(self, result_data, result_out, soft_name, ref, thread, part):
        """
        Write xml result file
        @param result_data: Interact object
        @param result_out: Result repertory
        @param soft_name: Software name
        @param ref: Number of the ref
        @param part: Number of the part
        """
        out=str(result_out)+soft_name+"_"+str(ref)+"_thread_"+str(thread)+"_part_"+str(part)+".xml"
        try:
            soft_result=open(out,"w")
            #soft_result.write(etree.tostring(result_data, method='xml', pretty_print=True))
            soft_result.write(result_data.getResult())
            soft_result.close()
        except IOError: sys.stderr.write("Error : can not open file %s"%out)
        except: sys.stderr.write("Something went wrong with %s"%out)

    #Processus esclave
    def Mpi_slave(self, result_out, buffer_size, compair):
        """
        Slave process
        @param compair: Comparison object
        """
        ar=Parse()
        status = MPI.Status()
        self.comm.send(1,dest=0, tag=1)
        root=Interaction()
        buffer=buffer_size
        flag_soft=True
        part=0
        while(status.tag!=0):
            #Mise en sommeil pour la reduction de consommation de ressource
            while not self.comm.Iprobe(source=0, tag=MPI.ANY_TAG):
                time.sleep(0.1)
            data = self.comm.recv(source=0, tag=MPI.ANY_TAG,status=status)
            if(status.tag==1):
                #Cree le noeud du soft
                if(flag_soft):
                    #compteur et nom du soft
                    root.setSoft(data[1],data[4])
                    flag_soft=False
                #Lancement de la comparaison
                result = compair.runComparison(data[0])
                #Parsing de la sortie
                ar.runParsing(data[1], data[2], data[3], data[5], data[6], result, root)                
                #Renvoyer le resultat de la ligne de commande
                self.comm.send(1,dest=status.source, tag=1)
                #Decrease the buffer
                buffer-=1
                #Vider le buffer
                if(buffer==0):
                    #Envoie des resultats
                    #self.comm.send(zlib.compress(root.getResult()), dest=1, tag=2)
                    self.Mpi_write_data(root, result_out, data[1], data[4], self.myrank, part)
                    #Vider le buffer
                    del(root)
                    #Remise a zero de l'arbre
                    root=Interaction()
                    #Reinitialisation du buffer
                    buffer=buffer_size
                    #Reinitialisation du flag soft
                    flag_soft=True
                    #Incrementation du numero de partie
                    part+=1
                    
            #Changement de soft
            elif(status.tag==2):
                #Envoie des resultats
                #self.comm.send(zlib.compress(root.getResult()), dest=1, tag=2)
                #Write data
                if(buffer!=buffer_size):
                    self.Mpi_write_data(root, result_out, data[0], data[1], self.myrank, part)
                #Reinitialisation du buffer
                buffer=buffer_size
                #Vider le buffer
                del(root)
                #Remise a zero de l'arbre
                root=Interaction()
                #Reinitialisation du flag soft
                flag_soft=True
                #Reinitialisation du numero de partie
                part=0
    
    def Mpi_exit_processus(self):
        """
        End of all process
        """
        #Reception de message
        for i in range(2,self.nprocs):
            #Processus libre
            self.comm.recv(source=MPI.ANY_SOURCE, tag=1)
        #Signal de fin
        for i in range(1,self.nprocs):
            self.comm.send(i,dest=i,tag=0)
    
    def Mpi_end_processus(self):
        """
        End of all process
        """
        for i in range(1, self.nprocs):
            self.comm.send(i,dest=i,tag=0)
        
    def Mpi_finalize(self):
        """
        Fin
        """
        #Synchronisation des threads
        self.comm.Barrier()
        #Finalisation
        MPI.Finalize()
         
    def Mpi_run(self,predict, compair):
        """
        Control Mpi procedure
        @param predict: List of data for comparison
        @param compair: Comparaison object
        """
        #Cas du Maitre
        if self.myrank == 0: self.Mpi_master(predict[0], predict[2],  predict[1], predict[5], compair)
        #Cas de l'esclave
        else: self.Mpi_slave(predict[1], predict[3], compair)
        #Finalisation    
        self.Mpi_finalize()
        