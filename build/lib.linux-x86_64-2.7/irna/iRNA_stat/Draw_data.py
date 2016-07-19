"""
 @brief: Draw and print statistic data
 @author: Amine Ghozlane
 @version: 1.0 
"""
import re, csv, sys
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import numpy as np
from Analysis import *
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
graphics = importr('graphics')
grdevices = importr('grDevices')
base = importr('base')
stats = importr('stats')
pROC = importr('pROC')
Cairo = importr('Cairo')

class draw_data(Analysis):
    
    def __init__(self,result):
        """
        Instanciate draw_data object 
        """
        Analysis.__init__(self)
        self.result=result
    
    def write_curve_param_data(self,listobj,dbmanage):
        """
        @param listobj: List of Computer object
        @param dbmanage: Access to the database
        """
        output=self.result+'soft_param.txt'
        listsRNA=dbmanage.getallRNA(0)
        try:
            with open(output, 'wt') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(["soft","sRNA","Slope","Intercept"])
                for obj in listobj:
                    if(len(obj.curve_param) !=0):                    
                        for i in xrange(len(obj.unique_sRNAidinint)):
                            if(obj.curve_param[i][0]!=None):
                                cuve_param_nan=np.isnan(obj.curve_param[i])
                                if(not True in cuve_param_nan):
                                    writer.writerow([obj.name+"_"+str(obj.softid),listsRNA[i][0],obj.curve_param[i][0],obj.curve_param[i][1]])
        except IOError: sys.exit("Error : the program can not open %s"%output)
            
#    def plot_linear_regression(self,obj, soft_inf, numsoft,dbmanage):
#        """
#        Plot linear regression of random data
#        """
#        listsRNA=dbmanage.getallRNA(0)
#        i=0
#        #Vecteur de parametre        
#        for srnaid in obj.unique_sRNAidinint:
#            #Score of one sRNA
#            numpair=np.where(obj.sRNAid_tab==srnaid)
#            cumnegtemp=obj.cumneg[i]
#            if cumnegtemp!=None:
#                normscoretemp=obj.norm_score[numpair]
#                #Nan values
#                cumneg_nanvalues=np.isnan(cumnegtemp)
#                normscore_nanvalues=np.isnan(normscoretemp)
#                #Filter commun nan values
#                commun=self.commun_values(cumneg_nanvalues,normscore_nanvalues)
#                if(False in commun):
#                    #only write values
#                    cumnegtemp=cumnegtemp[np.where(commun==False)]
#                    normscoretemp=normscoretemp[np.where(commun==False)]
#                    #If curve param available
#                    cuve_param_nan=np.isnan(obj.curve_param[i])
#                    if(not True in cuve_param_nan):
#                        #Print linear regression line
#                        plt.plot(normscoretemp,obj.curve_param[i][0]*normscoretemp+obj.curve_param[i][1])
#                        legend='y= %.2f x + %.2f'%(obj.curve_param[i][0],obj.curve_param[i][1])
#                        #plot legend
#                        plt.legend((legend,),'upper right')
#                    #Print data
#                    plt.plot(normscoretemp,cumnegtemp,".")              
#                    #Legend the graph
#                    
#                    xtext='normalised score'
#                    val=soft_inf.score_type[numsoft]
#                    if(val==2 or val==3):
#                        xtext="negative "+xtext
#                    plt.xlabel(xtext)
#                    plt.ylabel('log(-log(F(x)))')
#                    plt.title('Fitting extreme value distribution')
#                    #Save data
#                    plt.savefig(self.result+obj.name+"_"+str(obj.softid)+"_"+listsRNA[i][0]+"_norm.png")
#                    #Clear plot
#                    plt.clf()     
#            i+=1    

    def plot_linear_regression(self,obj, soft_inf, numsoft,dbmanage):
        """
        Plot linear regression of random data
        @param obj: Computer object
        @param soft_inf: Software information
        @param numsoft: Number of the soft
        @param dbmanage: Access to the database
        """
        listsRNA=dbmanage.getallRNA(0)
        i=0
        #Legend the graph
        val=soft_inf.score_type[numsoft]
        xtext='normalised score'
        if(val==2 or val==3):
            xtext="negative "+xtext
        #Vecteur de parametre        
        for srnaid in obj.unique_sRNAidinint:
            #Score of one sRNA
            numpair=np.where(obj.sRNAid_tab==srnaid)
            cumnegtemp=obj.cumneg[i]
            if cumnegtemp!=None:
                normscoretemp=obj.norm_score[numpair]
                #Nan values
                cumneg_nanvalues=np.isnan(cumnegtemp)
                normscore_nanvalues=np.isnan(normscoretemp)
                #Filter commun nan values
                commun=self.commun_values(cumneg_nanvalues,normscore_nanvalues)
                if(False in commun):
                    #only write values
                    cumnegtemp=cumnegtemp[np.where(commun==False)]
                    normscoretemp=normscoretemp[np.where(commun==False)]
                    #If curve param available
                    cuve_param_nan=np.isnan(obj.curve_param[i])
                    Cairo.CairoPNG(**{"height":1024,"width":1024,"filename":self.result+obj.name+"_"+str(obj.softid)+"_"+listsRNA[i][0]+"_norm.png"})
                    if(not True in cuve_param_nan):
                        kwargs={"xlab":xtext,"ylab":"log(-log(F(x)))","pch":18}
                        graphics.plot(robjects.FloatVector(normscoretemp),robjects.FloatVector(obj.curve_param[i][0]*normscoretemp+obj.curve_param[i][1]),**kwargs)
                        kwargs={"legend":'y= %.2f x + %.2f'%(obj.curve_param[i][0],obj.curve_param[i][1]),"fill":2,"border":"white"}
                        graphics.legend("topright",**kwargs)
                    else:
                        graphics.plot(robjects.FloatVector(normscoretemp),robjects.FloatVector(cumnegtemp))
                    grdevices.dev_off()   
            i+=1    

    def getnames(self,listobj):
        """
        name of sofwares
        @param listobj: List of Computer object
        """
        names=[]
        for obj in listobj:
            names+=[obj.name+"_"+str(obj.softid)]       
        return names
    
    def plot_roc_curves(self,listobj):
        """
        Plot roc curve
        @param listobj: List of Computer object
        """
        add=False
        i=0
        auc_result=np.arange(float(len(listobj)))
        palette = grdevices.rainbow(len(listobj))
        Cairo.CairoPNG(**{"height":1024,"width":1024,"filename":self.result+"roc_curve.png"})
        v = robjects.FloatVector([0,0,1,1,0.4,0,1,0.5])
        m = robjects.r['matrix'](v, nrow = 2)
        graphics.split_screen(fig=m)
        graphics.screen(2)
        kwargs = {"ty":"n","xaxt":"n","yaxt":"n","xlab":"","ylab":"","col":"white","col.lab":"white","fg":"white"}
        graphics.plot(robjects.IntVector([1,5]),robjects.IntVector([1,5]),**kwargs)
        kwargs = {"legend":robjects.StrVector(self.getnames(listobj)),"fill":robjects.StrVector(palette),"border":"white","ncol":3,"cex":1.5}
        graphics.legend(1.75,3,**kwargs)
        graphics.screen(1)
        regex=re.compile(".*([0-9]+\.[0-9]+).*")
        coords_roc=[]
        #Compute roc curve for each software
        for obj in listobj:
            #check values
            commun=np.isnan(obj.pValue)
            if(obj.pValue !=None and  False in commun):
                #"cex.lab":2.0,"cex.axis":2.0,"cex":2,
                kwargs ={"add":add, "plot":"TRUE", "na.rm":"TRUE","col":palette[i],"print.thres":"best","print.thres.best.method":"youden"}
                roc_data=pROC.roc(robjects.IntVector(obj.classtype),robjects.FloatVector(obj.pValue),**kwargs)
                #Get auc data
                b=regex.match(str(roc_data[12]))
                if(b): auc_result[i]=float(b.group(1))
                else: auc_result[i]= np.nan
                #Method de youden
                coords=pROC.coords(roc_data, "b", ret=robjects.StrVector(["se","sp","t"]))
                coords_roc+=[[coords[0],coords[1],coords[2]]]
                add=True
            else: coords_roc+=[None]
            i+=1
        graphics.close_screen()
        grdevices.dev_off()
        return auc_result,coords_roc
    
    def plot_roc_curves_statistics(self,listobj):
        """
        Statistical analysis of roc
        @param listobj: List of Computer object
        """
        i=0
        print("roc statistics:")
        #Compute statistics of roc curve for each software
        for obj in listobj:
            commun=np.isnan(obj.pValue)
            if(obj.pValue !=None and  False in commun):
                Cairo.CairoPNG(**{"height":1024,"width":1024,"filename":self.result+obj.name+"_"+str(obj.softid)+"_roc_curve_statistics.png"})
                kwargs ={"ci":True,"print.auc":True,"plot":True,"percent":True, "na.rm":True,"col":"black","print.thres":"best","print.thres.best.method":"youden"}
                roc_data=pROC.roc(robjects.IntVector(obj.classtype),robjects.FloatVector(obj.pValue),**kwargs)
                ciobj=pROC.ci_se(roc_data,**{"specificities":robjects.IntVector(range(0,101,5))}) 
                graphics.plot(ciobj, **{"type":"shape", "col":"#1c61b6AA"}) # plot as a blue shape    
                graphics.plot(pROC.ci(roc_data, **{"of":"thresholds", "thresholds":"best"})) # add one threshold
                grdevices.dev_off()
                i+=1
                print("(%d/%d)"%(i,len(listobj)))        
    
    def write_pValue(self,listobj,dbmanage):
        """
        Print pValue
        @param listobj: List of Computer object
        @param dbmanage: Access to the database
        """
        output=self.result+'pValue.txt'
        try:
            with open(output, 'wt') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(["soft","sRNA","mRNA","pValue","class"])
                for obj in listobj:
                    i=0
                    for interactid in obj.interactions:
                        writer.writerow([obj.name+"_"+str(obj.softid),dbmanage.getsRNAbyIntid(interactid[0]),dbmanage.getmRNAbyIntid(interactid[0]),obj.pValue[i],obj.classtype[i]])
                        i+=1
        except IOError: sys.exit("Error : the program can not open %s"%output)
    
    def write_rocthresValue(self,listobj,coords_roc,dbmanage):
        """
        @param listobj: List of Computer object
        @param coords_roc: Roc coordinates
        @param dbmanage: Access to the database
        """
        i=0
        output=self.result+'pValue_thres.txt'
        try:
            with open(output, 'wt') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(["soft","sensitivity","specificity","threshold"])
                for obj in listobj:
                    if(coords_roc[i]!=None):
                        writer.writerow([obj.name+"_"+str(obj.softid),coords_roc[i][0],coords_roc[i][1],coords_roc[i][2]])       
                    i+=1
        except IOError: sys.exit("Error : the program can not open %s"%output)
    
    def write_rocAuc(self,listobj,auc_result):
        """
        write AUC
        @param listobj: List of Computer object
        @param auc_result: AUC result
        """
        i=0
        output=self.result+'pValue_auc.txt'
        try:
            with open(output, 'wt') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(["soft","AUC"])
                for obj in listobj:
                    writer.writerow([obj.name+"_"+str(obj.softid),auc_result[i]])       
                    i+=1
        except IOError: sys.exit("Error : the program can not open %s"%output)
    
    def meanArray(self,array):
        """
        Mean array
        @param array: Matrix array
        """
        nanarray=np.where(np.isnan(array)==False)[0]
        if(len(nanarray)>1): return(np.ma.mean(np.ma.masked_array(array,np.isnan(array))))
        elif(len(nanarray)==1): return(array[nanarray])
        else: return None
    
    def getPredplotvect(self,listobj):
        """
        Get float vector from listobj
        @param listobj: List of Computer object
        """
        result_vect=[]
        for obj in listobj:
            sens=self.meanArray(obj.sensitivity)
            ppv=self.meanArray(obj.ppv)
            if sens==None: sens=0.0
            if ppv==None: ppv=0.0
            result_vect+=[sens, ppv]
        return robjects.FloatVector(result_vect)
    
    def pred_plot(self,listobj):
        """
        Plot sensitivity and ppv mean
        @param listobj: List of Computer object
        """
        data_matrix = robjects.r['matrix'](self.getPredplotvect(listobj), nrow = 2)
        Cairo.CairoPNG(**{"height":600,"width":2000,"filename":self.result+'pred_plot.png'})
        #kwargs={"beside":"TRUE"}
        kwargs={"beside":True,"names.arg":self.getnames(listobj),"col":robjects.StrVector(["black","red"]),"ylim":robjects.IntVector([0,1]),"ylab":"Value (%)","cex.lab":2.0,"cex.axis":2.0,"cex":1}
        graphics.barplot(data_matrix,**kwargs)
        #graphics.barplot(data_matrix)        
        kwargs={"legend":robjects.StrVector(["Sensitivity","PPV"]),"fill":robjects.StrVector(["black","red"]),"border":"white","cex":2}
        graphics.legend("topright",**kwargs)
        grdevices.dev_off()
    
    
    def write_pred_plot(self,listobj):
        """
        Write mean sensitivity and ppv of each software
        @param listobj: List of Computer object
        @param dbmanage: Access to the database
        """
        output=self.result+'pred_plot.txt'
        try:
            with open(output, 'wt') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(["soft","Mean sensitivity","Mean PPV"])
                for obj in listobj:
                    writer.writerow([obj.name+"_"+str(obj.softid),self.meanArray(obj.sensitivity), self.meanArray(obj.ppv)])
        except IOError: sys.exit("Error : the program can not open %s"%output)           
    
    def write_sens_ppv(self,listobj,dbmanage):
        """
        Write sensitivity and ppv
        @param listobj: List of Computer object
        @param dbmanage: Access to the database
        """
        i=0
        output=self.result+'sens_ppv.txt'
        try:
            with open(output, 'wt') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(["soft","sRNA","mRNA","Sensitivity","PPV"])
                for obj in listobj:
                    i=0
                    for interactid in obj.interactions:
                        writer.writerow([obj.name+"_"+str(obj.softid),dbmanage.getsRNAbyIntid(interactid[0]),dbmanage.getmRNAbyIntid(interactid[0]),obj.sensitivity[i],obj.ppv[i]])
                        i+=1
        except IOError: sys.exit("Error : the program can not open %s"%output) 
    
    def exec_plot(self,execution_time):
        """
        Plot execution time
        @param execution_time: Execution time
        """
        Cairo.CairoPNG(**{"height":600,"width":2000,"filename":self.result+'exec_plot.png'})
        kwargs={"col":"black","names.arg":robjects.StrVector(execution_time.soft),"ylab":"Time (s)","cex.lab":2.0,"cex.axis":2.0,"cex":1}
        graphics.barplot(robjects.FloatVector(execution_time.duration),**kwargs)
        grdevices.dev_off()
    
    def write_pValueSelect(self,listobj):
        """
        Write selected interaction
        @param listobj: List of Computer object
        """
        for obj in listobj:
            output=self.result+'%s_%d_multilist.txt'%(obj.name,obj.softid)
            try:
                with open(output, 'wb') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(obj.allsRNA)
                    writer.writerows(obj.matrix)
            except IOError: sys.exit("Error : the program can not open %s"%output) 
    
    def write_frequency(self,listobj):
        """
        Write mRNA frequency data
        @param listobj: List of Computer object
        """
        for obj in listobj:
            output=self.result+'%s_%d_frequency_mRNA.txt'%(obj.name,obj.softid)
            try:
                with open(output, 'wb') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(["mRNA","Frequency (%)"])
                    for i in obj.frequency:
                        writer.writerow([i["mRNA"],i["frequency"]])
            except IOError:
                sys.exit("Error : the program can not open %s"%output)
    
    def write_similarity(self,listobj):
        """
        Write group similarity data
        @param listobj: List of Computer object
        """
        for obj in listobj:
            output=self.result+'%s_%d_similarity_mRNA.txt'%(obj.name,obj.softid)
            try:
                with open(output, 'wb') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(["sRNA_1","sRNA_2","Similarity indice"])
                    for i in obj.similarity:
                        writer.writerow([i["sRNA_1"],i["sRNA_2"],i["similarity"]])
            except IOError: sys.exit("Error : the program can not open %s"%output) 
    
    def frequency_plot(self,listobj):
        """
        Plot frequency data
        @param listobj: List of Computer object
        """
        for obj in listobj:
            Cairo.CairoPNG(**{"height":1024,"width":1024,"filename":self.result+'%s_%d_frequency_mRNA.png'%(obj.name,obj.softid)})
            kwargs={"beside":True,"names.arg":robjects.StrVector([i["mRNA"] for i in obj.frequency]),"xlab":"mRNA","ylab":"Frequency"}
            graphics.barplot(robjects.FloatVector([i["frequency"] for i in obj.frequency]),**kwargs)
            grdevices.dev_off()

    def similarity_plot(self,listobj):
        """
        @param listobj: List of Computer object
        """
        for obj in listobj:
            Cairo.CairoPNG(**{"height":1024,"width":1024,"filename":self.result+'%s_%d_similarity_mRNA.png'%(obj.name,obj.softid)})
            kwargs={"beside":True,"names.arg":robjects.StrVector([i["sRNA_1"]+"_"+i["sRNA_2"] for i in obj.similarity]),"xlab":"mRNA","ylab":"Frequency"}
            graphics.barplot(robjects.FloatVector([i["similarity"] for i in obj.similarity]),**kwargs)
            grdevices.dev_off()         
    
    def write_pValueInteract(self,listobj,dbmanage):
        """
        @param listobj: List of Computer object
        @param dbmanage: Access to the database
        """
        for obj in listobj:
            output=self.result+'%s_%d_pValue.txt'%(obj.name,obj.softid)
            try:
                with open(output, 'wb') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(["sRNA","mRNA","pValue","selected"])
                    for i in xrange(obj.nbinteract):
                        val=0
                        if(obj.select[i]): val=1
                        writer.writerow([dbmanage.getsRNAbyIntid(obj.interactions[i][0]),dbmanage.getmRNAbyIntid(obj.interactions[i][0]),obj.pValue[i],val])
            except IOError: sys.exit("Error : the program can not open %s"%output) 
    
    def getBoolselect(self,selection):
        """
        Indicate if mRNA is selected as target
        @param selection: boolean of selection
        @return: Integer value
        """
        if(selection): return 1
        return 0
    
    def write_interaction(self,listobj,dbmanage,selection):
        """
        Write interaction data
        @param listobj: List of Computer object
        @param dbmanage: Access to the database
        @param selection: Bool for selection printing
        """
        for obj in listobj:
            output=self.result+'%s_%d_interaction.txt'%(obj.name,obj.softid)
            try:
                with open(output, 'wb') as f:
                    writer = csv.writer(f, delimiter='\t')
                    if(selection): writer.writerow(["sRNA","mRNA","sRNA_begin","sRNA_end","mRNA_begin","mRNA_end","sRNA_interaction_length","mRNA_interaction_length","sRNA_length","mRNA_length","selected"])
                    else: writer.writerow(["sRNA","mRNA","sRNA_begin","sRNA_end","mRNA_begin","mRNA_end","sRNA_interaction_length","mRNA_interaction_length","sRNA_length","mRNA_length","class"])
                    for i in xrange(obj.nbinteract):
                        if(obj.interact[i]!=None):
                            if(selection): writer.writerow([dbmanage.getsRNAbyIntid(obj.interactions[i][0]),dbmanage.getmRNAbyIntid(obj.interactions[i][0]),obj.interact[i][0],obj.interact[i][1],obj.interact[i][2],obj.interact[i][3],obj.interact[i][4],obj.interact[i][5],dbmanage.getsRNAlenbyIntid(obj.interactions[i][0]),dbmanage.getmRNAlenbyIntid(obj.interactions[i][0]),self.getBoolselect(obj.select[i])])
                            else: writer.writerow([dbmanage.getsRNAbyIntid(obj.interactions[i][0]),dbmanage.getmRNAbyIntid(obj.interactions[i][0]),obj.interact[i][0],obj.interact[i][1],obj.interact[i][2],obj.interact[i][3],obj.interact[i][4],obj.interact[i][5],dbmanage.getsRNAlenbyIntid(obj.interactions[i][0]),dbmanage.getmRNAlenbyIntid(obj.interactions[i][0]),obj.classtype[i]])
                        else:
                            if(selection): writer.writerow([dbmanage.getsRNAbyIntid(obj.interactions[i][0]),dbmanage.getmRNAbyIntid(obj.interactions[i][0]),np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,dbmanage.getsRNAlenbyIntid(obj.interactions[i][0]),dbmanage.getmRNAlenbyIntid(obj.interactions[i][0]),self.getBoolselect(obj.select[i])])
                            else: writer.writerow([dbmanage.getsRNAbyIntid(obj.interactions[i][0]),dbmanage.getmRNAbyIntid(obj.interactions[i][0]),np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,dbmanage.getsRNAlenbyIntid(obj.interactions[i][0]),dbmanage.getmRNAlenbyIntid(obj.interactions[i][0]),obj.classtype[i]])
            except IOError: sys.exit("Error : the program can not open %s"%output) 
    
    def write_general_data(self,listobj):
        """
        Write general data on prediction
        @param listobj: List of Computer object
        """
        for obj in listobj:
            output=self.result+'%s_%d_general_data.txt'%(obj.name,obj.softid)
            try:
                with open(output, 'wb') as f:
                    interaction=self.getRangePosit(obj.interact,4)
                    if(len(interaction)!=0):
                        f.write("Mean interaction length : %.0f nt\n"%np.mean(interaction))
                        f.write("Maximum interaction length : %.0f nt\n"%np.max(interaction))
                        f.write("Minimum interaction length : %.0f nt\n"%np.min(interaction))
            except IOError: sys.exit("Error : the program can not open %s"%output)     
    
    def plot(self,args,soft_inf,listobj,dbmanage,execution_time):
        """
        Plot main function
        """
        #Plot linear regression
        if(args.random):
            #print curve param
            self.write_curve_param_data(listobj,dbmanage)
            for obj in listobj:
                if(obj.cumneg !=None and len(obj.curve_param) !=0): self.plot_linear_regression(obj,soft_inf,obj.numsoft,dbmanage)
        #Roc curve analyse
        if(not args.random and args.exp_inf):
            #Plot roc curve 
            auc_result,coords_roc=self.plot_roc_curves(listobj)
            self.plot_roc_curves_statistics(listobj)
            #print
            self.write_pValue(listobj,dbmanage)
            self.write_rocthresValue(listobj,coords_roc,dbmanage)
            self.write_rocAuc(listobj,auc_result)
            self.write_interaction(listobj,dbmanage,False)
            self.write_general_data(listobj)
        #pValue_inf
        elif(not args.random and not args.exp_inf):
            self.write_pValueInteract(listobj,dbmanage)
            self.write_interaction(listobj,dbmanage,True)
            self.write_general_data(listobj)
        #Experiment interaction
        if(args.exp_inf):
            self.pred_plot(listobj)
            self.write_pred_plot(listobj)
            self.write_sens_ppv(listobj,dbmanage)
        #Software execution time
        if(execution_time):
            self.exec_plot(execution_time)
        #Print selected mRNA
        if(args.thres_inf):
            self.write_pValueSelect(listobj)
            self.write_frequency(listobj)
            self.write_similarity(listobj)
            self.frequency_plot(listobj)
            self.similarity_plot(listobj)