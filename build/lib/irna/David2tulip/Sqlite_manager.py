"""
 @brief: Handle database connection, disconnection and querying.
 @author: Amine Ghozlane
 @version: 1.0 
"""
import sqlite3, sys, time, numpy as np
from types import ListType
try: import sqlitebck
except ImportError: pass

class Sqlite_manager:
    """
    @brief: Handle database connection, disconnection and querying. 
    """
    
    def __init__(self, result,db_file,fastmode):
        """
        Instanciate Sqlite_manager object
        @param result: Path where to write the database
        @param db_file: SQlite db file
        @param fastmode: fastmode flag
        """
        if(result): self.db_path=result+"iRNA_result.db"
        self.db_file=db_file
        self.fastmode=fastmode
        
    def createSQLdb(self):
        """
        Create the database
        """
        #Create RNA database
        with self.conn:
            self.conn.executescript("""
            create table rna (rnaid integer primary key, name varchar unique, length integer, type integer);
            create table soft (softid integer primary key, name varchar);
            create table interact (interactid integer primary key, softid REFERENCES soft(softid), srnaid REFERENCES rna(rnaid), mrnaid REFERENCES rna(rnaid));
            create table contact (contactid integer primary key, interactid REFERENCES interact(interactid), srna_begin integer, srna_end integer, mrna_begin integer, mrna_end integer, score real);
            """)
    
    def connectDB(self):
        """
        Connect to the database
        """
        if(self.db_file): self.conn2=sqlite3.connect(self.db_file,check_same_thread = False)
        else: self.conn2=sqlite3.connect(self.db_path,check_same_thread = False)
        #connect
        with self.conn2: self.conn2.execute("PRAGMA synchronous = OFF")
        if(self.fastmode):
            self.conn=sqlite3.connect(':memory:')
            sqlitebck.copy(self.conn2,self.conn)
        else: self.conn=self.conn2
        self.cur=self.conn.cursor()
        
    def disconnectDB(self):
        """
        Disconnect the database
        """
        #close cursor 
        self.cur.close()
        if(self.fastmode): self.conn2.close()                    
        #Close connection
        self.conn.close()
    
    def disconnectDB2(self):
        """
        Disconnect the database and copy on to the disk
        """
        #close cursor 
        self.cur.close()
        #Copy in memory databse to file
        if(self.fastmode):
            sqlitebck.copy(self.conn,self.conn2)
            self.conn2.close() 
        #Close connection
        self.conn.close()
    
    def setSoft(self,name):
        """
        Insert soft name
        @param name: Soft name
        @return: Last insert id
        """
        try:
            self.cur.execute("insert into soft(name) values(?)",(name,))
            lastrowid=self.cur.lastrowid
            self.conn.commit()
        except sqlite3.OperationalError:
            commit_success=False
            while(not commit_success):
                time.sleep(1.0)
                self.commitRetry()
        return lastrowid
    
    def setRNA(self,RNAtab):
        """
        Insert RNA
        @param RNAtab: Table of RNA
        """
        try:
            if(RNAtab):
                if(type(RNAtab[0]) is ListType):
                    for i in RNAtab:
                        self.cur.execute("insert into rna(name,length,type) values(?,?,?)",(i[0],i[1],i[2]))
                else:
                    self.cur.execute("insert into rna(name,length,type) values(?,?,?)",(RNAtab[0],RNAtab[1],RNAtab[2]))
                self.conn.commit()
        except sqlite3.OperationalError:
            commit_success=False
            while(not commit_success):
                time.sleep(1.0)
                self.commitRetry()
        except sqlite3.IntegrityError:
            text="sRNA"
            if(type(RNAtab[0]) is ListType):
                if(i[2]): text="mRNA"
                sys.exit("Impossible insert %s : %s as it is already present in the database"%(text,i[0]))
            if(RNAtab[2]): text="mRNA"
            sys.exit("Impossible insert %s : %s as it is already present in the database"%(text,RNAtab[0]))

    def setInteract(self,sRNAid,mRNAid,softid):
        """
        Insert interaction
        @param sRNAid: sRNA key
        @param mRNAid: mRNA key
        @param softid: soft key
        @return: Last insert id
        """
        try:
            self.cur.execute("insert into interact(softid,srnaid,mrnaid) values(?,?,?)",(softid,sRNAid,mRNAid))
            lastrowid=self.cur.lastrowid
            self.conn.commit()
        except sqlite3.OperationalError:
            commit_success=False
            while(not commit_success):
                time.sleep(1.0)
                self.commitRetry()
        return lastrowid
    
    def commitRetry(self):
        """
        Retry the commit operation
        @return: State of the commit
        """
        try:
            print("Commit retry not good")
            self.conn.commit()
        #Fail
        except sqlite3.OperationalError: return False
        #success
        return True
    
    def setContact(self,interactid,tab):
        """
        Insert contacts
        @param interactid: Interaction key
        @param tab: Table of contact
        """ 
        try:
            if(tab):
                if(type(tab[0]) is ListType):
                    for i in tab:
                        t=[interactid]+i
                        self.cur.execute("insert into contact(interactid,srna_begin,srna_end,mrna_begin,mrna_end,score) values(?,?,?,?,?,?)",t)
                else:
                    tab=[interactid]+tab
                    self.cur.execute("insert into contact(interactid,srna_begin,srna_end,mrna_begin,mrna_end,score) values(?,?,?,?,?,?)",tab)
                self.conn.commit()
        except sqlite3.OperationalError:
            commit_success=False
            while(not commit_success):
                time.sleep(1.0)
                commit_success=self.commitRetry()
    
    def getSoft(self,softid):
        """
        Get if software exist by its id
        @param softid: Soft key
        @return: A softid
        """
        self.cur.execute("select softid from soft where softid=?",(softid,))
        result=self.cur.fetchone()
        if result: return result[0]
        else: return None
        
    def getSoftname(self,softid):
        """
        Get if software exist by its id
        @param softid: Soft key
        @return: A software name
        """
        self.cur.execute("select name from soft where softid=?",(softid,))
        result=self.cur.fetchone()
        if result: return result[0]
        else: return None
   
    def getidSoftsbyname(self,name):
        """
        Get id of softwares corresponding to one name
        @param name: name of a software
        @return: list of softid
        """
        self.cur.execute("select softid from soft where name=?",(name,))
        fetch=self.cur.fetchall()
        if(fetch): return np.array(self.listscore(fetch))
        else: return None
    
    def getallRNAlength(self,type_RNA):
        """
        Get all RNA length
        @param type_RNA: Type of RNA
        @return: list of RNA with their length
        """
        self.cur.execute("select name,length from rna where type=?",(type_RNA,))
        result=self.cur.fetchall()
        if(result): return result
        else: return None
    
    def getallRNA(self,type_RNA):
        """
        Get all RNA name from one type
        @param type_RNA: Type of RNA
        @return: list of RNA
        """
        self.cur.execute("select name from rna where type=?",(type_RNA,))
        result=self.cur.fetchall()
        if(result): return result
        else: return None
    
        
    def getRNA(self,name,type_RNA):
        """
        Get an RNA by its name and type
        @param name: RNA name
        @param type: RNA type
        @return: A rnaid
        """
        self.cur.execute("select rnaid from rna where name=? and type=?",(name,type_RNA))
        result=self.cur.fetchone()
        if result: return result[0]
        else: return None
        
    def getsRNAlenbyIntid(self,interactid):
        """
        Get an sRNA length from an interaction
        @param interactid: Interaction id 
        @return: sRNA length
        """
        if interactid:
            self.cur.execute("select length from interact as i, rna as a where i.interactid=? and a.rnaid=i.srnaid",(interactid,))
            result=self.cur.fetchone()
            if(result): return result[0]
            else: return None
        return None
    
    def getAllsRNAlenbyIntid(self,table_interactid):
        """
        Get an sRNA length from an interaction
        @param table_interactid: table of interaction id 
        @return: list of sRNA length
        """
        result=[]
        for interactid in table_interactid:
            self.cur.execute("select length from interact as i, rna as a where i.interactid=? and a.rnaid=i.srnaid",interactid)
            fetch=self.cur.fetchone()
            if(fetch): result+=[fetch[0]]
            else:
                raise ValueError("sRNA %s is missing in RNA description"%(self.getsRNAbyIntid(interactid[0])))
                sys.exit()
        return result
    
    def getsRNAidbyIntid(self,interactid):
        """
        Get an sRNAid from an interaction
        @param interactid: Interaction id
        """
        if(interactid):
            self.cur.execute("select srnaid from interact where interactid=?",(interactid,))
            result=self.cur.fetchone()
            if(result): return result[0]
            else: return None
        return None
    
    def getAllsRNAidbyIntid(self,table_interactid):
        """
        Get an sRNAid from an interaction
        @param interactid: Interaction id
        @return: List of sRNAid
        """
        i=0
        result=np.arange(len(table_interactid))
        for interactid in table_interactid:
            self.cur.execute("select srnaid from interact where interactid=?",interactid)
            fetch=self.cur.fetchone()
            if(fetch): result[i]=fetch[0]
            else: result[i]=-1
            i+=1
        return result
    
    def getmRNAlenbyIntid(self,interactid):
        """
        Get an mRNA from an interaction
        @param interactid: Interaction id 
        @return: mRNA length
        """
        if interactid:
            self.cur.execute("select length from interact as i, rna as a where i.interactid=? and a.rnaid=i.mrnaid",(interactid,))
            result=self.cur.fetchone()
            if(result): return result[0]
            else: return None
        return None
    
    def getAllmRNAlenbyIntid(self,table_interactid):
        """
        Get an mRNA from an interaction
        @param table_interactid: table of interaction id 
        @return: List of mRNA length
        """
        result=[]
        for interactid in table_interactid:
            self.cur.execute("select length from interact as i, rna as a where i.interactid=? and a.rnaid=i.mrnaid",interactid)
            fetch=self.cur.fetchone()
            if(fetch): result+=[fetch[0]]
            else:
                raise ValueError("mRNA %s is missing in RNA description"%(self.getmRNAbyIntid(interactid[0])))
                sys.exit()
        return result
    
    def getNbInteract(self,softid):
        """
        Count the number of interaction for one software
        @param softid: soft key
        @return: Number of interaction
        """
        self.cur.execute("select count(*) from interact where softid=?",(softid,))
        result=self.cur.fetchone()
        if(result): return result[0]
        else: return None 
        
        
    def getInteract(self,softid):
        """
        Get interactions of one software
        @param softid: soft key
        @return: List of interactid
        """
        self.cur.execute("select interactid from interact where softid=?",(softid,))
        result=self.cur.fetchall()
        if(result): return result
        else: return None
        
    def getScore(self,interactid):
        """
        Get the score of an interaction
        @param interactid: interaction key
        @return: List of score
        """
        self.cur.execute("select score from contact where interactid=?",(interactid,))
        fetch=self.cur.fetchall()
        if(fetch): return np.array(self.listscore(fetch))
        else: return None
        
    def getAllScore(self,table_interactid):
        """
        Get the score for several interactions
        @param table_interactid: list of interaction key
        @return: List of score
        """
        result=[]
        for interactid in table_interactid:
            self.cur.execute("select score from contact where interactid=?",interactid)
            fetch=self.cur.fetchall()
            if(fetch): result+=[np.array(self.listscore(fetch))]
            else: result+=[None]
        return result
    
    def getAllSoft(self):
        """
        Get all the software
        @return: List of software
        """
        self.cur.execute("select * from soft")
        result=self.cur.fetchall()
        if(result): return result
        else: return None
    
    def getsRNAbyIntid(self,interactid):
        """
        Get the name of an sRNA from interaction key
        @param interactid: interaction key
        @return: name of the sRNA
        """
        self.cur.execute("select name from interact as i, rna as a where i.interactid=? and a.rnaid=i.srnaid",(interactid,))
        result=self.cur.fetchone()
        if(result): return result[0]
        else: return None
        
    def getmRNAbyIntid(self,interactid):
        """
        Get the name of an mRNA from interaction key
        @param interactid: interaction key
        @return: name of the mRNA
        """
        self.cur.execute("select name from interact as i, rna as a where i.interactid=? and a.rnaid=i.mrnaid",(interactid,))
        result=self.cur.fetchone()
        if(result): return result[0]
        else: return None
    
    def listit(self,t):
        """
        Convert fetchall result into list of list
        @return: List of list
        """
        return list(map(self.listit, t)) if isinstance(t, (list, tuple)) else t
    
    def listscore(self,score):
        """
        Unlist score
        @return: List of score
        """
        if(score): return [i[0] for i in score]
        return None
            
    
    def getPositionsbyIntid(self,interactid):
        """
        Get position for one interaction
        @param interactid: interact key
        @return: List of position
        """
        self.cur.execute("select srna_begin, srna_end, mrna_begin, mrna_end from contact where interactid=?",(interactid,))
        result=self.cur.fetchall()
        if(result): return result
        else: return None
        
    def getallsRNA(self):
        """
        Get all sRNA
        @return: List of rnaid
        """
        self.cur.execute("select rnaid from rna where type=0")
        result=self.cur.fetchall()
        if result: return result
        else: return None
    
    def getallsRNAname(self):
        """
        Get an sRNA name and unlist
        @return: List of rnaid
        """
        self.cur.execute("select name from rna where type=0")
        result=self.cur.fetchall()
        if result: return self.listscore(result)
        else: return None
        
    def getIntbysRNAid(self,srnaid):
        """
        Get the interactid for an sRNA
        @param srnaid: rna key
        @return: List of interactid
        """
        self.cur.execute("select interactid from interact where srnaid=?",(srnaid,))
        result=self.cur.fetchall()
        if(result): return result
        else: return None
        
    def IndexInteractidOnContact(self):
        """
        Create index on interactid for contact
        """
        self.cur.execute("create index if not exists contact_interactid on contact(interactid)")
        self.conn.commit()
        
    def IndexSrnaidOnInteract(self):   
        """
        Create index on srnaid for interact
        """
        self.cur.execute("create index if not exists interact_srnaid on interact(srnaid)")
        self.conn.commit()
        
    def IndexMrnaidOnInteract(self):
        """
        Create index on mrnaid for interact
        """    
        self.cur.execute("create index if not exists interact_mrnaid on interact(mrnaid)")
        self.conn.commit()
            
    def createIndexes(self):
        """
        Create index for most use connection
        """
        self.IndexInteractidOnContact()
        self.IndexSrnaidOnInteract()
        self.IndexMrnaidOnInteract()