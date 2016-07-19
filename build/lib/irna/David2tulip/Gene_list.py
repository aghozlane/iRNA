"""
 @brief: Handle DAVID enrichment data
 @author: Amine Ghozlane
 @version: 1.0 
"""
class Gene_list:
    
    def __init__(self,sudsobject,sRNA):
        """
        Instanciate Gene_list object
        @param sudsobject: suds object
        @param sRNA: Name of the sRNA 
        """
        self.sRNA=sRNA
        #self.EASEBonferroni=sudsobject.EASEBonferroni
        #self.afdr=sudsobject.afdr
        #self.benjamini=sudsobject.benjamini
        #self.bonferroni=sudsobject.bonferroni
        self.categoryName=sudsobject.categoryName
        #self.ease=sudsobject.ease
        #self.fisher=sudsobject.fisher
        #self.foldEnrichment=sudsobject.foldEnrichment
        self.geneIds=sudsobject.geneIds
        #self.id=sudsobject.id
        #self.listHits=sudsobject.listHits
        #self.listName=sudsobject.listName
        #self.listTotals=sudsobject.listTotals
        #self.percent=sudsobject.percent
        #self.popHits=sudsobject.popHits
        #self.popTotals=sudsobject.popTotals
        #self.rfdr=sudsobject.rfdr
        #self.scores=sudsobject.scores
        self.termName=sudsobject.termName