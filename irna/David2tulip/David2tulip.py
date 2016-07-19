#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 @brief: Create tulip files
 @author: Amine Ghozlane
 @version: 1.0 
"""
import os, sys, csv
from Graph import *
from Similarity import *
from pValue import *
from Multilist import *
from Interact import *
from David import *
from Davidconfig import *
from Database_data import *
from Files import *
from Name import *
#import argparse
### Main program ###

#def usage(parser,args):
#    """
#    Test correct usage of arguments
#    @param parser: Parser object
#    @param args: Arguments
#    """
#    #Run without arguments
#    if len(sys.argv)== 1:
#        parser.print_usage()
#        sys.exit()
#    #Test multilist file argument
#    if args.multilist:
#        if(not os.path.isfile(args.multilist)):
#            print("Error with \"%s\" : -m required a multilist file\n"%args.multilist)
#            parser.print_help()
#            sys.exit()            
#    #Test result repertory argument
#    if args.results:
#        if(not os.path.isdir(args.results)):
#            print("Error with \"%s\" : -r required a repertory\n"%args.results)
#            parser.print_help()
#            sys.exit()
#    #Test pValue file argument
#    if args.pValue:
#        if(not os.path.isfile(args.pValue)):
#            print("Error with \"%s\" : -p required a file\n"%args.pValue)
#            parser.print_help()
#            sys.exit()
#    #Test enrichment repertory argument
#    if args.enrichment:
#        if(not os.path.isfile(args.enrichment)):
#            print("Error with \"%s\" : -e required a file\n"%args.enrichment)
#            parser.print_help()
#            sys.exit()
#    #Test filter file
#    if args.filter:
#        if(not os.path.isfile(args.filter)):
#            print("Error with \"%s\" : -f required a file\n"%args.filter)
#            parser.print_help()
#            sys.exit()
#    #Test filter file
#    if args.similarity:
#        if(not os.path.isfile(args.similarity)):
#            print("Error with \"%s\" : -s required a file\n"%args.similarity)
#            parser.print_help()
#            sys.exit()
#    #Test filter file
#    if args.interact:
#        if(not os.path.isfile(args.interact)):
#            print("Error with \"%s\" : -s required a file\n"%args.interact)
#            parser.print_help()
#            sys.exit()
#    #Test david arguments
#    if args.enrichment and args.david:
#        print("Warning : enrichment file \"%s\" will be ignored\n"%args.enrichment)
#    #Test iRNA db
#    if args.iRNA:
#        if(not os.path.isfile(args.iRNA)):
#            print("Error with \"%s\" : -s required a file\n"%args.iRNA)
#            parser.print_help()
#            sys.exit()
#            
##Determine les fichiers fournis en arguments
#def getArgument():
#    """
#    Determine the argument
#    @return: arguments
#    """
#    #Parsing arguments
#    parser = argparse.ArgumentParser(description='Create nodes and edge csv for Tulip from multiple datasets.')
#    parser.add_argument('-m', '--multilist',help='Multilist file from iRNA',required=True)
#    parser.set_defaults(pValue=None)
#    parser.add_argument('-l', '--filter',help='Filter multilist for only known metabolites and use it for analysis')
#    parser.add_argument('-s','--similarity',help='Set similarity of sRNA groups')
#    parser.add_argument('-p', '--pValue',help='pValue of selected files')
#    parser.add_argument('-r', '--results',help='Path to result repertory',required=True)
#    parser.add_argument('-d', '--david',help='Submit multilist to DAVID',action='store_true')
#    parser.add_argument('-i', '--interact',help='Path to interaction file')
#    parser.add_argument('-e', '--enrichment',help='Path to DAVID enrichment file')
#    parser.add_argument('-y','--iRNA',help='Path to iRNA db file')
#    parser.add_argument('-f', '--fastmode',help='Use sqlitebck package to get faster sqlite implementation',action='store_true')
#    args = parser.parse_args()
#    
#    #Verify usage
#    usage(parser,args)
#    return args

def main():
    """
    Main program function
    """
    #Get the arguments
    #args=getArgument()
    args=Files()
    #Initiate graph object
    mygraph=Graph()
    #If the uniprot2go filter is done
    if(args.filter):
        #Intanciate a multilist object
        mymultilist=Multilist(args.multilist)
        #filter data
        mymultilist.GOfilter(args.filter,args.results)
        #Change multilist file
        args.multilist=mymultilist.filtered
    #If multilist file is given
    if(args.multilist):
        config=Davidconfig(args.config)    
        #Parse multilist data
        mygraph.parse(Multilist(args.multilist))  
        #Parse pValue data
        if(args.pValue): mygraph.parse(pValue(args.pValue))
        #Parse similarity data
        if(args.similarity): mygraph.parse(Similarity(args.similarity))
        #Parse interaction data
        if(args.interact): mygraph.parse(Interact(args.interact))
        if(args.iRNA_db): mygraph.parse(Database_data(args.iRNA_db,args.fastmode))
        #Analysis david results from file
        if(args.enrichment):
            #Get enrichment data
            enrichment=config.readconfig(David(args.enrichment))
            #Set data on the graph
            mygraph.parse(enrichment)
        #Analyse multilist with david
        elif(args.david):
            #Set config data
            enrichment=config.readconfig(David())
            #Get david data
            enrichment.analysis(mygraph)
            #Write data
            enrichment.writefile(args.results)
            #Set data on the graph
            mygraph.parse(enrichment)
        if(args.name): mygraph.parse(Name(args.name))
        #Generate nodes csv
        mygraph.writefile(args.results)
    
if __name__ == "__main__":
    main()