#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 @brief: Create sRNA and mRNA multifasta files
 @author: Amine Ghozlane
 @version: 1.0 
"""
from Files import *
from Genbank import *
from Sequences import *
from RNA import *

def main():
    """
    Main program function
    """
    srna_data=None
    seq=Sequences()
    #Get the arguments
    args=Files()
    #Add genbank data
    seq.parse(Genbank(args.genbank))
    #Get sRNA
    if(args.sRNA_file):
        srna_data=RNA()
        srna_data.readRNA(args.sRNA_file)
        seq.writesRNA(args.results, srna_data)
    #Write output
    if(args.mRNA_file):
        mrna_data=RNA()
        mrna_data.readRNA(args.mRNA_file)
        seq.list_genes+=mrna_data.RNA_obj
    seq.writemRNA(args.begin,args.end,args.results,args.complete)
    seq.writeCorresponding(args.results)
    
if __name__ == "__main__":
    main()
