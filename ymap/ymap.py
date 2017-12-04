#!/usr/bin/env python
#
# year of release ! 2016
#
# With the exponential growth of Posttranslational modifications (PTMs) 
# data and and lack of characterisation of all the PTM-types. Its important
# to undersand properly the functions and experimental relevence of PTMs by 
# creating the tools that facilitate the PTMs based analyses. 
# And to understand the importance of PTMs in Yeast genome, its important 
# to make it easier to map experimental mutations to PTM positional data. 
# it's also important and relevent to translate genetic abrretions to understand 
# the phenotype. 
# We architect a python (yMap) library to help users to understand which parts of
# mutated proteins are affected during the yeast experimentation.
# This facilitation not only would help bioligists to interpret their data 
# efficiently but also gives freedom to save time by mapping data to mutations 
# easily 
#
# The yMap program is a python based fast and robust automated method to map 
# large yeast variants to proteins post-translational modifications, proteins domains,
# proteins-DNA binding domains, proteins structural regions, proteins active and 
# binding sites, proteins networks visualisation. 
# For Usage see README file
#
# Dependencies:
# Orange Bioinformatics 
# see README file 
#

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals
try:
    from builtins import next
    from builtins import str
    from builtins import range
    from builtins import object
    from builtins import bytes
except ImportError:
    pass
import os
import sys
import math
import zipfile
from itertools import groupby
import shutil
import time
import urllib
from decimal import Decimal # ADDED
from collections import OrderedDict
import webbrowser
try:
    import Orange
except ImportError:
    import Orange3
from orangecontrib.bio import go    
from six.moves import range
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
import pkg_resources
from pkg_resources import resource_stream
import re

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ Mutation type (Synon | Non-Synon | Stop codon) module (see exmple data) \\\\\\\\\\\\\\\\\\\\\
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


def translate_dna(dna):
    """ calculate the start position for the final codon """
    last_codon_start = len(dna) - 2 
    protein = "" 
    # process the dna sequence in three base chunks
    for start in range(0,last_codon_start,3): 
        codon = dna[start:start+3]
        aa = genetic_code.get(codon.upper(), 'X') 
        protein = protein + aa 
    return protein 


def revcomp(dna, reverse=True, complement=True):
    """ reverse complement of a protein in negative strand"""
    bases = 'ATGCTACG'
    complement_dict = {bases[i]:bases[i+4] for i in range(4)}
    if reverse:
        dna = reversed(dna)
        result_as_list = None
    if complement:
        result_as_list = [complement_dict[base] for base in dna]
    else:
        result_as_list = [base for base in dna]
    return ''.join(result_as_list)


def mutation_file(mut_gene_input, gff_input, d_id_input, mut_prot_output):
        """ defines the mutation types; either Non-Synonmous or Stop Codon"""
        with open(mut_prot_output, 'wb') as t:   
            with open(mut_gene_input, 'r') as mut: 
                for m in mut:
                    m = m.rstrip().split()
                    with open(d_id_input, 'r') as id:    
                        for i in id:
                            i = i.rstrip().split()
                            if not m[0].startswith('c'.upper()):
                                if len(m) != 5  or not m[0].startswith('c'.lower()):
                                    raise StopIteration('Please correct the format of input mutation file')
                                else:
                                    if m[4] == i[2]:
                                        take = m[4]+'\t'+m[0]+'\t'+i[3]+'\t'+m[1]+'\t'+i[4]+'\t'+m[2]+'\t'+m[3]+'\t'+i[5]
                                        take1= take.rstrip().split()       
                                        with open(gff_input, 'r') as orf:     
                                            linee = orf.readlines()[23078:]
                                            up = (x[1] for x in groupby(linee, lambda line: line[0] == ">")) 
                                            for head in up:
                                                head = next(head)[1:].strip()
                                                seq = "".join(s.strip() for s in next(up))
                                                if head == take1[1] and take1[0] == i[2] and take1[7] == '-':   
                                                    cod = 1 + (int(take1[4])-int(take1[3]))                       
                                                    cc = math.ceil(int(cod)/float(3))
                                                    c = str(cc).split('.')                         
                                                    cn = int(c[0])-1    
                                                    sli_n = seq[int(take1[2]):int(take1[4])]                  
                                                    rev_sli_n = revcomp(sli_n, reverse=True, complement=True)  
                                                    sli_m_n = sli_n[:int(-cod)]+take1[6]+sli_n[int(-cod)+1:] 
                                                    rev_sli_m_n = revcomp(sli_m_n, reverse=True, complement=True)   
                                                    wild_type_rev_n = translate_dna(rev_sli_n)                
                                                    mut_type_n = translate_dna(rev_sli_m_n)
                                                    try:
                                                        if wild_type_rev_n[cn] != mut_type_n[cn] and mut_type_n[cn] == '_':
                                                            pic = take1[0]+'\t'+str(c[0])+'\t'+wild_type_rev_n[cn]+'\t'+mut_type_n[cn]+'\t'+'Stop' +'\t'+take1[1]+'\t'+take1[3]
                                                            if pic > str(0): 
                                                                t = open(mut_prot_output, 'a')
                                                                t.write(pic+'\n')
                                                    except IndexError as e:
                                                        pic1 =  take1[0]+ '\t'+ 'Error:'+'\t'+ str(e)
                                                        t = open(mut_prot_output, 'a+')
                                                        t.write(pic1+'\n')
                                                        continue
                                                    try:
                                                        if wild_type_rev_n[cn] != mut_type_n[cn] and mut_type_n[cn] != '_':
                                                            pic = take1[0]+'\t'+str(c[0])+'\t'+wild_type_rev_n[cn]+'\t'+mut_type_n[cn]+'\t'+'Non-Synonymous' +'\t'+take1[1]+'\t'+take1[3]                                                                                                                    
                                                            if pic > str(0):
                                                                t = open(mut_prot_output, 'a+')
                                                                t.write(pic+'\n')
                                                    except IndexError as e:
                                                        pic1 =  take1[0]+ '\t'+ 'Error:'+'\t'+ str(e)
                                                        t = open(mut_prot_output, 'a+')
                                                        t.write(pic1+'\n')
                                                        continue
                                                if head == take1[1] and take1[0]==i[2] and take1[7] == '+':
                                                    code = int(take1[3])-int(take1[2])
                                                    code1 = 1 + (int(take1[3])-int(take1[2])) 
                                                    cce = math.ceil(int(code1)/float(3)) 
                                                    ce = str(cce).split('.') 
                                                    cp = int(ce[0])-1                  
                                                    pos = int(take1[2]) - 1                               
                                                    sli_p = seq[int(pos):int(take1[4])]                   
                                                    sli_m_p = sli_p[:int(code)]+take1[6]+sli_p[int(code)+1:]  
                                                    wild_type_p = translate_dna(sli_p)                    
                                                    mut_type_p = translate_dna(sli_m_p)
                                                    try:                   
                                                        if wild_type_p[cp] != mut_type_p[cp] and mut_type_p[cp] != '_': 
                                                            pick = take1[0]+'\t'+str(ce[0])+'\t'+wild_type_p[cp]+'\t'+mut_type_p[cp]+'\t'+'Non-Synonymous'+'\t'+take1[1]+'\t'+take1[3]
                                                            if pick > str(0):
                                                                with open(mut_prot_output, 'a+') as t:
                                                                    t.write(pick+'\n')
                                                    except IndexError as e:
                                                        pic1 =  take1[0]+ '\t'+ 'Error:'+'\t'+ str(e)
                                                        t = open(mut_prot_output, 'a+')
                                                        t.write(pic1+'\n')
                                                        continue
                                                    try:
                                                        if wild_type_p[cp] != mut_type_p[cp] and mut_type_p[cp]=='_':
                                                            pick = take1[0]+'\t'+str(ce[0])+'\t'+wild_type_p[cp]+'\t'+mut_type_p[cp]+'\t'+'Stop' +'\t'+take1[1]+'\t'+take1[3]
                                                            if pick > str(0):
                                                                with open(mut_prot_output, 'a+') as t:
                                                                    t.write(pick+'\n')
                                                    except IndexError as e:
                                                        pic1 =  take1[0]+ '\t'+ 'Error:'+'\t'+ str(e)
                                                        t = open(mut_prot_output, 'a+')
                                                        t.write(pic1+'\n')


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# //////////////////                     UniProt data                 /////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def parse_gene_names(yeastID_input):
    """Parse the yeastID_input file (mapping UniProt IDs to gene names) into a dictionary."""
    gene_names = {}
    with open(yeastID_input, 'r') as proteins:
        next(proteins) # skip the header
        for line in proteins:
            common_names = []
            sgd_names = []
            uniprot_id, sgd_name, common_name = line.rstrip('\n').split('\t')
            if sgd_name == '': # Replace empty strings with a missing value string e.g. NA
                sgd_name = 'NA'
            if common_name == '':
                common_name = 'NA'
            if ';' in common_name: # Check for multiple common names for a given UniProt ID
                common_names = common_name.split('; ')
                common_names = ['NA' if name == '' else name for name in common_names] # Replace empty strings with a missing value string e.g. NA
            if ';' in sgd_name: # Check for multiple locus names for a given UniProt ID
                sgd_names = sgd_name.split('; ')
                sgd_names = ['NA' if name == '' else name for name in sgd_names] # Replace empty strings with a missing value string e.g. NA
            if common_names or sgd_names:
                gene_names[uniprot_id] = list(zip(common_names, sgd_names)) # Create a list of tuples, each a (common name, sgd name) pair
            else:        
                gene_names[uniprot_id] = [(common_name, sgd_name)]
    return gene_names


def parse_mutations(mut_prot_input):
    """Parse the mut_prot_input file (mapping mutated proteins to mutation positions) into a dictionary."""
    mutated_proteins = {}
    with open(mut_prot_input, 'r') as mutations:
        #TODO: Write a check for a header - perhaps specify parameter header=T, like some R functions (give user flexibility) 
        next(mutations) # Skip the header
        for line in mutations:
            common_name, mutation_pos = line.rstrip('\n').split('\t')[:2] # Take index to work with mutation input file with >2 columns
            if common_name not in mutated_proteins:
                mutated_proteins[common_name] = {mutation_pos} # Mutation positions put into set to remove duplicates
            else:
                mutated_proteins[common_name].add(mutation_pos)
    return mutated_proteins

def parse_biogrid(uniprot_biogrid_input): 
    """Return dictionary mapping UniProt IDs to BioGrid IDs.""" 
    uniprot_biogrid_map = {}
    with open(uniprot_biogrid_input, 'r') as uniprot_biogrid:
        for line in uniprot_biogrid:
            uniprot_id, biogrid_ids = line.rstrip(';\n').split('\t')
            biogrid_ids = biogrid_ids.split(';')
            uniprot_biogrid_map[uniprot_id] = biogrid_ids
    return uniprot_biogrid_map

def gff(gff_output):
    """Downloads the current General Feature Format (GFF) file for the Saccharomyces cerevisiae genome"""
    yeast_gff = urlopen('http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff')
    with open(gff_output,'wb') as output_file:
        output_file.write(yeast_gff.read())
    

def frmt(gff_input, frmt_output):
    """Reformat a General Feature Format (GFF) file into a tab-delimited file (tsv) with columns gene id, start and end with strand orientation."""
    id_regex = re.compile(r'ID=(.*?);')
    lines = []
    with open(gff_input, 'r') as gff:
        for line in gff:
            if line.startswith('###'):
                break
            elif line.startswith(('##', '#')):
                continue
            chr, source, feature, start, end, score, strand, frame, attribute = line.split()
            if feature == 'gene':
                gene_id = id_regex.match(attribute).group(1)
                new_line = '\t'.join([gene_id, start, end, strand]) + '\n'
                lines.append(new_line)
    with open(frmt_output,'w') as parsed_gff:
        parsed_gff.writelines(lines)


# IMPORTANT TO NOTE - UniProt IDs are not unique in the output file of this function. See P02994 for example.
# This is because some UniProt IDs map to more than one locus. 
def id_map(yeastID_input, frmt_input, d_id_map_output):
    """Map Uniprot IDs to gene names and genomic loci (start, end positions and strand orientation)."""
    gene_names = {}
    lines = []
    for uniprot_id, names in yeastID_input.items():
        for common_name, sgd_name in names:
            gene_names[sgd_name] = [uniprot_id, sgd_name, common_name]
    with open(frmt_input, 'r') as gff:
        #TODO: Write a check for a header - perhaps specify parameter header=T, like some R functions (give user flexibility) 
        next(gff) # Skip the header
        for line in gff:
            sgd_name, start, end, strand = line.rstrip('\n').split('\t')
            new_line = gene_names[sgd_name].copy()
            new_line.extend([start, end, strand])
            new_line = '\t'.join(new_line) + '\n'
            lines.append(new_line)
    with open(d_id_map_output, 'w') as protein_loci_mapped:
        protein_loci_mapped.writelines(lines)


def pTMdata(uniprot_output):
    """Download up-to-date UniProt data file using the following query settings: Species: Saccharomyces cerevisiae, ...""" #TODO: Update documentation with query settings
    #TODO: Make the uniprot URL cleaner/more understandable. See http://www.uniprot.org/help/api_queries
    uniprot = urlopen('http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=gff&columns=id,feature(MODIFIED%20RESIDUE)')
    with open(uniprot_output,'wb') as output_file:
        output_file.write(uniprot.read())


def make_ptms_file(uniprot_input, ptms_output):
    """Extract Post-Translational Modification (PTM) data from a downloaded Uniprot file and write to tab-delimited (tsv) file."""
    annotation_regex = re.compile(r'Note=(.*?)(;|$)') # Alternation needed for entries that only contain a "Note" in the attribute string (see below)
    lines = []
    with open(uniprot_input, 'r') as uniprot:
        for line in uniprot:
            if not line.startswith('##'): #TODO: Decide whether this implementation is better than "if line.startswith('##'): continue"
                uniprot_id, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split('\t')
                if feature in ('Lipidation', 'Glycosylation', 'Modified residue', 'Cross-link'):
                    modification = annotation_regex.match(attribute).group(1)
                    new_line = '\t'.join([uniprot_id, start, modification]) + '\n'
                    lines.append(new_line)
    with open(ptms_output, 'w') as ptms:
        ptms.writelines(lines)


def iD(yeastID_output):
    """Download up-to-date UniProt data file using the following query settings: Species: Saccharomyces cerevisiae, ...""" #TODO: Update documentation with query settings
    #TODO: Make the uniprot URL cleaner/more understandable. See http://www.uniprot.org/help/api_queries
    uniprot = urlopen('http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=tab&columns=id,genes(OLN),%2Cgenes(PREFERRED)')
    with open(yeastID_output,'wb') as output_file:
        output_file.write(uniprot.read())

    
def pmap(yeastID_input, ptms_input, ptm_id_output):          
    """Map UniProt IDs to protein names and post-translational modifications (PTMs).
    
    Arguments:
    yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names
    ptms_input -- file path, mapping UniProt ID to a PTM and its position in polypeptide
    ptm_id_output -- file path, mapping UniProt ID, ordered locus name, common name, PTM position and PTM 
    """
    lines = [] 
    with open(ptms_input, 'r') as ptms:
        for line in ptms:
            uniprot_id, ptm_pos, modification = line.rstrip('\n').split('\t')
            for common_name, sgd_name in yeastID_input[uniprot_id]:
                new_line = [uniprot_id, sgd_name, common_name, ptm_pos, modification]
                new_line = '\t'.join(new_line) + '\n'
                lines.append(new_line)
    with open(ptm_id_output, 'w') as ptm_id:
        ptm_id.writelines(lines)
                                         
#TODO: Still don't like this implementation, but it's ok, I guess
def ptm_map(mut_prot_input, ptm_id_input, mapped_ptms_output, summary_output):
    """Map the positions of mutations in mutated proteins to post-translational modifications (PTMs).
    
    Arguments:
    mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein
    ptm_id_input -- file path, mapping UniProt ID, ordered locus name, common name, PTM position and PTM
    mapped_ptms_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, PTM position and PTM
    summary_output -- file path, recording a summary of all features to which mutations have been mapped
    """
    mapped_ptms_lines = []
    summary_lines = []
    with open(ptm_id_input, 'r') as ptms:
        for line in ptms:
            uniprot_id, sgd_name, common_name, ptm_pos, ptm = line.rstrip('\n').split('\t')
            for mutated_protein, mutated_positions in mut_prot_input.items():
                if mutated_protein == common_name:
                    for mut_pos in mutated_positions:
                        if mut_pos == ptm_pos:
                            mapped_ptm_line = '\t'.join([uniprot_id, sgd_name, common_name, ptm_pos, ptm, 'UniProt']) + '\n'
                            summary_line = '\t'.join([uniprot_id, sgd_name, common_name, ptm_pos, ptm, 'PTMs', 'UniProt']) + '\n'
                            mapped_ptms_lines.append(mapped_ptm_line)
                            summary_lines.append(summary_line)
    with open(mapped_ptms_output, 'w') as mapped_ptms:
        mapped_ptms.writelines(mapped_ptms_lines)
    with open(summary_output, 'a') as summary:
        summary.writelines(summary_lines)


#TODO Find out how important PROSITE references are? What are they used for?
def make_domains_file(uniprot_input, domains_output): 
    """Extract domain data from a downloaded UniProt file and write to tab-delimited (tsv) file."""
    annotation_regex = re.compile(r'Note=(.*?)(;|$)') # Alternation needed for entries that only contain a "Note" in the attribute string (see below)
    prosite_regex = re.compile(r'PROSITE(.*?)(\||$)') # Note: sometimes 'PROSITE:' and sometimes 'PROSITE-ProRule:' (Only last entry 'E9PAE3' apparently)
    lines = []
    with open(uniprot_input, 'r') as uniprot:
        for line in uniprot:
            if not line.startswith('##'): #TODO: Decide whether this implementation is better than "if line.startswith('##'): continue"
                uniprot_id, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split('\t')
                if feature == 'Domain':
                    domain = annotation_regex.match(attribute).group(1) #TODO: Check whether this is actually a domain name or if it is the common name of the protein
                    prosite_ref = ''
                    prosites_matched = prosite_regex.finditer(attribute) # Find all PROSITE references in the attribute column (sometimes there are more than one)
                    for prosite_match in prosites_matched:
                        prosite_ref += prosite_match.group() #TODO: Is this what we want or do we only want group(1) e.g. something like 'PRU00457'
                    prosite_ref = prosite_ref.rstrip('|')
                    new_line = '\t'.join([uniprot_id, start, end, domain, prosite_ref]) + '\n'
                    lines.append(new_line)
    with open(domains_output, 'w') as domains:
        domains.writelines(lines)
            

def d_map(yeastID_input, domains_input, domain_id_output):
    """Map UniProt IDs to protein names and domains and PROSITE references
    
    Arguments:
    yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names
    domains_input -- file path, mapping UniProt ID to domains, their start and end positions in the polypeptide and PROSITE references 
    domain_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and domain data 
    """
    lines = [] 
    with open(domains_input, 'r') as domains:
        for line in domains:
            uniprot_id, start, end, domain, prosite_ref = line.rstrip('\n').split('\t')
            for common_name, sgd_name in yeastID_input[uniprot_id]:
                new_line = [uniprot_id, sgd_name, common_name, start, end, domain, prosite_ref]
                new_line = '\t'.join(new_line) + '\n'
                lines.append(new_line)
    with open(domain_id_output, 'w') as domain_id:
        domain_id.writelines(lines)


def dmap(mut_prot_input, domain_id_input, mapped_domains_output, summary_output):
    """Map the positions of mutations in mutated proteins to domains.
    
    Arguments:
    mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein
    domain_id_input -- file path, mapping UniProt ID, ordered locus name, common name, domain position (start and end) and domain name
    mapped_domains_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, domain position and domain name
    summary_output -- file path, recording a summary of all features to which mutations have been mapped
    """
    mapped_domains_lines = []
    summary_lines = []
    with open(domain_id_input, 'r') as domains:
        for line in domains:
            uniprot_id, sgd_name, common_name, domain_start, domain_end, domain_name, prosite_ref = line.rstrip('\n').split('\t')
            for mutated_protein, mutated_positions in mut_prot_input.items():
                if mutated_protein == common_name:
                    for mut_pos in mutated_positions:
                        if int(domain_start) <= int(mut_pos) <= int(domain_end): # Check if the mutation lies in a domain of the mutated protein
                            mapped_domain_line = '\t'.join([uniprot_id, sgd_name, common_name, domain_start, mut_pos, domain_end, domain_name, prosite_ref, 'UniProt']) + '\n'
                            summary_line = '\t'.join([uniprot_id, sgd_name, common_name, domain_start, mut_pos, domain_end, domain_name, prosite_ref, 'Domains', 'UniProt']) + '\n'
                            mapped_domains_lines.append(mapped_domain_line)
                            summary_lines.append(summary_line)
    with open(mapped_domains_output, 'w') as mapped_domains:
        mapped_domains.writelines(mapped_domains_lines)
    with open(summary_output, 'a') as summary:
        summary.writelines(summary_lines)

#TODO: The bottleneck in the program is the GO enrichment analysis. We have only a 9 genes in the test file 'mapped_domains.txt' and it takes 
# between 14s and 20s to perform the enrichment and about the same time to load the ontology and annotations.
# Perhaps we should have an option on the website where the user can select whether or not they want to perform GO enrichment.
def enrich(mapped_mut_input):
    """Perform Gene Ontology (GO) enrichment on a set of genes (obtained from the second column of an tab-delimited input file).
    
    Write the results of the enrichment to a file. If a GO term has been found to be overrepresented, then its ID, name, the 
    significance of the enrichment (for further details on how p-value is calculated see ?????????????), a list of the genes 
    in the gene set that are annotated with the enriched GO term, and a reference count (???????) are written as a line in 
    the file.
    """
    lines = []
    p_value_output = os.path.join(os.path.dirname(mapped_mut_input), p_value_file)
    with open(mapped_mut_input, 'r') as mapped_mutations:
        sgd_gene_names = [line.split('\t')[1] for line in mapped_mutations]
    enriched_go_terms = annotations.get_enriched_terms(sgd_gene_names)
    if len(enriched_go_terms) == 0:
        line = 'No enriched GO terms found.'
    else:
        for go_id, (genes, p_value, ref) in enriched_go_terms.items(): #TODO: Find out what the reference count is
            if p_value < 0.05:
                term = ontology[go_id]
                formatted_p_value = '%.2E' % p_value
                genes = ', '.join(genes)
                line_elements = [term.id, term.name, formatted_p_value, genes, str(ref)]
                line = '\t'.join(line_elements) + '\n'
    lines.append(line)
    with open(p_value_output, 'w') as out:
        out.writelines(lines)
                           


def make_bact_file(uniprot_input, bact_output):
    """Extract binding site and active site data from a downloaded UniProt file and write to tab-delimited (tsv) file."""
    annotation_regex = re.compile(r'Note=(.*?)(;|$)') # Alternation needed for entries that only contain a "Note" in the attribute string (see below)
    lines = []
    with open(uniprot_input, 'r') as uniprot:
        for line in uniprot:
            if not line.startswith('##'):
                uniprot_id, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split('\t')
                if feature in ('Binding site', 'Active site'):
                    match = annotation_regex.match(attribute) # Only need to search for regex at the beginning of the string
                    binds_or_activity = match.group(1) if match else ''
                    new_line = '\t'.join([uniprot_id, feature, start, binds_or_activity]) + '\n'
                    lines.append(new_line)
    with open(bact_output, 'w') as bact:
        bact.writelines(lines)
                        
                         
def id(yeastID_input, bact_input, sites_id_output): 
    """Map UniProt IDs to protein names and binding sites and active sites.
    
    Arguments:
    yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names
    bact_input -- file path, mapping UniProt ID to binding sites and active site, their position in the polypeptide and their binding partner/activity
    sites_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and binding/active site data
    """
    lines = [] 
    with open(bact_input, 'r') as bact:
        for line in bact:
            uniprot_id, feature, position, binds_or_activity = line.rstrip('\n').split('\t')
            for common_name, sgd_name in yeastID_input[uniprot_id]:
                new_line = [uniprot_id, sgd_name, common_name, feature, position, binds_or_activity]
                new_line = '\t'.join(new_line) + '\n'
                lines.append(new_line)
    with open(sites_id_output, 'w') as sites_id:
        sites_id.writelines(lines)


def mmap(mut_prot_input, sites_id_input, mapped_sites_output, summary_output):
    """Map the positions of mutations in mutated proteins to binding sites and active sites.
    
    Arguments:
    mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein
    sites_id_input -- file path, mapping UniProt ID, ordered locus name, common name, binding/active site, their position and binding partner/activity 
    mapped_sites_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, binding/active site data
    summary_output -- file path, recording a summary of all features to which mutations have been mapped
    """
    mapped_sites_lines = []
    summary_lines = []
    with open(sites_id_input, 'r') as sites:
        for line in sites:
            uniprot_id, sgd_name, common_name, feature, bact_pos, binds_or_activity = line.rstrip('\n').split('\t')
            for mutated_protein, mutated_positions in mut_prot_input.items():
                if mutated_protein == common_name:
                    for mut_pos in mutated_positions:
                        if mut_pos == bact_pos: # Check if the mutation coincides with the position of a binding/active site
                            mapped_sites_line = '\t'.join([uniprot_id, sgd_name, common_name, feature, mut_pos, binds_or_activity, 'UniProt']) + '\n'
                            summary_line = '\t'.join([uniprot_id, sgd_name, common_name, feature, mut_pos, binds_or_activity, 'Active/Binding sites', 'UniProt']) + '\n'
                            mapped_sites_lines.append(mapped_sites_line)
                            summary_lines.append(summary_line)
    with open(mapped_sites_output, 'w') as mapped_sites:
        mapped_sites.writelines(mapped_sites_lines)
    with open(summary_output, 'a') as summary:
        summary.writelines(summary_lines)


def make_nucleotide_file(uniprot_input, nucleotide_output):
    """Extract nucleotide binding site data from a downloaded UniProt file and write to tab-delimited (tsv) file."""
    annotation_regex = re.compile(r'Note=(.*?)(;|$)') # Alternation needed for entries that only contain a "Note" in the attribute string (see below)
    lines = []
    with open(uniprot_input, 'r') as uniprot:
        for line in uniprot:
            if not line.startswith('##'):
                uniprot_id, source, feature, start, end, score, strand, frame, attribute = line.rstrip().split('\t')
                if feature == 'Nucleotide binding':
                    match = annotation_regex.match(attribute) # Only need to search for regex at the beginning of the string
                    binding_nucleotide = match.group(1) if match else ''
                    new_line = '\t'.join([uniprot_id, feature, start, end, binding_nucleotide]) + '\n'
                    lines.append(new_line)
    with open(nucleotide_output, 'w') as nucleotide:
        nucleotide.writelines(lines)


def n_map(yeastID_input, nucleotide_input, nucleotide_id_output): 
    """Map UniProt IDs to protein names and nucleotide binding sites.
    
    Arguments:
    yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names
    nucleotide_input -- file path, mapping UniProt ID to nucleotide binding sites, their position (start and end) in the polypeptide and their binding nucleotide
    nucleotide_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and nucleotide binding site data
    """
    lines = [] 
    with open(nucleotide_input, 'r') as nucleotide:
        for line in nucleotide:
            uniprot_id, feature, start, end, binding_nucleotide = line.rstrip('\n').split('\t')
            for common_name, sgd_name in yeastID_input[uniprot_id]:
                new_line = [uniprot_id, sgd_name, common_name, feature, start, end, binding_nucleotide]
                new_line = '\t'.join(new_line) + '\n'
                lines.append(new_line)
    with open(nucleotide_id_output, 'w') as nucleotide_id:
        nucleotide_id.writelines(lines)


def nucleotide_map(mut_prot_input, nucleotide_id_input, mapped_nucleotide_output, summary_output):
    """Map the positions of mutations in mutated proteins to nucleotide binding sites.
    
    Arguments:
    mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein
    nucleotide_id_input -- file path, mapping UniProt ID, ordered locus name, common name, the start and end position of nucleotide binding sites and the binding nucleotide
    mapped_nucleotide_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, nucleotide binding site data
    summary_output -- file path, recording a summary of all features to which mutations have been mapped
    """
    mapped_nucleotide_lines = []
    summary_lines = []
    with open(nucleotide_id_input, 'r') as sites:
        for line in sites:
            uniprot_id, sgd_name, common_name, feature, nuc_bind_start, nuc_bind_end, binding_nucleotide = line.rstrip('\n').split('\t')
            for mutated_protein, mutated_positions in mut_prot_input.items():
                if mutated_protein == common_name:
                    for mut_pos in mutated_positions:
                        if int(nuc_bind_start) <= int(mut_pos) <= int(nuc_bind_end): # Check if the mutation coincides with the position of a binding/active site
                            mapped_sites_line = '\t'.join([uniprot_id, sgd_name, common_name, feature, nuc_bind_start, mut_pos, nuc_bind_end, binding_nucleotide, 'UniProt']) + '\n'
                            summary_line = '\t'.join([uniprot_id, sgd_name, common_name, feature, nuc_bind_start, mut_pos, nuc_bind_end, binding_nucleotide, 'Nucleotide binding sites', 'UniProt']) + '\n'
                            mapped_nucleotide_lines.append(mapped_sites_line)
                            summary_lines.append(summary_line)
    with open(mapped_nucleotide_output, 'w') as mapped_nucleotides:
        mapped_nucleotides.writelines(mapped_nucleotide_lines)
    with open(summary_output, 'a') as summary:
        summary.writelines(summary_lines)


def bioGrid(uniprot_biogrid_output):
    """Download a list of UniProt IDs and their associated BioGrid IDs."""
    #TODO: Make the uniprot URL cleaner/more understandable. See http://www.uniprot.org/help/api_queries
    uniprot = urlopen('http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=tab&columns=id,database(BioGrid)')
    with open(uniprot_biogrid_output,'wb') as output_file:
        output_file.write(uniprot.read())
        

def preWeb(uniprot_biogrid_input, mapped_mut_input): 

    """ maps mutations to BioGrid ids """ 
    biog_output = os.path.join(os.path.dirname(mapped_mut_input), biog_file) 
    with open(biog_output, 'w') as out:
        with open(uniprot_biogrid_input, 'r') as fl:
            for f in fl:
                f = f.rstrip().split()
                if len(f) > 1:
                    i = f[1].split(';')
                    take = f[0]+'\t'+i[0]
                    take = take.split()
                    with open(mapped_mut_input, 'r') as pro:
                        for p in pro:
                            p = p.split()
                            if take[0] == p[0]:
                                take2 = take[0]+'\t'+take[1]+'\t'+'UniProt'
                                if take2 > str(0):
                                    with open(biog_output, 'a') as out:
                                        out.write(take2+'\n')


def bweb(biog_input): 

    """ opens the BioGrid db in browser with as many tabs as mutated proteins""" 

    url = 'http://thebiogrid.org/'
    fl = open(biog_input, 'r')
    for f in OrderedDict.fromkeys(fl):
        f = f.split()
        webbrowser.open(url + f[1])


def make_pdb_file(uniprot_input, pdb_output):

    """ Structure data filtration from UniProt"""

    with open(pdb_output, 'w') as stru:
        with open(uniprot_input, 'r') as raw:
            for r in raw:
                if not r.startswith('##'):
                    line = r.split()
                    if line[2] == 'Beta' and len(line[9]) > 1:
                        take = line[9].split('|')
                        take3 = line[0]+'\t'+line[2]+'\t'+line[4]+'\t'+line[5]+'\t'+take[1]
                        if take3 > str(0):
                            with open(pdb_output, 'a') as stru:
                                stru.write(take3+'\n')
                                continue
                    if len(line) > 7 and line[2] == 'Helix' or line[2]=='Turn':
                        if len(line[8])>1:
                            tak = line[8].split('|')
                            tak3 = line[0]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+take[1]
                            if tak3 > str(0):
                                with open(pdb_output, 'a+') as stru:
                                    stru.write(tak3+'\n')


def mu_map(yeastID_input, mut_prot_input, mapped_mutation_pos_output):

    """ mutations proteins mapped to the yeastID file"""

    with open(mapped_mutation_pos_output, 'w') as f:
        with open(mut_prot_input) as mutation_file:
            for a in mutation_file:
                a = a.split()
                with open(yeastID_input) as id:
                    for i in id:
                        i = i.split()
                        if len(i) > 2:
                            if a[0] == i[2]:
                                take = i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+a[1]
                                if take > str(0):
                                    with open(mapped_mutation_pos_output, 'a') as f:
                                        f.write(take+'\n')


def pdb(pdb_input, mut_pos_input, mapped_struct_output, summary_output):

    """ This code maps mutations to the proteins structural regions"""

    with open(mapped_struct_output, 'w') as s:
        with open(pdb_input, 'r') as raw:
            for i in raw:
                i = i.split()
                with open(mut_pos_input) as mu: 
                    for m in mu:
                        m = m.split()
                        if i[0] == m[0]:
                            try:
                                if m[3] == 'Error:':
                                    with open(mapped_struct_output, 'a') as s:
                                        s.write("input file contains error position for" +m[2]+ "protein"+'\n')
                                        continue
                                if int(i[2]) <= int(m[3]) and int(i[3]) >= int(m[3]):
                                    take = m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+i[1]+'\t'+i[2]+'\t'+m[3]+'\t'+i[3]+'\t'+i[4]+'\t'+'UniProt'
                                    if take > str(0):
                                        with open(mapped_struct_output, 'a+') as s:
                                            s.write(take+'\n')
                                        with open(summary_output, 'a+') as summary:
                                            summary.write(m[0]+'\t'+m[2]+'\t'+m[3]+'\t'+i[4]+'\t'+i[1]+'\t'+'UniProt'+'\n')
                            except IndexError:
                                pass


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#/////////////////// Annotated PTMs data from other resources than UniProt (know to play role in PPI and cross-talk) /////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#To get the mutational effects on PPI, PTM based crosstalk, Protein domains, we need to run the following data files from one local dict.; the data
#retrieved from PTMcode 2.0 and PTMfunc, for this reason. To run your own lis tagainst this program, all you need to do is to change the file name in
#the test varible and then you get to go, the output contains the pvalue of the GO terms effected by the mutations and also step wise protein output data
#to interpret your experiment."""
#This frame work contains the advance stage of mapping, where same one code can be used for the mapping to the different 
#PTM types, present at interface  and/or ppi.


def interface(yeastID_input, interface_sites_input, mut_prot_input, mapped_interface_output, summary_output):

    """PTM present at the interface of two proteins and known to play role in interaction (Beltrao et al. Cell 2012)"""
    
    with open(mapped_interface_output, 'w') as out:
        with open(interface_sites_input, 'r') as f:
            for l in f:
                line = l.split()
                if len(line) > 5:
                    take = line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[5]
                    take = take.split()
                    with open(mut_prot_input) as mu:
                        for m in mu:
                            m = m.split()
                            if m[0] == take[1] and m[1] == take[2]:
                                take2 = take[0]+'\t'+take[1]+'\t'+take[2]+'\t'+take[3]+'\t'+'PTMfunc'
                                fi = take2.split()
                                with open(yeastID_input) as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[1]:
                                            with open(summary_output, 'a+') as summary:
                                                summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'Interface'+'\t'+'PTMfunc'+'\n')
                                            if take2 > str(0):
                                                with open(mapped_interface_output, 'a') as out:
                                                    out.write(take2+'\n')
                                                
         

def ppi(yeastID_input, ppi_input, mut_prot_input, mapped_ppi_output, summary_output):

    """ PTM present at the interface of two proteins and known to play role in interaction (PTMfunc; Beltrao et al. Cell 2012)"""

    with open(mapped_ppi_output, 'w') as out:
        with open(ppi_input, 'r') as f:
            for ls in f:
                line = ls.split()
                with open(mut_prot_input) as mu:
                    for m in mu:
                        m = m.split()
                        if len(line) > 7:
                            if m[0] == line[1] and m[1] == line[3]:
                                take = line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+'PTMfunc'
                                fi = take.split()
                                with open(yeastID_input) as i:
                                    for di in i:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            with open(summary_output, 'a+') as summary:
                                                summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+'\t'+'PPI'+'\t'+'PTMfunc'+'\n')
                                            if take > str(0):
                                                with open(mapped_ppi_output, 'a') as out:
                                                    out.write(take+'\n')
                                                    continue
                            if m[0] == line[6] and m[1] == line[3]:
                                take2 = line[6]+'\t'+line[2]+'\t'+line[3]+'\t'+'PTMfunc'
                                fi = take2.split()
                                with open(yeastID_input) as i:
                                    for di in i:
                                        di=di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            with open(summary_output, 'a+') as summary:
                                                summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+'\t'+'PPI'+'\t'+'PTMfunc'+'\n')
                                            if take2 > str(0):
                                                with open(mapped_ppi_output, 'a+') as out:
                                                    out.write(take2+'\n')
                    
    
def withinPro(yeastID_input, within_prot_input, mut_prot_input, mapped_within_prot_output, summary_output):

    """ PTMs (predicted) involved in the crosstalk within a given protein at baker's years (Minguez el 2012)"""

    with open(mapped_within_prot_output, 'w') as file1:
        with open(within_prot_input, 'r') as f:
            for l in f:
                line = l.split()
                if len(line)>19:
                    take = line[15]+'\t'+line[16]+'\t'+line[3]+'\t'+line[17]+'\t'+line[7]+'\t'+line[19]
                    take = take.split()
                    with open(mut_prot_input, 'r') as mu:
                        for m in mu:
                            m = m.split()
                            if m[0] == take[1] and m[1]==take[3]:
                                take2 = take[0]+'\t'+take[1]+'\t'+take[2]+'\t'+take[3]+'\t'+'PTMcode'
                                fi=take2.split()
                                with open(yeastID_input) as id:
                                    for di in id:
                                        di=di.split()
                                        if len(di) > 2 and di[2] == fi[1]:
                                            with open(summary_output, 'a+') as summary:
                                                summary.write(di[0]+'\t'+di[2]+'\t'+fi[3]+'\t'+fi[2]+'\t'+'WithinProtein'+'\t'+'PTMcode'+'\n')
                                            if take2 > str(0):
                                                with open(mapped_within_prot_output, 'a') as file1:
                                                    file1.write(take2+'\n')
                                                    continue
                            if m[0] == take[1] and m[1] == take[5]:
                                take3 = take[0]+'\t'+take[1]+'\t'+take[4]+'\t'+take[5]+'\t'+'PTMcode'
                                fi = take3.split()
                                with open(yeastID_input) as id:
                                    for di in id:
                                        di=di.split()
                                        if len(di) > 2 and di[2] == fi[1]:
                                            with open(summary_output, 'a+') as summary:
                                                summary.write(di[0]+'\t'+di[2]+'\t'+fi[3]+'\t'+fi[2]+'\t'+'WithinProtein'+'\t'+'PTMcode'+'\n')
                                            if take3 > str(0):
                                                with open(mapped_within_prot_output, 'a+') as file1:
                                                    file1.write(take3+'\n')
                         

def betweenPro(yeastID_input, between_prot_input, mut_prot_input, mapped_between_prot_output, summary_output):

    """ PTMs (predicted) involved in the crosstalk in different proteins at baker's years (PTMcode 2.0; Minguez el 2012) """

    with open(mapped_between_prot_output, 'w') as file1:
        with open(between_prot_input, 'r') as f:
            for l in f:
                line = l.split()
                if len(line)>20:
                    take = line[16]+'\t'+line[18]+'\t'+line[15]+'\t'+line[17]+'\t'+line[19]+'\t'+line[21]+'\t'+line[4]+'\t'+line[8]
                    take = take.split()
                    with open(mut_prot_input, 'r') as mu:
                        for m in mu:
                            m = m.split()
                            if m[0] == take[0] and m[1]==take[4]:
                                take2 = take[0]+'\t'+take[2]+'\t'+take[4]+'\t'+take[6]+'\t'+'PTMcode'
                                fi = take2.split()
                                with open(yeastID_input) as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            with open(summary_output, 'a+') as summary:
                                                summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'BetweenProteins'+'\t'+'PTMcode'+'\n')
                                            if take2 > str(0):
                                                with open(mapped_between_prot_output, 'a') as file1:
                                                    file1.write(take2+'\n')
                                                    continue
                            if m[0] == take[1] and m[1] == take[5]:
                                take3 = take[1]+'\t'+take[3]+'\t'+take[5]+'\t'+take[7]+'\t'+'PTMcode'
                                fi=take3.split()
                                with open(yeastID_input) as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            with open(summary_output, 'a+') as summary:
                                                summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'BetweenProteins'+'\t'+'PTMcode'+'\n')
                                            if take3 > str(0):
                                                with open(mapped_between_prot_output, 'a+') as file1:
                                                    file1.write(take3+'\n')


# TODO: yeastID.txt is openen each line of mutation, suggestion to open it once store the di[0] and d1[2] in a
# dictionary di[2] as keys list of di[0]'s as  values, for a given fi[1] you can then retreive all corresponging di[0]'s
# don't open hotspot again and again it is already
# there are many other cases where something similar hapens
def hotspot(yeastID_input, regulatory_hotspots_input, mut_prot_input, mapped_hotspot_output, summary_output):

    """ PTMs containing motifs in a close proximity are named hotspots (Beltrao et al. Cell 2012)"""

    with open(mapped_hotspot_output, 'w') as hotspot:
        with open(regulatory_hotspots_input, 'r') as f:
            for l in f:
                line = l.split()
                with open(mut_prot_input, 'r') as mu:
                    for m in mu:
                        m = m.split()
                        if len(line) > 6:
                            if m[0] == line[2] and m[1] == line[3]:
                                take = line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[5]+'\t'+line[6]+'\t'+'PTMfunc'
                                fi = take.split()
                                with open(yeastID_input) as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[1]:
                                            with open(summary_output, 'a+') as summary:
                                                summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'HotSpot'+'\t'+'PTMFunc'+'\n')
                                            if take > str(0):
                                                with open(mapped_hotspot_output, 'a') as hotspot:
                                                    hotspot.write(take+'\n')
                          

def sum_file_map(summary_input, final_report_output):  

    """ reports all the results in a 'final-report' file """

    with open(final_report_output, 'w') as x:
        with open(summary_input) as fil1:
            for fi in OrderedDict.fromkeys(fil1):
                x.write(fi)


def resc(output_dir):
    """
    documentation needed....
    """

    try:
        r = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/3DID_aceksites_interfaceRes_sc.txt").read().decode()
        with open(interface_acet_file_path,'w') as h:
            h.write(r+'\n')
    except IOError:
        pass
    try:
        ri = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/3DID_phosphosites_interfaceRes_sc.txt").read().decode()
        with open(interface_phos_file_path,'w') as hi:
            hi.write(ri+'\n')
    except IOError:
        pass
    try:
        riu = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/3DID_ubisites_interfaceRessc_sc.txt").read().decode()
        with open(interface_ubiq_file_path,'w') as hiu:
            hiu.write(riu+'\n')
    except IOError:
        pass
    try:
        rac = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/SC_acet_interactions.txt").read().decode()
        with open(interact_acet_file_path,'w') as hia:
            hia.write(rac+'\n')
    except IOError:
        pass
    try:
        t = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/sc_btw_proteins.txt.zip").read()
        with open(between_prot_zip_file_path,'wb') as ht:
            ht.write(t)
    except IOError:
        pass
    try:
        zipfile.ZipFile(between_prot_zip_file_path, 'r').extractall()
    except IOError:
        pass
    try:
        rps = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/SC_psites_interactions_sc.txt").read().decode()
        with open(interact_phos_file_path,'w') as hip:
            hip.write(rps+'\n')
    except IOError:
        pass
    try:
        rui = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/SC_ubi_interactions_sc.txt").read().decode()
        with open(interact_ubiq_file_path,'w') as hui:
            hui.write(rui+'\n')
    except IOError:
        pass
    try:
        rin = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/sc_within_proteins.txt").read().decode()
        with open(within_prot_file_path,'w') as hin:
            hin.write(rin+'\n')
    except IOError:
        pass
    try:
        rsc = resource_stream("ymap", "/data/PTMcode+PTMfunc_data/schotspot_updated.txt").read().decode()
        with open(regulatory_hotspots_file_path,'w') as his:
            his.write(rsc+'\n')
    except IOError:
        pass
    return 


#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                               USEAGE       (Optional) 
#------------------------------------------------------------------------------------------------------------
#This usage strategy is optional, and a user can use above written codes in any convenient way as
#required by experiemental settings and data interpretation (see README for proper use)
#////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Create Ontology and Annotation objects for use in GO enrichment (see enrich() function)
ontology = go.Ontology()
annotations = go.Annotations('sgd')

wd = os.getcwd() 

# Setup filenames
# Input files
mutation_prot_file = 'mutation.txt'
mutation_gene_file = 'mutated_proteins.txt'

# Intermediate processing files
uniprot_file = 'uniprot_mod_raw.txt' # downloaded
bact_file = 'bact.txt'
ptms_file = 'PTMs.txt'
domains_file = 'domains.txt'
nucleotide_file = 'nucleotide.txt'
pdb_file = 'pdb.txt'

gff_file = 'gff.txt' # downloaded
frmt_file = 'frmt.txt'

yeastID_file = 'yeastID.txt' # downloaded
d_id_map_file = 'd_id_map.txt' # Links Uniprot ID to Gene ID, ORF start, ORF stop, strand orientation 
sites_id_file = 'sites_id.txt'
ptm_id_file = 'PTM_id_file.txt'
domain_id_file = 'id_domain.txt'
nucleotide_id_file = 'id_nucleotide.txt' # Nucleotide binding sites

# Downloaded and copied PTMfunc and PTMcode data files
interface_acet_file = '3DID_aceksites_interfaceRes_sc.txt'
interface_phos_file = '3DID_phosphosites_interfaceRes_sc.txt'
interface_ubiq_file = '3DID_ubisites_interfaceRessc_sc.txt'
regulatory_hotspots_file = 'schotspot_updated.txt'
interact_acet_file = 'SC_acet_interactions.txt'
interact_phos_file = 'SC_psites_interactions_sc.txt'
interact_ubiq_file = 'SC_ubi_interactions_sc.txt'
within_prot_file = 'sc_within_proteins.txt'
between_prot_zip_file = 'sc_btw_proteins.txt.zip'
between_prot_file = 'sc_btw_proteins.txt'

uniprot_biogrid_file = 'uniprot_bioGrid.txt' # downloaded

# Output files
summary_file = 'summary.txt'

mapped_sites_file = 'ab_mutation_file.txt'
mapped_ptms_file = 'mapped_ptms.txt'
mapped_domains_file = 'domains_mapped.txt'
mapped_nucleotide_file = 'nucleotide_map.txt'
mapped_mutation_pos_file = 'mutation_id.txt'
mapped_struct_file = 'stru_mutation.txt' # Structural regions of protein e.g. beta sheet, alpha helix, turn 

general_mapped_interface_file = 'interface_mutation.txt' # This variable refers to the name of the file written by ymap.interface().
# The mapped_interface variables are used in parts of this program for semantic clarity only
mapped_interface_acet_file = general_mapped_interface_file
mapped_interface_phos_file = general_mapped_interface_file
mapped_interface_ubiq_file = general_mapped_interface_file
general_mapped_ppi_file = 'ppi_mutation.txt'
# The mapped_interact variables are used in parts of this program for semantic clarity only
mapped_interact_acet_file = general_mapped_ppi_file
mapped_interact_phos_file = general_mapped_ppi_file
mapped_interact_ubiq_file = general_mapped_ppi_file
mapped_hotspot_file = 'hotspot.txt'
mapped_within_prot_file = 'within_protein.txt'
mapped_between_prot_file = 'ptm_between_proteins.txt'

p_value_file = 'pvalue.txt'
biog_file = 'biog.txt'
final_report_file = 'final_report.txt'

# Setup directory tree paths 
# Data, input, output directories
data_dir_path = os.path.join(wd, 'ymap_data')
input_dir_path = os.path.join(wd, 'input')
output_dir_path = os.path.join(wd, 'output') 

# Paths for output directories
domains_dir_path = os.path.join(output_dir_path, 'Domains') 
ptms_dir_path = os.path.join(output_dir_path, 'PTMs')
nuc_bind_dir_path = os.path.join(output_dir_path, 'Nucleotide_binding')
ab_sites_dir_path = os.path.join(output_dir_path, 'A-B-sites')
pdb_dir_path = os.path.join(output_dir_path, 'PDB')
interface_dir_path = os.path.join(output_dir_path, 'Interface')
interface_acet_dir_path = os.path.join(interface_dir_path, 'Acetylation')
interface_phos_dir_path = os.path.join(interface_dir_path, 'Phosphorylation')
interface_ubiq_dir_path = os.path.join(interface_dir_path, 'Ubiquitination')
ppi_dir_path = os.path.join(output_dir_path, 'PPI')
ppi_acet_dir_path = os.path.join(ppi_dir_path, 'Acetylation')
ppi_phos_dir_path = os.path.join(ppi_dir_path, 'Phosphorylation')
ppi_ubiq_dir_path = os.path.join(ppi_dir_path, 'Ubiquitination')
ptm_within_dir_path = os.path.join(output_dir_path, 'PTMs_within_Proteins')
ptm_between_dir_path = os.path.join(output_dir_path, 'PTMs_between_Proteins')
ptm_hotspot_dir_path = os.path.join(output_dir_path, 'PTMs_hotSpots')

dirs = [
    domains_dir_path,  
    ptms_dir_path, 
    nuc_bind_dir_path, 
    ab_sites_dir_path, 
    pdb_dir_path, 
    interface_dir_path, 
    interface_acet_dir_path, 
    interface_phos_dir_path, 
    interface_ubiq_dir_path, 
    ppi_dir_path, 
    ppi_acet_dir_path, 
    ppi_phos_dir_path, 
    ppi_ubiq_dir_path, 
    ptm_within_dir_path, 
    ptm_between_dir_path, 
    ptm_hotspot_dir_path
]

def makeDirTree(dirs):
    for directory in dirs:
        os.makedirs(directory)

# Setup paths to data, input and output files
# Input files
mutation_prot_file_path = os.path.join(input_dir_path, mutation_prot_file)
mutation_gene_file_path = os.path.join(input_dir_path, mutation_gene_file)

# Intermediate processing files
uniprot_file_path = os.path.join(data_dir_path, uniprot_file)
bact_file_path = os.path.join(data_dir_path, bact_file)
ptms_file_path = os.path.join(data_dir_path, ptms_file)
domains_file_path = os.path.join(data_dir_path, domains_file)
nucleotide_file_path = os.path.join(data_dir_path, nucleotide_file)
pdb_file_path = os.path.join(data_dir_path, pdb_file)

gff_file_path = os.path.join(data_dir_path, gff_file)
frmt_file_path = os.path.join(data_dir_path, frmt_file)

yeastID_file_path = os.path.join(data_dir_path, yeastID_file)
d_id_map_file_path = os.path.join(data_dir_path, d_id_map_file) 
sites_id_file_path = os.path.join(data_dir_path, sites_id_file)
ptm_id_file_path = os.path.join(data_dir_path, ptm_id_file)
domain_id_file_path = os.path.join(data_dir_path, domain_id_file)
nucleotide_id_file_path = os.path.join(data_dir_path, nucleotide_id_file)

# Downloaded and copied PTMfunc and PTMcode data files
interface_acet_file_path = os.path.join(data_dir_path, interface_acet_file)
interface_phos_file_path = os.path.join(data_dir_path, interface_phos_file)
interface_ubiq_file_path = os.path.join(data_dir_path, interface_ubiq_file)
regulatory_hotspots_file_path = os.path.join(data_dir_path, regulatory_hotspots_file)
interact_acet_file_path = os.path.join(data_dir_path, interact_acet_file)
interact_phos_file_path = os.path.join(data_dir_path, interact_phos_file)
interact_ubiq_file_path = os.path.join(data_dir_path, interact_ubiq_file)
within_prot_file_path = os.path.join(data_dir_path, within_prot_file)
between_prot_zip_file_path = os.path.join(data_dir_path, between_prot_zip_file)
between_prot_file_path = os.path.join(data_dir_path, between_prot_file)

uniprot_biogrid_file_path = os.path.join(data_dir_path, uniprot_biogrid_file)

# Paths to output files
summary_file_path = os.path.join(output_dir_path, summary_file)
final_report_file_path = os.path.join(output_dir_path, final_report_file)

mapped_sites_file_path = os.path.join(output_dir_path, mapped_sites_file)
mapped_ptms_file_path = os.path.join(output_dir_path, mapped_ptms_file)
mapped_domains_file_path = os.path.join(output_dir_path, mapped_domains_file)
mapped_nucleotide_file_path = os.path.join(output_dir_path, mapped_nucleotide_file)
mapped_mutation_pos_file_path = os.path.join(output_dir_path, mapped_mutation_pos_file)
mapped_struct_file_path = os.path.join(output_dir_path, mapped_struct_file)

mapped_interface_acet_file_path = os.path.join(output_dir_path, mapped_interface_acet_file)
mapped_interface_phos_file_path = os.path.join(output_dir_path, mapped_interface_phos_file)
mapped_interface_ubiq_file_path = os.path.join(output_dir_path, mapped_interface_ubiq_file)
mapped_interact_acet_file_path = os.path.join(output_dir_path, mapped_interact_acet_file)
mapped_interact_phos_file_path = os.path.join(output_dir_path, mapped_interact_phos_file)
mapped_interact_ubiq_file_path = os.path.join(output_dir_path, mapped_interact_ubiq_file)
mapped_hotspot_file_path = os.path.join(output_dir_path, mapped_hotspot_file)
mapped_within_prot_file_path = os.path.join(output_dir_path, mapped_within_prot_file)
mapped_between_prot_file_path = os.path.join(output_dir_path, mapped_between_prot_file)


def data(): 

    """ this function will download and clean required data to run ymap methods smoothly """

    start_time = time.time()
    os.mkdir(data_dir_path)
    os.chdir(data_dir_path)
    try:
        resc(data_dir_path)
    except IOError:
        pass
    try:
        pTMdata(uniprot_file_path)
    except IOError:
        pass
    try:
        make_ptms_file(uniprot_file_path, ptms_file_path)
    except IOError:
        pass
    try:
        iD(yeastID_file_path)  # this method does not return anything, all the rest below also don't
    except IOError:
        pass
    try:
        pmap(yeastID_file_path, ptms_file_path, ptm_id_file_path)
    except IOError:
        pass
    try:
        make_domains_file(uniprot_file_path, domains_file_path)
    except IOError:
        pass
    try:
        d_map(yeastID_file_path, domains_file_path, domain_id_file_path)
    except IOError:
        pass
    try:
        make_bact_file(uniprot_file_path, bact_file_path)
    except IOError:
            pass
    try:
        id(yeastID_file_path, bact_file_path, sites_id_file_path)
    except IOError:
            pass
    try:
        bioGrid(uniprot_biogrid_file_path)
    except IOError:
            pass
    try:
        make_pdb_file(uniprot_file_path, pdb_file_path)
    except IOError:
        pass
    try:
        gff(gff_file_path)
    except IOError:
        pass
    try:
        frmt(gff_file_path, frmt_file_path)
    except IOError:
        pass
    try:
        id_map(yeastID_file_path, frmt_file_path, d_id_map_file_path)
    except IOError:
        pass
    try:
        make_nucleotide_file(uniprot_file_path, nucleotide_file_path)
    except IOError:
        pass
    try:
        n_map(yeastID_file_path, nucleotide_file_path, nucleotide_id_file_path)
    except IOError:
        pass
    try:
        z = zipfile.ZipFile(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+between_prot_zip_file, 'r')
        z.extractall()
    except IOError:
        pass
    try:
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+interface_acet_file, wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+interface_phos_file, wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+interface_ubiq_file, wd) 
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+interact_acet_file, wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+interact_phos_file, wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+interact_ubiq_file, wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+within_prot_file, wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+regulatory_hotspots_file, wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+between_prot_file, wd)    
    except IOError:
        pass
    os.chdir(wd)
    return "All required data downloaded in %s seconds" % (time.time() - start_time)


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#//////////////////////////////// Following two codes are used for return the mutations at proteins level \\\\\\\\\\\\\\\\\\
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def mutation_types_file(): 

    """ mutation type and amino acid change calculation where ref. and mutant base known """

    start_time = time.time()
    try:
        mutation_file(mutation_gene_file_path, gff_file_path, d_id_map_file_path, mutation_prot_file_path)
    except IOError:
        pass
    return "Mutations with mutations types are available to map on functional entities"


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#////////////////////////////////// Following series of codes will return three files - mapped mutations, pvalue and biog.txt - for each type of data types \\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


def ptm():

    """ PTMs mapping to mutations """

    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            a = ptm_map(mutation_prot_file_path, ptm_id_file_path, mapped_ptms_file_path, summary_file_path)
        except IOError:
            pass
        try:    
            p = enrich(mapped_ptms_file_path)
        except IOError:
            pass
        try:
            preWeb(uniprot_biogrid_file_path, mapped_ptms_file_path)
        except IOError:
            pass
        try:
            os.mkdir(ptms_dir_path)
            shutil.move(output_dir_path+"/"+mapped_ptms_file, ptms_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, ptms_dir_path)
            shutil.move(output_dir_path+"/"+biog_file, ptms_dir_path)
        except IOError:
            pass
    return "PTMs mapped in %s seconds" % (time.time() - start_time)


def domain():

    """ protein domain mapping """

    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            dom = dmap(mutation_prot_file_path, domain_id_file_path, mapped_domains_file_path, summary_file_path)
        except IOError:
            pass
        try:
           p = enrich(mapped_domains_file_path)  
        except IOError:
            pass
        try:
            preWeb(uniprot_biogrid_file_path, mapped_domains_file_path)
        except IOError:
            pass
        try:
            os.mkdir(domains_dir_path)
        except IOError:
            pass
        try:
            shutil.move(output_dir_path+"/"+mapped_domains_file, domains_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, domains_dir_path)
            shutil.move(output_dir_path+"/"+biog_file, domains_dir_path)
        except IOError:
            pass
    return "Domains mapped in %s seconds" % (time.time() - start_time)
    

def nucleo():

    """ DNA-protein binding motif mapping """

    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            nucleotide_map(mutation_prot_file_path, nucleotide_id_file_path, mapped_nucleotide_file_path, summary_file_path)
        except IOError:
            pass
        try:
           p = enrich(mapped_nucleotide_file_path)  
        except IOError:
            pass
        try:
            preWeb(uniprot_biogrid_file_path, mapped_nucleotide_file_path)
        except IOError:
            pass
        try:
            os.mkdir(nuc_bind_dir_path)
        except IOError:
            pass
        try:
            shutil.move(output_dir_path+"/"+mapped_nucleotide_file, nuc_bind_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, nuc_bind_dir_path)
            shutil.move(output_dir_path+"/"+biog_file, nuc_bind_dir_path)
        except IOError:
            pass
    return "Nucleotide_binding domains mapped in %s seconds" % (time.time() - start_time)


def ab():

    """ active and binding site mapping """

    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            mm = mmap(mutation_prot_file_path, sites_id_file_path, mapped_sites_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_sites_file_path)
        except IOError:
            pass
        try:
            preWeb(uniprot_biogrid_file_path, mapped_sites_file_path)
        except IOError:
            pass
        try:
            os.mkdir(ab_sites_dir_path)
            shutil.move(output_dir_path+"/"+mapped_sites_file, ab_sites_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, ab_sites_dir_path)
            shutil.move(output_dir_path+"/"+biog_file, ab_sites_dir_path)
        except IOError:
            pass
    return "Active-Binding proteins sites mapped in %s seconds" % (time.time() - start_time)


def struc_map():

    """ structural regions mapping """

    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            mu_map(yeastID_file_path, mutation_prot_file_path, mapped_mutation_pos_file_path)
        except IOError:
            pass
        try:
            pd = pdb(pdb_file_path, mapped_mutation_pos_file_path, mapped_struct_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_struct_file_path)
        except IOError:
            pass
        try:
            preWeb(uniprot_biogrid_file_path, mapped_struct_file_path)
        except IOError:
            pass
        try:
            os.mkdir(pdb_dir_path)
            shutil.move(output_dir_path+"/"+mapped_struct_file, pdb_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, pdb_dir_path)
            shutil.move(output_dir_path+"/"+biog_file, pdb_dir_path)
        except IOError:
            pass
        return "Mutations are mapped to structural features in %s seconds" % (time.time() - start_time)

def intf():

    """ east = effective data which shows PTMs present at interface, ppi and 
        domain (hotspot) this analaysis could lead to an effective way to interpret
        user's mutational data on Yeast proteins from PTMfunc (also 3DID db) and PTMcode 2.0"""
    
    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            interface(yeastID_file_path, interface_acet_file_path, mutation_prot_file_path, mapped_interface_acet_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_interface_acet_file_path)
        except IOError:
            pass 
        try:
            os.mkdir(interface_dir_path)
            os.mkdir(interface_acet_dir_path)
            shutil.move(output_dir_path+"/"+mapped_interface_acet_file, interface_acet_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, interface_acet_dir_path)
        except IOError:
            pass
        try:
            interface(yeastID_file_path, interface_phos_file_path, mutation_prot_file_path, mapped_interface_phos_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_interface_phos_file_path)
        except IOError:
            pass
        try:
            os.mkdir(interface_phos_dir_path)
            shutil.move(output_dir_path+"/"+mapped_interface_phos_file, interface_phos_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, interface_phos_dir_path)
        except IOError:
            pass
        try:   
            interface(yeastID_file_path, interface_ubiq_file_path, mutation_prot_file_path, mapped_interface_ubiq_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_interface_ubiq_file_path)
        except IOError:
            pass
        try:
            os.mkdir(interface_ubiq_dir_path)
            shutil.move(output_dir_path+"/"+mapped_interface_ubiq_file, interface_ubiq_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, interface_ubiq_dir_path)
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)

def pi():
    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            ppi(yeastID_file_path, interact_acet_file_path, mutation_prot_file_path, mapped_interact_acet_file_path, summary_file_path)
        except IOError:
            pass
        try:
            enrich(mapped_interact_acet_file_path)
        except IOError:
            pass
        try:
            os.mkdir(ppi_dir_path)
            os.mkdir(ppi_acet_dir_path)
            shutil.move(output_dir_path+"/"+mapped_interact_acet_file, ppi_acet_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, ppi_acet_dir_path)
        except IOError:
            pass
        try:
            ppi(yeastID_file_path, interact_phos_file_path, mutation_prot_file_path, mapped_interact_phos_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_interact_phos_file_path)
        except IOError:
            pass
        try:
            os.mkdir(ppi_phos_dir_path)
            shutil.move(output_dir_path+"/"+mapped_interact_phos_file, ppi_phos_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, ppi_phos_dir_path)
        except IOError:
            pass
    
        try:
            ppi(yeastID_file_path, interact_ubiq_file_path, mutation_prot_file_path, mapped_interact_ubiq_file_path, summary_file_path)
        except IOError:
            pass
        try:
            enrich(mapped_interact_ubiq_file_path)
        except IOError:
            pass 
        try:
            os.mkdir(ppi_ubiq_dir_path)
            shutil.move(output_dir_path+"/"+mapped_interact_ubiq_file, ppi_ubiq_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, ppi_ubiq_dir_path)
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)

def withP():
    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            withinPro(yeastID_file_path, within_prot_file_path, mutation_prot_file_path, mapped_within_prot_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_within_prot_file_path)
        except IOError:
            pass
        try:
            os.mkdir(ptm_within_dir_path)
            shutil.move(output_dir_path+"/"+mapped_within_prot_file, ptm_within_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, ptm_within_dir_path)
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)

def betweenP():
    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            betweenPro(yeastID_file_path, between_prot_file_path, mutation_prot_file_path, mapped_between_prot_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_between_prot_file_path)
        except IOError:
            pass
        try:
            os.mkdir(ptm_between_dir_path)
            shutil.move(output_dir_path+"/"+mapped_between_prot_file, ptm_between_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, ptm_between_dir_path)
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)


def hotS():
    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            hotspot(yeastID_file_path, regulatory_hotspots_file_path, mutation_prot_file_path, mapped_hotspot_file_path, summary_file_path)
        except IOError:
            pass
        try:
            p = enrich(mapped_hotspot_file_path)
        except IOError:
            pass
        try:
            os.mkdir(ptm_hotspot_dir_path)
            shutil.move(output_dir_path+"/"+mapped_hotspot_file, ptm_hotspot_dir_path)
            shutil.move(output_dir_path+"/"+p_value_file, ptm_hotspot_dir_path)
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ Following two codes with perform all the codes on all the data /////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def uniprot_data(): 

    """ to perform all functions on UniProt(like ptm, domain and ab () functions) all together """ 

    try:
        ptm()
    except IOError:
        pass
    try:
        domain()
    except IOError:
        pass
    try:
        ab()
    except IOError:
        pass
    try:
        struc_map()
    except IOError:
        pass
    try:
        nucleo()
    except IOError:
        pass
    return "The Uniprot data is resolved into functional for interpretation"


def functional_data():

    """ to perform all functions on UniProt(like ptm, domain and ab () functions) all together """

    if not os.path.exists(mutation_prot_file_path):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            intf()
        except IOError:
            pass
        try:
            pi()
        except IOError:
            pass
        try:
            withP()
        except IOError:
            pass
        try:
            betweenP()
        except IOError:
            pass
        try:
            hotS()
        except IOError:
            pass
        return "The data from PTMcode and PTMfunc on PTMs functional biasedness is resolved into functional for interpretation"


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////  Final module of ymap package for executing all the modules \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


def ymap_genes():

    """ returns all the results of all the codes of yMap; starting from genetics coordinates of proteins """

    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
        if not os.path.exists(output_dir_path):
            os.mkdir(output_dir_path)
            os.chdir(output_dir_path)
        else:
            os.chdir(output_dir_path)
        try:
            mutation_types_file()
        except IOError:
            pass
        try:
            uniprot_data()
        except IOError:
            pass
        try:
            functional_data()
        except IOError:
            pass
        try:
            sum_file_map(summary_file_path, final_report_file_path)
        except IOError:
            pass
        try:
            y = (time.time() - start_time)
            os.makedirs('yMap-results'+str(y))
        except IOError:
            pass
        try:    
            p = enrich(final_report_file_path)
        except IOError:
            pass
        try:
            preWeb(uniprot_biogrid_file_path, final_report_file_path)
        except IOError:
            pass
        try:
            shutil.move(output_dir_path+"/"+'PTMs', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move('Domains', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move('Nucleotide_binding',output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'A-B-sites', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PDB', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'Interface', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PPI', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PTMs_within_Proteins', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PTMs_between_Proteins',output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PTMs_hotSpots',output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(mutation_prot_file_path, output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(mutation_gene_file_path, output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+final_report_file, output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+p_value_file, output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+biog_file, output_dir_path+"/"+'yMap-results'+str(y))
            os.remove(mapped_mutation_pos_file_path)          
            os.remove(summary_file_path)
        except IOError: 
            pass
        os.chdir(wd)
        return "All functional data from genomic coordinates is ready in about %s seconds" % (time.time() - start_time)


def ymap_proteins():

    """ returns all the results of all the codes of yMap; starting from proteins level mutation positions """

    start_time = time.time()
    if not os.path.exists(mutation_prot_file_path):
        raise StopIteration('because of missing mutation file')
    else:
        os.makedirs(output_dir_path, exist_ok = True)
        os.chdir(output_dir_path)
        try:
            uniprot_data()
        except IOError:
            pass
        try:
            functional_data()
        except IOError:
            pass
        try:
            sum_file_map(summary_file_path, final_report_file_path)
        except IOError:
            pass
        try:
            y = (time.time() - start_time)
            os.makedirs('yMap-results'+str(y))
        except IOError:
            pass
        try:    
            p = enrich(final_report_file_path)
        except IOError:
            pass
        try:
            preWeb(uniprot_biogrid_file_path, final_report_file_path)
        except IOError:
            pass
        try:
            shutil.move(output_dir_path+"/"+'PTMs', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move('Domains', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move('Nucleotide_binding',output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'A-B-sites', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PDB', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'Interface', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PPI', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PTMs_within_Proteins', output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PTMs_between_Proteins',output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+'PTMs_hotSpots',output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(mutation_prot_file_path, output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+final_report_file, output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+p_value_file, output_dir_path+"/"+'yMap-results'+str(y))
            shutil.move(output_dir_path+"/"+biog_file, output_dir_path+"/"+'yMap-results'+str(y))
            os.remove(mapped_mutation_pos_file_path)
            os.remove(summary_file_path)
        except IOError: 
            pass
        os.chdir(wd)
        return "All functional data from proteins mutation-positions is ready in about %s seconds" % (time.time() - start_time)


def web(): 
    """ NOTE: to use the following function change to dir to respective folder to run web based analysis """   
    os.chdir(input('specify biog.txt path:'))   # specify biog.txt path:/yMap-results78.50792193412781
    bweb(biog_file)
    return "Web is ready for networks exploration of mutated proteins"

def path():
    "Path to the BioGrid ids path for visualisation"
    try:
        os.chdir(raw_input("paste here path to biog.txt file:"))
    except IOError:
        pass
    return "you need to provide path/to/biog.txt"

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-d', '--ydata', help='downloads required data to run yMap successfully')
    parser.add_argument('-g', '--ygenes', help='performs the yMap on genes level mutation positions')
    parser.add_argument('-p', '--yproteins', help='performs the yMap on proteins level mutation positions')
    parser.add_argument('-w', '--yweb', help='generates BioGrid web pages for interactome visualisation; '
                                             'paste the path to biog.txt file')

    args = parser.parse_args()
    if args.ydata:
        try:
            data()
        except IOError:
            pass
    elif args.ygenes:
        try:
            ymap_genes()
        except IOError:
            pass
    elif args.yproteins:
        try:
            ymap_proteins()
        except IOError:
            pass
    elif args.yweb:
        try:
            web()
        except IOError:
            pass
    else:
        print ("to run a function seek help")
