from __future__ import print_function, absolute_import, division, unicode_literals
import unittest
from tempfile import mkdtemp
from shutil import copy, rmtree
from os import chdir, mkdir, makedirs, path
from time import clock
from filecmp import cmp
import ymap.ymap as ymap

# Set this to True when you want to regenerate reference files for future testing
generate_ref_files = False

# Reference directory paths
ref_dir = path.abspath('test_files')

# Fake files - to test parsing files into data structures 
fake_yeastID = 'fake_yeastID.txt'
ref_fake_yeastID = path.join(ref_dir, fake_yeastID)
fake_uniprot_biogrid = 'fake_uniprot_biogrid.txt'
ref_fake_uniprot_biogrid = path.join(ref_dir, fake_uniprot_biogrid)
fake_gff = 'fake_gff.txt'
ref_fake_gff = path.join(ref_dir, fake_gff)
fake_d_id_map = 'fake_d_id_map.txt'
ref_fake_d_id_map = path.join(ref_dir, fake_d_id_map)
fake_mutation_gene = 'fake_mutated_genes.txt'
ref_fake_mutation_gene = path.join(ref_dir, fake_mutation_gene)
short_d_id_map = 'short_d_id_map.txt'
ref_short_d_id_map = path.join(ref_dir, short_d_id_map)
short_frmt = 'short_frmt.txt'
ref_short_frmt = path.join(ref_dir, short_frmt)

# Input file paths
ref_mutation_prot = path.join(ref_dir, ymap.mutation_prot_file)
ref_mutation_prot_converted = path.join(ref_dir, 'converted_' + ymap.mutation_prot_file)
ref_mutation_gene = path.join(ref_dir, ymap.mutation_gene_file)

# Downloaded files from UniProt and SGD (Saccharomyces Genome Database) and copied PTMfunc and PTMcode data files
ref_uniprot = path.join(ref_dir, ymap.uniprot_file)
ref_gff = path.join(ref_dir, ymap.gff_file)
ref_yeastID = path.join(ref_dir, ymap.yeastID_file)

ref_interface_acet = path.join(ref_dir, ymap.interface_acet_file)
ref_interface_phos = path.join(ref_dir, ymap.interface_phos_file)
ref_interface_ubiq = path.join(ref_dir, ymap.interface_ubiq_file)
ref_hotspot = path.join(ref_dir, ymap.regulatory_hotspots_file)
ref_interact_acet = path.join(ref_dir, ymap.interact_acet_file)
ref_interact_phos = path.join(ref_dir, ymap.interact_phos_file)
ref_interact_ubiq = path.join(ref_dir, ymap.interact_ubiq_file)
ref_within_prot = path.join(ref_dir, ymap.within_prot_file)
ref_between_prot = path.join(ref_dir, ymap.between_prot_file)

ref_uniprot_biogrid = path.join(ref_dir, ymap.uniprot_biogrid_file)

# Intermediate processing files
ref_bact = path.join(ref_dir, ymap.bact_file)
ref_ptms = path.join(ref_dir, ymap.ptms_file)
ref_domains = path.join(ref_dir, ymap.domains_file)
ref_nucleotide = path.join(ref_dir, ymap.nucleotide_file)
ref_pdb = path.join(ref_dir, ymap.pdb_file)

ref_frmt = path.join(ref_dir, ymap.frmt_file)

ref_d_id_map = path.join(ref_dir, ymap.d_id_map_file)
ref_sites_id = path.join(ref_dir, ymap.sites_id_file)
ref_ptm_id = path.join(ref_dir, ymap.ptm_id_file)
ref_domain_id = path.join(ref_dir, ymap.domain_id_file)
ref_nucleotide_id = path.join(ref_dir, ymap.nucleotide_id_file)
ref_struct_id = path.join(ref_dir, ymap.struct_id_file)

ref_interface_acet_id = path.join(ref_dir, ymap.interface_acet_id_file)
ref_interface_phos_id = path.join(ref_dir, ymap.interface_phos_id_file)
ref_interface_ubiq_id = path.join(ref_dir, ymap.interface_ubiq_id_file)
ref_interact_acet_id = path.join(ref_dir, ymap.interact_acet_id_file)
ref_interact_phos_id = path.join(ref_dir, ymap.interact_phos_id_file)
ref_interact_ubiq_id = path.join(ref_dir, ymap.interact_ubiq_id_file)
ref_hotspot_id = path.join(ref_dir, ymap.regulatory_hotspots_id_file)
ref_within_prot_id = path.join(ref_dir, ymap.within_prot_id_file)
ref_between_prot_id = path.join(ref_dir, ymap.between_prot_id_file)

# Output files
ref_mapped_sites = path.join(ref_dir, ymap.mapped_sites_file)
ref_mapped_ptms = path.join(ref_dir, ymap.mapped_ptms_file)
ref_mapped_domains = path.join(ref_dir, ymap.mapped_domains_file)
ref_mapped_nucleotide = path.join(ref_dir, ymap.mapped_nucleotide_file)
ref_mapped_struct = path.join(ref_dir, ymap.mapped_struct_file)

ref_mapped_interface_acet = path.join(ref_dir, 'acet_' + ymap.mapped_interface_acet_file)
ref_mapped_interface_phos = path.join(ref_dir, 'phos_' + ymap.mapped_interface_phos_file)
ref_mapped_interface_ubiq = path.join(ref_dir, 'ubiq_' + ymap.mapped_interface_ubiq_file)
ref_mapped_interact_acet = path.join(ref_dir, 'acet_' + ymap.mapped_interact_acet_file)
ref_mapped_interact_phos = path.join(ref_dir, 'phos_' + ymap.mapped_interact_phos_file)
ref_mapped_interact_ubiq = path.join(ref_dir, 'ubiq_' + ymap.mapped_interact_ubiq_file)
ref_mapped_hotspot = path.join(ref_dir, ymap.mapped_hotspot_file)
ref_mapped_within_prot = path.join(ref_dir, ymap.mapped_within_prot_file)
ref_mapped_between_prot = path.join(ref_dir, ymap.mapped_between_prot_file)

ref_uniprot_summary = path.join(ref_dir, 'uniprot_' + ymap.summary_file)
ref_func_summary = path.join(ref_dir, 'func_' + ymap.summary_file)
ref_summary = path.join(ref_dir, ymap.summary_file)
ref_final_report = path.join(ref_dir, ymap.final_report_file)
ref_errors = path.join(ref_dir, ymap.errors_file)
ref_p_value = path.join(ref_dir, ymap.p_value_file)
ref_biog = path.join(ref_dir, ymap.biog_file)


class Timer:
    '''Timer (to be used as context manager) to test speed of code section.'''
    
    def __enter__(self):
        self.start = clock()
        return self

    def __exit__(self, *args):
        self.end = clock()
        self.interval = self.end - self.start


class DataDownloadTest(unittest.TestCase):
    '''Test ydata downloads.
    
    Note: These tests may fail because the newly downloaded files may include new/modified data
    that is not in the reference files. Hence file comparison is probably not appropriate.
    Indeed, is testing even necessary for these functions?
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        # Setup temporary file paths
        self.uniprot = path.join(self.test_dir, ymap.uniprot_file)
        self.yeastID = path.join(self.test_dir, ymap.yeastID_file)
        self.gff = path.join(self.test_dir, ymap.gff_file)
        self.uniprot_biogrid = path.join(self.test_dir, ymap.uniprot_biogrid_file)
        
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
        
    def test_gff(self):
        ymap.gff(self.gff)
        # Perform checks
        self.assertTrue(path.isfile(self.gff))
        if generate_ref_files:
            copy(self.gff, ref_dir)
        self.assertTrue(cmp(self.gff, ref_gff))
        
    def test_iD(self):
        ymap.iD(self.yeastID)
        # Perform checks
        self.assertTrue(path.isfile(self.yeastID))
        if generate_ref_files:
            copy(self.yeastID, ref_dir)
        self.assertTrue(cmp(self.yeastID, ref_yeastID))
    
    def test_pTMdata(self):
        ymap.pTMdata(self.uniprot)
        # Perform checks
        self.assertTrue(path.isfile(self.uniprot))
        if generate_ref_files:
            copy(self.uniprot, ref_dir)
        self.assertTrue(cmp(self.uniprot, ref_uniprot))
    
    def test_bioGrid(self):
        ymap.bioGrid(self.uniprot_biogrid)
        # Perform checks
        self.assertTrue(path.isfile(self.uniprot_biogrid))
        if generate_ref_files:
            copy(self.uniprot_biogrid, ref_dir)
        self.assertTrue(cmp(self.uniprot_biogrid, ref_uniprot_biogrid))
    
    def test_resc(self):
        # To test copying text and zip files - UNNECESSARY?
        pass
    
    def test_extractZips(self):
        pass


# ydata processing
class UniprotRawProcessingTest(unittest.TestCase):
         
    '''Test UniProt file processing.'''
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        # Copy reference files to temporary directory
        self.uniprot_file = copy(ref_uniprot, self.test_dir)
        # Setup temporary file paths
        self.bact = path.join(self.test_dir, ymap.bact_file)
        self.ptms = path.join(self.test_dir, ymap.ptms_file)
        self.domains = path.join(self.test_dir, ymap.domains_file)
        self.nucleotide = path.join(self.test_dir, ymap.nucleotide_file)
        self.pdb = path.join(self.test_dir, ymap.pdb_file)
        
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
    
    def test_make_bact_file(self):
        ymap.make_bact_file(self.uniprot_file, self.bact)
        # Perform checks
        self.assertTrue(path.isfile(self.bact))
        if generate_ref_files:
            copy(self.bact, ref_dir)
        self.assertTrue(cmp(self.bact, ref_bact))
        
    def test_make_ptms_file(self):
        ymap.make_ptms_file(self.uniprot_file, self.ptms)
        # Perform checks
        self.assertTrue(path.isfile(self.ptms))
        if generate_ref_files:
            copy(self.ptms, ref_dir)
        self.assertTrue(cmp(self.ptms, ref_ptms))
    
    def test_make_domains_file(self):
        ymap.make_domains_file(self.uniprot_file, self.domains)
        # Perform checks
        self.assertTrue(path.isfile(self.domains))
        if generate_ref_files:
            copy(self.domains, ref_dir)
        self.assertTrue(cmp(self.domains, ref_domains))
    
    def test_make_nucleotide_file(self):
        ymap.make_nucleotide_file(self.uniprot_file, self.nucleotide)
        # Perform checks
        self.assertTrue(path.isfile(self.nucleotide))
        if generate_ref_files:
            copy(self.nucleotide, ref_dir)
        self.assertTrue(cmp(self.nucleotide, ref_nucleotide))
    
    def test_make_pdb_file(self):
        ymap.make_pdb_file(self.uniprot_file, self.pdb)
        # Perform checks
        self.assertTrue(path.isfile(self.pdb))
        if generate_ref_files:
            copy(self.pdb, ref_dir)
        self.assertTrue(cmp(self.pdb, ref_pdb))

# gff.txt processing
class GffProcessingTest(unittest.TestCase):
       
    '''Test gff.txt processing.'''
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        # Copy reference files to temporary directory
        self.gff = copy(ref_gff, self.test_dir)
        # Setup temporary file paths
        self.frmt = path.join(self.test_dir, ymap.frmt_file)
    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
        
    def test_frmt(self):
        ymap.frmt(self.gff, self.frmt)
        # Perform checks
        self.assertTrue(path.isfile(self.frmt))
        if generate_ref_files:
            copy(self.frmt, ref_dir)
        self.assertTrue(cmp(self.frmt, ref_frmt))

# yeastID.txt processing
class YeastIDProcessingTest(unittest.TestCase):        
    '''Test yeastID.txt processing.'''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        # Copy reference files to temporary directory
        self.yeastID = copy(ref_yeastID, ymap.yeastID_file)
        self.fake_yeastID = copy(ref_fake_yeastID, fake_yeastID)
        self.frmt = copy(ref_frmt, ymap.frmt_file)
        self.short_frmt = copy(ref_short_frmt, ymap.frmt_file)
        self.bact = copy(ref_bact, ymap.bact_file)
        self.ptms = copy(ref_ptms, ymap.ptms_file)
        self.domains = copy(ref_domains, ymap.domains_file)
        self.nucleotide = copy(ref_nucleotide, ymap.nucleotide_file)
        self.struct = copy(ref_pdb, ymap.pdb_file)
        # Setup reference data structures
        self.gene_names = {
            'P25491': [('YDJ1','YNL064C')],
            'P08521': [('NA', 'YOR302W')],
            'P10081': [('TIF1', 'YKR059W'), ('TIF2', 'YJL138C')],
            'E9PAE3': [('truncated TYB', 'NA')],
            'E2QC18': [('NA', 'NA')]
            }
        self.gene_names_by_locus = {
            'YNL064C': ['P25491', 'YNL064C', 'YDJ1'],
            'YOR302W': ['P08521', 'YOR302W', 'NA'],
            'YKR059W': ['P10081', 'YKR059W', 'TIF1'],
            'YJL138C': ['P10081', 'YJL138C', 'TIF2'],
            'NA': ['E2QC18', 'NA', 'NA'] # Not really sure how to cope with this difficulty (genes that have locus id '' will overwrite items in dictionary)
            }
        self.gene_names_complete = ymap.parse_gene_names(self.yeastID) #TODO: I know this is bad practice in testing
        # For an example of how to avoid this bad bad practice see what has been done for test_id_map() method.
        # This solution involves making a load of custom files including only keys that you expect to find
        # in the above gene_names dictionary. This would be tedious and I don't have the patience.
        # Perhaps this whole testing file is bad practice? I mean, nice idea, but badly implemented.
        self.gene_names_by_locus_complete = ymap.parse_names_by_locus(self.gene_names_complete) # Again bad practice...
        # Setup temporary file paths                         
        self.d_id_map = path.join(self.test_dir, ymap.d_id_map_file)
        self.short_d_id_map = path.join(self.test_dir, short_d_id_map)
        self.sites_id = path.join(self.test_dir, ymap.sites_id_file)
        self.ptm_id = path.join(self.test_dir, ymap.ptm_id_file)
        self.domain_id = path.join(self.test_dir, ymap.domain_id_file)
        self.nucleotide_id = path.join(self.test_dir, ymap.nucleotide_id_file)
        self.struct_id = path.join(self.test_dir, ymap.struct_id_file)
        # Try the Timer
        self.timer = Timer()
        
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
    
    def test_parse_gene_names(self):
        gene_names = ymap.parse_gene_names(self.fake_yeastID)
        # Perform checks
        self.assertEqual(gene_names, self.gene_names)
        
    def test_parse_gene_names_by_locus(self):
        gene_names_by_locus = ymap.parse_names_by_locus(self.gene_names)
        # Perform checks
        self.assertEqual(gene_names_by_locus, self.gene_names_by_locus)
        
    def test_id_map(self):
        ymap.id_map(self.gene_names_by_locus, self.short_frmt, self.short_d_id_map)
        # Perform checks
        self.assertTrue(path.isfile(self.short_d_id_map))
        if generate_ref_files:
            copy(self.short_d_id_map, ref_dir)
        self.assertTrue(cmp(self.short_d_id_map, ref_short_d_id_map))
    
    def test_id(self):
        ymap.id(self.gene_names_complete, self.bact, self.sites_id)
        # Perform checks
        self.assertTrue(path.isfile(self.sites_id))
        if generate_ref_files:
            copy(self.sites_id, ref_dir)
        self.assertTrue(cmp(self.sites_id, ref_sites_id))
    
    def test_pmap(self):
        ymap.pmap(self.gene_names_complete, self.ptms, self.ptm_id)
        # Perform checks
        self.assertTrue(path.isfile(self.ptm_id))
        if generate_ref_files:
            copy(self.ptm_id, ref_dir)
        self.assertTrue(cmp(self.ptm_id, ref_ptm_id))
    
    def test_d_map(self):
        ymap.d_map(self.gene_names_complete, self.domains, self.domain_id)
        # Perform checks
        self.assertTrue(path.isfile(self.domain_id))
        if generate_ref_files:
            copy(self.domain_id, ref_dir)
        self.assertTrue(cmp(self.domain_id, ref_domain_id))
    
    def test_n_map(self):
        ymap.n_map(self.gene_names_complete, self.nucleotide, self.nucleotide_id)
        # Perform checks
        self.assertTrue(path.isfile(self.nucleotide_id))
        if generate_ref_files:
            copy(self.nucleotide_id, ref_dir)
        self.assertTrue(cmp(self.nucleotide_id, ref_nucleotide_id))
    
    def test_mu_map(self):
        ymap.mu_map(self.gene_names_complete, self.struct, self.struct_id)
        # Perform checks
        self.assertTrue(path.isfile(self.struct_id))
        if generate_ref_files:
            copy(self.struct_id, ref_dir)
        self.assertTrue(cmp(self.struct_id, ref_struct_id))
        
    def test_interface_map(self):
        # Copy reference files to test directory
        self.interface_acet = copy(ref_interface_acet, ymap.interface_acet_file)
        self.interface_phos = copy(ref_interface_phos, ymap.interface_phos_file) # Could implement, but maybe a test for interface_acet_id is sufficient?
        self.interface_ubiq = copy(ref_interface_ubiq, ymap.interface_ubiq_file)
        # Setup temporary files
        self.interface_acet_id = path.join(self.test_dir, ymap.interface_acet_id_file)
        self.interface_phos_id = path.join(self.test_dir, ymap.interface_phos_id_file)
        self.interface_ubiq_id = path.join(self.test_dir, ymap.interface_ubiq_id_file)
        # Run ymap function under test
        ymap.interface_map(self.gene_names_by_locus_complete, self.interface_acet, self.interface_acet_id)
        # Perform checks
        self.assertTrue(path.isfile(self.interface_acet_id))
        if generate_ref_files:
            copy(self.interface_acet_id, ref_dir)
        self.assertTrue(cmp(self.interface_acet_id, ref_interface_acet_id))
        
    def test_ppi_map(self):
        # Copy reference files to test directory
        self.interact_acet = copy(ref_interact_acet, ymap.interact_acet_file)
        self.interact_phos = copy(ref_interact_phos, ymap.interact_phos_file)
        self.interact_ubiq = copy(ref_interact_ubiq, ymap.interact_ubiq_file)
        # Setup temporary files
        self.interact_acet_id = path.join(self.test_dir, ymap.interact_acet_id_file)
        self.interact_phos_id = path.join(self.test_dir, ymap.interact_phos_id_file)
        self.interact_ubiq_id = path.join(self.test_dir, ymap.interact_ubiq_id_file)
        # Run ymap function under test
        ymap.ppi_map(self.gene_names_by_locus_complete, self.interact_acet, self.interact_acet_id)
        # Perform checks
        self.assertTrue(path.isfile(self.interact_acet_id))
        if generate_ref_files:
            copy(self.interact_acet_id, ref_dir)
        self.assertTrue(cmp(self.interact_acet_id, ref_interact_acet_id))
        
    def test_hotspot_map(self):
        # Copy reference files to test directory
        self.hotspot = copy(ref_hotspot, ymap.regulatory_hotspots_file)
        # Setup temporary files
        self.hotspot_id = path.join(self.test_dir, ymap.regulatory_hotspots_id_file)     
        # Run ymap function under test
        ymap.hotspot_map(self.gene_names_by_locus_complete, self.hotspot, self.hotspot_id)
        # Perform checks
        self.assertTrue(path.isfile(self.hotspot_id))
        if generate_ref_files:
            copy(self.hotspot_id, ref_dir)
        self.assertTrue(cmp(self.hotspot_id, ref_hotspot_id))
        
    def test_withinPro_map(self):
        # Copy reference files to test directory
        self.within_prot = copy(ref_within_prot, ymap.within_prot_file)
        # Setup temporary files
        self.within_prot_id = path.join(self.test_dir, ymap.within_prot_id_file)
        # Run ymap function under test
        ymap.withinPro_map(self.gene_names_by_locus_complete, self.within_prot, self.within_prot_id)
        # Perform checks
        self.assertTrue(path.isfile(self.within_prot_id))
        if generate_ref_files:
            copy(self.within_prot_id, ref_dir)
        self.assertTrue(cmp(self.within_prot_id, ref_within_prot_id))
        
    def test_betweenPro_map(self):
        # Copy reference files to test directory
        self.between_prot = copy(ref_between_prot, ymap.between_prot_file)
        # Setup temporary files
        self.between_prot_id = path.join(self.test_dir, ymap.between_prot_id_file)
        # Run ymap function under test
        ymap.betweenPro_map(self.gene_names_by_locus_complete, self.between_prot, self.between_prot_id)
        # Perform checks
        self.assertTrue(path.isfile(self.between_prot_id))
        if generate_ref_files:
            copy(self.between_prot_id, ref_dir)
        self.assertTrue(cmp(self.between_prot_id, ref_between_prot_id))
        
    
# yproteins
# uniprot_data()
class UniprotMappingTest(unittest.TestCase):
    '''Test mapping of mutations to UniProt annotations.'''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        # Copy reference files to temporary directory
        self.mutation_prot = copy(ref_mutation_prot, ymap.mutation_prot_file) # Input file
        self.yeastID = copy(ref_yeastID, ymap.yeastID_file)
        self.d_id_map = copy(ref_d_id_map, ymap.d_id_map_file)
        self.sites_id = copy(ref_sites_id, ymap.sites_id_file)
        self.ptm_id = copy(ref_ptm_id, ymap.ptm_id_file)
        self.domain_id = copy(ref_domain_id, ymap.domain_id_file)
        self.nucleotide_id = copy(ref_nucleotide_id, ymap.nucleotide_id_file)
        self.struct_id = copy(ref_struct_id, ymap.struct_id_file)
        # Setup reference data structures
        self.mutations = {
            'YDJ1': {'116'},
            'DPH5': {'172'},
            'IMD3': {'121'},
            'LRG1': {'155'},
            'MYO5': {'719'},
            'RRP5': {'1003'},
            'TKL1': {'335', '469'},
            'SKN7': {'412'},
            'NTH1': {'478'},
            'RPL4B': {'2'},
            'TSC11': {'84'},
            'HTB2': {'37', '46'},
            'ACT1': {'290', '60', '198'},
            'BDH1': {'62'},
            'YPK1': {'357'},
            'YEF3': {'467'},
            'HTA1': {'96'},
            'CDC28': {'254'},
            'CDC24': {'740'},
            'CKS1': {'168'},
            'FUS3': {'179'}
            }
        self.gene_names_complete = ymap.parse_gene_names(self.yeastID) # Again, bad practice because tests in this class may fail if ymap.parse_gene_names fails i.e. tests lose their isolation
        # Setup temporary file paths
        self.summary = path.join(self.test_dir, ymap.summary_file)         
        self.mapped_sites = path.join(self.test_dir, ymap.mapped_sites_file)
        self.mapped_ptms = path.join(self.test_dir, ymap.mapped_ptms_file)
        self.mapped_domains = path.join(self.test_dir, ymap.mapped_domains_file)
        self.mapped_nucleotide = path.join(self.test_dir, ymap.mapped_nucleotide_file)
        self.mapped_struct = path.join(self.test_dir, ymap.mapped_struct_file)
                    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
    
    def test_parse_mutations(self):
        mutations = ymap.parse_mutations(self.mutation_prot)
        # Perform checks
        self.assertEqual(mutations, self.mutations)
    
    def test_pdb(self):
        ymap.pdb(self.mutations, self.struct_id, self.mapped_struct, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_struct))
        if generate_ref_files:
            copy(self.mapped_struct, ref_dir)
        self.assertTrue(cmp(self.mapped_struct, ref_mapped_struct))

    def test_mmap(self):
        ymap.mmap(self.mutations, self.sites_id, self.mapped_sites, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_sites))
        if generate_ref_files:
            copy(self.mapped_sites, ref_dir)
        self.assertTrue(cmp(self.mapped_sites, ref_mapped_sites))
    
    def test_ptm_map(self):
        ymap.ptm_map(self.mutations, self.ptm_id, self.mapped_ptms, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_ptms))
        if generate_ref_files:
            copy(self.mapped_ptms, ref_dir)
        self.assertTrue(cmp(self.mapped_ptms, ref_mapped_ptms))
    
    def test_dmap(self):
        ymap.dmap(self.mutations, self.domain_id, self.mapped_domains, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_domains))
        if generate_ref_files:
            copy(self.mapped_domains, ref_dir)
        self.assertTrue(cmp(self.mapped_domains, ref_mapped_domains))
    
    def test_nucleotide_map(self):
        ymap.nucleotide_map(self.mutations, self.nucleotide_id, self.mapped_nucleotide, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_nucleotide))
        if generate_ref_files:
            copy(self.mapped_nucleotide, ref_dir)
        self.assertTrue(cmp(self.mapped_nucleotide, ref_mapped_nucleotide))
    
    def test_summary(self):
        ymap.pdb(self.mutations, self.struct_id, self.mapped_struct, self.summary)
        ymap.mmap(self.mutations, self.sites_id, self.mapped_sites, self.summary)
        ymap.ptm_map(self.mutations, self.ptm_id, self.mapped_ptms, self.summary)
        ymap.dmap(self.mutations, self.domain_id, self.mapped_domains, self.summary)
        ymap.nucleotide_map(self.mutations, self.nucleotide_id, self.mapped_nucleotide, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.summary))
        if generate_ref_files:
            copy(self.summary, ref_uniprot_summary)
        self.assertTrue(cmp(self.summary, ref_uniprot_summary))

        
# functional_data()
class FunctionalMappingTest(unittest.TestCase):
    '''Test mapping of mutations to PTMFunc and PTMcode annotations.'''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        # Copy reference files to temporary directory
        self.yeastID = copy(ref_yeastID, ymap.yeastID_file)        
        self.interface_acet_id = copy(ref_interface_acet_id, ymap.interface_acet_id_file)
#         self.interface_phos_id = copy(ref_interface_phos_id, ymap.interface_phos_id_file)
#         self.interface_ubiq_id = copy(ref_interface_ubiq_id, ymap.interface_ubiq_id_file)
        self.interact_acet_id = copy(ref_interact_acet_id, ymap.interact_acet_id_file)
#         self.interact_phos_id = copy(ref_interact_phos_id, ymap.interact_phos_id_file)
#         self.interact_ubiq_id = copy(ref_interact_ubiq_id, ymap.interact_ubiq_id_file)
        self.hotspot_id = copy(ref_hotspot_id, ymap.regulatory_hotspots_id_file)        
        self.within_prot_id = copy(ref_within_prot_id, ymap.within_prot_id_file)
        self.between_prot_id = copy(ref_between_prot_id, ymap.between_prot_id_file)
        # Setup reference data structures
        self.mutations = { # NOTE: test for parsing mutation file into dictionary is in UniprotMappingTest Class
            'YDJ1': {'116'},
            'DPH5': {'172'},
            'IMD3': {'121'},
            'LRG1': {'155'},
            'MYO5': {'719'},
            'RRP5': {'1003'},
            'TKL1': {'335', '469'},
            'SKN7': {'412'},
            'NTH1': {'478'},
            'RPL4B': {'2'},
            'TSC11': {'84'},
            'HTB2': {'37', '46'},
            'ACT1': {'290', '60', '198'},
            'BDH1': {'62'},
            'YPK1': {'357'},
            'YEF3': {'467'},
            'HTA1': {'96'},
            'CDC28': {'254'},
            'CDC24': {'740'},
            'CKS1': {'168'},
            'FUS3': {'179'}
            }
        # Setup temporary file paths
        self.summary = path.join(self.test_dir, ymap.summary_file)         
        self.mapped_interface_acet = path.join(self.test_dir, ymap.mapped_interface_acet_file)
        self.mapped_interface_phos = path.join(self.test_dir, ymap.mapped_interface_phos_file)
        self.mapped_interface_ubiq = path.join(self.test_dir, ymap.mapped_interface_ubiq_file)
        self.mapped_interact_acet = path.join(self.test_dir, ymap.mapped_interact_acet_file)
        self.mapped_interact_phos = path.join(self.test_dir, ymap.mapped_interact_phos_file)
        self.mapped_interact_ubiq = path.join(self.test_dir, ymap.mapped_interact_ubiq_file)
        self.mapped_hotspot = path.join(self.test_dir, ymap.mapped_hotspot_file)
        self.mapped_within_prot = path.join(self.test_dir, ymap.mapped_within_prot_file)
        self.mapped_between_prot = path.join(self.test_dir, ymap.mapped_between_prot_file)
                    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
        
    # intf()    
    def test_interface(self):
        # For each 3DID_<site>_interfaceRes_sc.txt file
        ymap.interface(self.mutations, self.interface_acet_id, self.mapped_interface_acet, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_interface_acet))
        if generate_ref_files:
            copy(self.mapped_interface_acet, ref_mapped_interface_acet)
        self.assertTrue(cmp(self.mapped_interface_acet, ref_mapped_interface_acet))
        
#         ymap.interface(self.mutations, self.interface_phos_id, self.mapped_interface_phos, self.summary)
#         
#         self.assertTrue(path.isfile(self.mapped_interface_phos))
#         if generate_ref_files:
#             copy(self.mapped_interface_phos, ref_mapped_interface_phos)
#         self.assertTrue(cmp(self.mapped_interface_acet, ref_mapped_interface_phos))
#         
#         ymap.interface(self.mutations, self.interface_ubiq_id, self.mapped_interface_ubiq, self.summary)
#         
#         self.assertTrue(path.isfile(self.mapped_interface_ubiq))
#         if generate_ref_files:
#             copy(self.mapped_interface_ubiq, ref_mapped_interface_ubiq)
#         self.assertTrue(cmp(self.mapped_interface_acet, ref_mapped_interface_ubiq))
    
    # pi()
    def test_ppi(self):
        # For each SC_<site>_interactions_symap.txt
        ymap.ppi(self.mutations, self.interact_acet_id, self.mapped_interact_acet, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_interact_acet))
        if generate_ref_files:
            copy(self.mapped_interact_acet, ref_mapped_interact_acet)
        self.assertTrue(cmp(self.mapped_interact_acet, ref_mapped_interact_acet))
        
#         ymap.ppi(self.mutations, self.interact_phos_id, self.mapped_interact_phos, self.summary)
#         
#         self.assertTrue(path.isfile(self.mapped_interact_phos))
#         if generate_ref_files:
#             copy(self.mapped_interact_phos, ref_mapped_interact_phos)
#         self.assertTrue(cmp(self.mapped_interact_acet, ref_mapped_interact_phos))
#         
#         ymap.ppi(self.mutations, self.interact_ubiq_id, self.mapped_interact_ubiq, self.summary)
#         
#         self.assertTrue(path.isfile(self.mapped_interact_ubiq))
#         if generate_ref_files:
#             copy(self.mapped_interact_ubiq, ref_mapped_interact_ubiq)
#         self.assertTrue(cmp(self.mapped_interact_acet, ref_mapped_interact_ubiq))
    
    # withP()    
    def test_withinPro(self):
        ymap.withinPro(self.mutations, self.within_prot_id, self.mapped_within_prot, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_within_prot))
        if generate_ref_files:
            copy(self.mapped_within_prot, ref_dir)
        self.assertTrue(cmp(self.mapped_within_prot, ref_mapped_within_prot))
    
    # betweenP()
    def test_betweenPro(self):
        ymap.betweenPro(self.mutations, self.between_prot_id, self.mapped_between_prot, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_between_prot))
        if generate_ref_files:
            copy(self.mapped_between_prot, ref_dir)
        self.assertTrue(cmp(self.mapped_between_prot, ref_mapped_between_prot))
    
    # hotS()
    def test_hotspot(self):
        ymap.hotspot(self.mutations, self.hotspot_id, self.mapped_hotspot, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.mapped_hotspot))
        if generate_ref_files:
            copy(self.mapped_hotspot, ref_dir)
        self.assertTrue(cmp(self.mapped_hotspot, ref_mapped_hotspot))
    
    def test_summary(self):
        ymap.interface(self.mutations, self.interface_acet_id, self.mapped_interface_acet, self.summary)
#         ymap.interface(self.mutations, self.interface_phos, self.mapped_interface_phos, self.summary)
#         ymap.interface(self.mutations, self.interface_ubiq, self.mapped_interface_ubiq, self.summary)
        ymap.ppi(self.mutations, self.interact_acet_id, self.mapped_interact_acet, self.summary)
#         ymap.ppi(self.mutations, self.interact_phos, self.mapped_interact_phos, self.summary)
#         ymap.ppi(self.mutations, self.interact_ubiq, self.mapped_interact_ubiq, self.summary)
        ymap.withinPro(self.mutations, self.within_prot_id, self.mapped_within_prot, self.summary)
        ymap.betweenPro(self.mutations, self.between_prot_id, self.mapped_between_prot, self.summary)
        ymap.hotspot(self.mutations, self.hotspot_id, self.mapped_hotspot, self.summary)
        # Perform checks
        self.assertTrue(path.isfile(self.summary))
        if generate_ref_files:
            copy(self.summary, ref_func_summary)
        self.assertTrue(cmp(self.summary, ref_func_summary))


# ygenes
class GeneToProteinTest(unittest.TestCase):
    '''Test conversion of gene-level mutation file to protein-level mutation file.'''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        # Copy reference files to temporary directory
        self.mutation_gene = copy(ref_mutation_gene, ymap.mutation_gene_file) # Input file
        self.fake_mutation_gene = copy(ref_fake_mutation_gene, fake_mutation_gene)
        self.gff = copy(ref_gff, ymap.gff_file)
        self.fake_gff = copy(ref_fake_gff, fake_gff)
        self.d_id_map = copy(ref_d_id_map, ymap.d_id_map_file)
        self.fake_d_id_map = copy(ref_fake_d_id_map, fake_d_id_map)
        # Setup reference data structures
        self.mutations = {
            'YPK1': {('chrXI', '205907', 'C', 'T')},
            'INTERGENIC': {('chrIV', '1334251', 'A', 'T')},
            'YMR134W': {('chrXIII', '537942', 'G', 'A')},
            'HEM3': {('chrIV', '93192', 'C', 'T'), ('chrIV', '93537', 'A', 'C')}
            }
        self.genome = {
            'chrI': '''CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT''',
            'chrII': '''AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGATGTTCAACCAAAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATAGGTTGTCTTTTTATCCCACTTCTTCGCACTTGTCTCTCGCTACTGCCGTGCAACAAACAC'''
            }
        self.gene_loci = {
            'NA': ['538', '792', '+', 'chrI'],
            'PAU8': ['1807', '2169', '-', 'chrI'],
            'BI2': ['36540', '38579', '+', 'chrmt']
            }
        # Setup temporary file paths
        self.mutation_prot = path.join(self.test_dir, ymap.mutation_prot_file)
        self.errors = path.join(self.test_dir, ymap.errors_file)     
                    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
        
    def test_mutation_file(self):
        ymap.mutation_file(self.mutation_gene, self.gff, self.d_id_map, self.mutation_prot, self.errors)
        # Perform checks
        self.assertTrue(path.isfile(self.mutation_prot))
        self.assertTrue(path.isfile(self.errors))
        if generate_ref_files:
            copy(self.mutation_prot, ref_mutation_prot_converted)
            copy(self.errors, ref_errors)
        self.assertTrue(cmp(self.mutation_prot, ref_mutation_prot_converted))
        self.assertTrue(cmp(self.errors, ref_errors))
    
    def test_revcomp(self):
        pass
    
    def test_translate_dna(self):
        pass
    
    def test_parse_gene_mutations(self):
        mutations = ymap.parse_gene_mutations(self.fake_mutation_gene)
        # Perform checks
        self.assertEqual(mutations, self.mutations)
        
    def test_parse_genome(self):
        genome = ymap.parse_genome(self.fake_gff)
        # Perform checks
        self.assertEqual(genome, self.genome)
    
    def test_parse_gene_loci(self):
        gene_loci = ymap.parse_gene_loci(self.fake_d_id_map)
        # Perform checks
        self.assertEqual(gene_loci, self.gene_loci)

    
# GO enrichment, Biogrid network and final report
class FinalAnalysisTest(unittest.TestCase):
    '''Test production of files containing p-values for GO enrichment and files containing Biogrid IDs.
    Test conversion of summary file into final report file.
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        # Copy reference files to temporary directory
        self.uniprot_biogrid = copy(ref_uniprot_biogrid, ymap.uniprot_biogrid_file)
        self.fake_uniprot_biogrid = copy(ref_fake_uniprot_biogrid, fake_uniprot_biogrid)
        # Setup reference data structures
        self.uniprot_biogrid_map = {
            'P25491': ['35759'],
            'P47036': [''],
            'P0CZ17': ['31428', '31431', '31433'],
            'P0CX80': ['36485', '36487']
            }
        self.uniprot_biogrid_map_complete = ymap.parse_biogrid(self.uniprot_biogrid) # Again bad practice
        # Setup temporary file paths
        self.final_report = path.join(self.test_dir, ymap.final_report_file)
        self.p_value = path.join(self.test_dir, ymap.p_value_file)
        self.biog = path.join(self.test_dir, ymap.biog_file)
                    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)

    def test_parse_biogrid(self):
        uniprot_biogrid_map = ymap.parse_biogrid(self.fake_uniprot_biogrid)
        # Perform checks
        self.assertEqual(uniprot_biogrid_map, self.uniprot_biogrid_map)
    
    def test_sum_file_map(self):
        self.summary = copy(ref_uniprot_summary, ymap.summary_file)
        ymap.sum_file_map(self.summary, self.final_report)
        # Perform checks
        self.assertTrue(path.isfile(self.final_report))
        if generate_ref_files:
            copy(self.final_report, ref_final_report)
        self.assertTrue(cmp(self.final_report, ref_final_report))
    
    def test_enrich(self):
        self.final_report = copy(ref_final_report, ymap.final_report_file)
        ymap.enrich(self.final_report)
        # Perform checks
        self.assertTrue(path.isfile(self.p_value))
        if generate_ref_files:
            copy(self.p_value, ref_p_value)
        self.assertTrue(cmp(self.p_value, ref_p_value))
    
    def test_preWeb(self):
        self.final_report = copy(ref_final_report, ymap.final_report_file)
        ymap.preWeb(self.uniprot_biogrid_map_complete, self.final_report)
        # Perform checks
        self.assertTrue(path.isfile(self.biog))
        if generate_ref_files:
            copy(self.biog, ref_biog)
        self.assertTrue(cmp(self.biog, ref_biog))
    
    
    # yweb
    def test_bweb(self):
        pass


# NOTE: Run from command line. cd to root ymap directory and run 'python -m unittest tests.tests' 
if __name__ == '__main__':
    unittest.main()
