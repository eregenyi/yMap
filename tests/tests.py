from __future__ import print_function, absolute_import, division, unicode_literals
import unittest
from tempfile import mkdtemp
from shutil import copy, rmtree
from os import chdir, mkdir, makedirs, path
from time import clock
from filecmp import cmp
import ymap.ymap as ymap

# Set this to True when you want to regenerate reference files for future testing
generate_ref_files = True

# Reference directory paths
ref_dir = path.abspath('test_files')

ref_mutation_prot = path.join(ref_dir, ymap.mutation_prot_file)
ref_mutation_prot_converted = path.join(ref_dir, 'converted_' + ymap.mutation_prot_file)
ref_mutation_gene = path.join(ref_dir, ymap.mutation_gene_file)


ref_uniprot_file = path.join(ref_dir, ymap.uniprot_file)
ref_bact = path.join(ref_dir, ymap.bact_file)
ref_ptms = path.join(ref_dir, ymap.ptms_file)
ref_domains = path.join(ref_dir, ymap.domains_file)
ref_nucleotide = path.join(ref_dir, ymap.nucleotide_file)
ref_pdb = path.join(ref_dir, ymap.pdb_file)

ref_gff = path.join(ref_dir, ymap.gff_file)
ref_frmt = path.join(ref_dir, ymap.frmt_file)

ref_yeastID = path.join(ref_dir, ymap.yeastID_file)
ref_d_id_map = path.join(ref_dir, ymap.d_id_map_file)
ref_sites_id = path.join(ref_dir, ymap.sites_id_file)
ref_ptm_id = path.join(ref_dir, ymap.ptm_id_file)
ref_domain_id = path.join(ref_dir, ymap.domain_id_file)
ref_nucleotide_id = path.join(ref_dir, ymap.nucleotide_id_file)

ref_uniprot_biogrid = path.join(ref_dir, ymap.uniprot_biogrid_file)


ref_uniprot_summary = path.join(ref_dir, 'uniprot_' + ymap.summary_file)
ref_func_summary = path.join(ref_dir, 'func_' + ymap.summary_file)
ref_summary = path.join(ref_dir, ymap.summary_file)
ref_final_report = path.join(ref_dir, ymap.final_report_file)

ref_mapped_sites = path.join(ref_dir, ymap.mapped_sites_file)
ref_mapped_ptms = path.join(ref_dir, ymap.mapped_ptms_file)
ref_mapped_domains = path.join(ref_dir, ymap.mapped_domains_file)
ref_mapped_nucleotide = path.join(ref_dir, ymap.mapped_nucleotide_file)
ref_mapped_mutation_pos = path.join(ref_dir, ymap.mapped_mutation_pos_file)
ref_mapped_struct = path.join(ref_dir, ymap.mapped_struct_file)

ref_interface_acet = path.join(ref_dir, ymap.interface_acet_file)
ref_interface_phos = path.join(ref_dir, ymap.interface_phos_file)
ref_interface_ubiq = path.join(ref_dir, ymap.interface_ubiq_file)
ref_hotspot = path.join(ref_dir, ymap.regulatory_hotspots_file)
ref_interact_acet = path.join(ref_dir, ymap.interact_acet_file)
ref_interact_phos = path.join(ref_dir, ymap.interact_phos_file)
ref_interact_ubiq = path.join(ref_dir, ymap.interact_ubiq_file)
ref_within_prot = path.join(ref_dir, ymap.within_prot_file)
ref_between_prot = path.join(ref_dir, ymap.between_prot_file)

ref_mapped_interface_acet = path.join(ref_dir, 'acet_' + ymap.mapped_interface_acet_file)
ref_mapped_interface_phos = path.join(ref_dir, 'phos_' + ymap.mapped_interface_phos_file)
ref_mapped_interface_ubiq = path.join(ref_dir, 'ubiq_' + ymap.mapped_interface_ubiq_file)
ref_mapped_interact_acet = path.join(ref_dir, 'acet_' + ymap.mapped_interact_acet_file)
ref_mapped_interact_phos = path.join(ref_dir, 'phos_' + ymap.mapped_interact_phos_file)
ref_mapped_interact_ubiq = path.join(ref_dir, 'ubiq_' + ymap.mapped_interact_ubiq_file)
ref_mapped_hotspot = path.join(ref_dir, ymap.mapped_hotspot_file)
ref_mapped_within_prot = path.join(ref_dir, ymap.mapped_within_prot_file)
ref_mapped_between_prot = path.join(ref_dir, ymap.mapped_between_prot_file)

ref_p_value = path.join(ref_dir, ymap.p_value_file)
ref_biog = path.join(ref_dir, ymap.biog_file)

# Create object to allow function calls
#TODO: Remove the class from ymap.py - because it is pointless
c = ymap.YGtPM()


class Timer:
    def __enter__(self):
        self.start = clock()
        return self

    def __exit__(self, *args):
        self.end = clock()
        self.interval = self.end - self.start
   
@unittest.skip("showing class skipping")
class DataDownloadTest(unittest.TestCase):
    '''
    Test ydata downloads - these tests are low priority
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
        
    def test_gff(self):
        pass
    
    def test_iD(self):
        # This test may fail because yeastID.txt is downloaded and data may have changed from that in test_files folder.
#         c.iD()
#         self.assertTrue(path.isfile('yeastID.txt'))
#         self.assertTrue(cmp('yeastID.txt', path.join(ref_dir, 'yeastID.txt')))
        pass
    
    def test_pTMdata(self):
        pass
    
    def test_bioGrid(self):
        pass
    
    def test_resc(self):
        # To test copying text and zip files - UNNECESSARY?
        pass

# ydata processing
@unittest.skip("showing class skipping")
class UniprotRawProcessingTest(unittest.TestCase):        
    '''
    Test uniprot_mod_raw.txt processing
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.uniprot_file = copy(ref_uniprot_file, self.test_dir)
        
        # Setup temporary file paths
        self.bact = path.join(self.test_dir, ymap.bact_file)
        self.ptms = path.join(self.test_dir, ymap.ptms_file)
        self.domains = path.join(self.test_dir, ymap.domains_file)
        self.nucleotide = path.join(self.test_dir, ymap.nucleotide_file)
        self.pdb = path.join(self.test_dir, ymap.pdb_file)
        
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
    
    def test_ab(self):
        c.ab(self.uniprot_file, self.bact)
        
        self.assertTrue(path.isfile(self.bact))
        if generate_ref_files:
            copy(self.bact, ref_dir)
        self.assertTrue(cmp(self.bact, ref_bact))
        
    def test_clean(self):
        c.clean(self.uniprot_file, self.ptms)
        
        self.assertTrue(path.isfile(self.ptms))
        if generate_ref_files:
            copy(self.ptms, ref_dir)
        self.assertTrue(cmp(self.ptms, ref_ptms))
    
    def test_dclean(self):
        c.dclean(self.uniprot_file, self.domains)
        
        self.assertTrue(path.isfile(self.domains))
        if generate_ref_files:
            copy(self.domains, ref_dir)
        self.assertTrue(cmp(self.domains, ref_domains))
    
    def test_nucleotide(self):
        c.nucleotide(self.uniprot_file, self.nucleotide)
        
        self.assertTrue(path.isfile(self.nucleotide))
        if generate_ref_files:
            copy(self.nucleotide, ref_dir)
        self.assertTrue(cmp(self.nucleotide, ref_nucleotide))
    
    def test_pdb_c(self):
        c.pdb_c(self.uniprot_file, self.pdb)
        
        self.assertTrue(path.isfile(self.pdb))
        if generate_ref_files:
            copy(self.pdb, ref_dir)
        self.assertTrue(cmp(self.pdb, ref_pdb))

    # gff.txt processing
@unittest.skip("showing class skipping")
class GffProcessingTest(unittest.TestCase):        
    '''
    Test gff.txt processing
    '''
    
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
        c.frmt(self.gff, self.frmt)
        
        self.assertTrue(path.isfile(self.frmt))
        if generate_ref_files:
            copy(self.frmt, ref_dir)
        self.assertTrue(cmp(self.frmt, ref_frmt))

    # yeastID.txt processing
@unittest.skip("showing class skipping")
class YeastIDProcessingTest(unittest.TestCase):        
    '''
    Test yeastID.txt processing
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.yeastID = copy(ref_yeastID, ymap.yeastID_file)
        self.frmt = copy(ref_frmt, ymap.frmt_file)
        self.bact = copy(ref_bact, ymap.bact_file)
        self.ptms = copy(ref_ptms, ymap.ptms_file)
        self.domains = copy(ref_domains, ymap.domains_file)
        self.nucleotide = copy(ref_nucleotide, ymap.nucleotide_file)
        
        # Setup temporary file paths                         
        self.d_id_map = path.join(self.test_dir, ymap.d_id_map_file)
        self.sites_id = path.join(self.test_dir, ymap.sites_id_file)
        self.ptm_id = path.join(self.test_dir, ymap.ptm_id_file)
        self.domain_id = path.join(self.test_dir, ymap.domain_id_file)
        self.nucleotide_id = path.join(self.test_dir, ymap.nucleotide_id_file)                        
    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
    
    def test_id_map(self):
        c.id_map(self.yeastID, self.frmt, self.d_id_map)
        
        self.assertTrue(path.isfile(self.d_id_map))
        if generate_ref_files:
            copy(self.d_id_map, ref_dir)
        self.assertTrue(cmp(self.d_id_map, ref_d_id_map))
    
    def test_id(self):
        c.id(self.yeastID, self.bact, self.sites_id)
        
        self.assertTrue(path.isfile(self.sites_id))
        if generate_ref_files:
            copy(self.sites_id, ref_dir)
        self.assertTrue(cmp(self.sites_id, ref_sites_id))
    
    def test_pmap(self):
        c.pmap(self.yeastID, self.ptms, self.ptm_id)
        
        self.assertTrue(path.isfile(self.ptm_id))
        if generate_ref_files:
            copy(self.ptm_id, ref_dir)
        self.assertTrue(cmp(self.ptm_id, ref_ptm_id))
    
    def test_d_map(self):
        c.d_map(self.yeastID, self.domains, self.domain_id)
        
        self.assertTrue(path.isfile(self.domain_id))
        if generate_ref_files:
            copy(self.domain_id, ref_dir)
        self.assertTrue(cmp(self.domain_id, ref_domain_id))
    
    def test_n_map(self):
        c.n_map(self.yeastID, self.nucleotide, self.nucleotide_id)
        
        self.assertTrue(path.isfile(self.nucleotide_id))
        if generate_ref_files:
            copy(self.nucleotide_id, ref_dir)
        self.assertTrue(cmp(self.nucleotide_id, ref_nucleotide_id))
    
    
    
    # yproteins
    # uniprot_data()
@unittest.skip("showing class skipping")
class UniprotMappingTest(unittest.TestCase):
    '''
    Test mapping of mutations to UniProt annotations
    '''
    
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
        self.pdb = copy(ref_pdb, ymap.pdb_file)
        
        # Setup temporary file paths
        self.summary = path.join(self.test_dir, ymap.summary_file)         
        self.mapped_sites = path.join(self.test_dir, ymap.mapped_sites_file)
        self.mapped_ptms = path.join(self.test_dir, ymap.mapped_ptms_file)
        self.mapped_domains = path.join(self.test_dir, ymap.mapped_domains_file)
        self.mapped_nucleotide = path.join(self.test_dir, ymap.mapped_nucleotide_file)
        self.mapped_mutation_pos = path.join(self.test_dir, ymap.mapped_mutation_pos_file)
        self.mapped_struct = path.join(self.test_dir, ymap.mapped_struct_file)
                    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)    
        
    def test_mu_map(self):
        c.mu_map(self.yeastID, self.mutation_prot, self.mapped_mutation_pos)
        
        self.assertTrue(path.isfile(self.mapped_mutation_pos))
        if generate_ref_files:
            copy(self.mapped_mutation_pos, ref_dir)
        self.assertTrue(cmp(self.mapped_mutation_pos, ref_mapped_mutation_pos))
    
    def test_pdb(self):
        # Setup temporary file (Note: generated by mu_map so shouldn't go in setup?)
        self.mapped_mutation_pos = copy(ref_mapped_mutation_pos, ymap.mapped_mutation_pos_file)
        
        c.pdb(self.pdb, self.mapped_mutation_pos, self.mapped_struct, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_struct))
        if generate_ref_files:
            copy(self.mapped_struct, ref_dir)
        self.assertTrue(cmp(self.mapped_struct, ref_mapped_struct))

    def test_mmap(self):
        c.mmap(self.mutation_prot, self.sites_id, self.mapped_sites, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_sites))
        if generate_ref_files:
            copy(self.mapped_sites, ref_dir)
        self.assertTrue(cmp(self.mapped_sites, ref_mapped_sites))
    
    def test_ptm_map(self):
        c.ptm_map(self.mutation_prot, self.ptm_id, self.mapped_ptms, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_ptms))
        if generate_ref_files:
            copy(self.mapped_ptms, ref_dir)
        self.assertTrue(cmp(self.mapped_ptms, ref_mapped_ptms))
    
    def test_dmap(self):
        c.dmap(self.mutation_prot, self.domain_id, self.mapped_domains, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_domains))
        if generate_ref_files:
            copy(self.mapped_domains, ref_dir)
        self.assertTrue(cmp(self.mapped_domains, ref_mapped_domains))
    
    def test_nucleotide_map(self):
        c.nucleotide_map(self.mutation_prot, self.nucleotide_id, self.mapped_nucleotide, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_nucleotide))
        if generate_ref_files:
            copy(self.mapped_nucleotide, ref_dir)
        self.assertTrue(cmp(self.mapped_nucleotide, ref_mapped_nucleotide))
    
    def test_summary(self):
        self.skipTest('Save time until custom test files reduce test execution time') #TODO: Remove when ready
        
        # Setup temporary file (Note: generated by mu_map so shouldn't go in setup?)
        self.mapped_mutation_pos = copy(ref_mapped_mutation_pos, ymap.mapped_mutation_pos_file)
        
        c.pdb(self.pdb, self.mapped_mutation_pos, self.mapped_struct, self.summary)
        c.mmap(self.mutation_prot, self.sites_id, self.mapped_sites, self.summary)
        c.ptm_map(self.mutation_prot, self.ptm_id, self.mapped_ptms, self.summary)
        c.dmap(self.mutation_prot, self.domain_id, self.mapped_domains, self.summary)
        c.nucleotide_map(self.mutation_prot, self.nucleotide_id, self.mapped_nucleotide, self.summary)
        
        self.assertTrue(path.isfile(self.summary))
        if generate_ref_files:
            copy(self.summary, ref_uniprot_summary)
        self.assertTrue(cmp(self.summary, ref_uniprot_summary))

        
    # functional_data()
@unittest.skip("showing class skipping")
class FunctionalMappingTest(unittest.TestCase):
    '''
    Test mapping of mutations to PTMFunc and PTMcode annotations
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.mutation_prot = copy(ref_mutation_prot, ymap.mutation_prot_file) # Input file
        self.yeastID = copy(ref_yeastID, ymap.yeastID_file)        
        self.interface_acet = copy(ref_interface_acet, ymap.interface_acet_file)
        self.interface_phos = copy(ref_interface_phos, ymap.interface_phos_file)
        self.interface_ubiq = copy(ref_interface_ubiq, ymap.interface_ubiq_file)
        self.interact_acet = copy(ref_interact_acet, ymap.interact_acet_file)
        self.interact_phos = copy(ref_interact_phos, ymap.interact_phos_file)
        self.interact_ubiq = copy(ref_interact_ubiq, ymap.interact_ubiq_file)
        self.hotspot = copy(ref_hotspot, ymap.regulatory_hotspots_file)        
        self.within_prot = copy(ref_within_prot, ymap.within_prot_file)
        self.between_prot = copy(ref_between_prot, ymap.between_prot_file)

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
        ymap.interface(self.yeastID, self.interface_acet, self.mutation_prot, self.mapped_interface_acet, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_interface_acet))
        if generate_ref_files:
            copy(self.mapped_interface_acet, ref_mapped_interface_acet)
        self.assertTrue(cmp(self.mapped_interface_acet, ref_mapped_interface_acet))
        
        
        ymap.interface(self.yeastID, self.interface_phos, self.mutation_prot, self.mapped_interface_phos, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_interface_phos))
        if generate_ref_files:
            copy(self.mapped_interface_phos, ref_mapped_interface_phos)
        self.assertTrue(cmp(self.mapped_interface_acet, ref_mapped_interface_phos))
        
        
        ymap.interface(self.yeastID, self.interface_ubiq, self.mutation_prot, self.mapped_interface_ubiq, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_interface_ubiq))
        if generate_ref_files:
            copy(self.mapped_interface_ubiq, ref_mapped_interface_ubiq)
        self.assertTrue(cmp(self.mapped_interface_acet, ref_mapped_interface_ubiq))
    
    # pi()
    def test_ppi(self):
        # For each SC_<site>_interactions_sc.txt
        ymap.ppi(self.yeastID, self.interact_acet, self.mutation_prot, self.mapped_interact_acet, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_interact_acet))
        if generate_ref_files:
            copy(self.mapped_interact_acet, ref_mapped_interact_acet)
        self.assertTrue(cmp(self.mapped_interact_acet, ref_mapped_interact_acet))
        
        ymap.ppi(self.yeastID, self.interact_phos, self.mutation_prot, self.mapped_interact_phos, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_interact_phos))
        if generate_ref_files:
            copy(self.mapped_interact_phos, ref_mapped_interact_phos)
        self.assertTrue(cmp(self.mapped_interact_acet, ref_mapped_interact_phos))
        
        ymap.ppi(self.yeastID, self.interact_ubiq, self.mutation_prot, self.mapped_interact_ubiq, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_interact_ubiq))
        if generate_ref_files:
            copy(self.mapped_interact_ubiq, ref_mapped_interact_ubiq)
        self.assertTrue(cmp(self.mapped_interact_acet, ref_mapped_interact_ubiq))
    
    # withP()    
    def test_withinPro(self):
        ymap.withinPro(self.yeastID, self.within_prot, self.mutation_prot, self.mapped_within_prot, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_within_prot))
        if generate_ref_files:
            copy(self.mapped_within_prot, ref_dir)
        self.assertTrue(cmp(self.mapped_within_prot, ref_mapped_within_prot))
    
    # betweenP()
    def test_betweenPro(self):
        ymap.betweenPro(self.yeastID, self.between_prot, self.mutation_prot, self.mapped_between_prot, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_between_prot))
        if generate_ref_files:
            copy(self.mapped_between_prot, ref_dir)
        self.assertTrue(cmp(self.mapped_between_prot, ref_mapped_between_prot))
    
    # hotS()
    def test_hotspot(self):
        ymap.hotspot(self.yeastID, self.hotspot, self.mutation_prot, self.mapped_hotspot, self.summary)
        
        self.assertTrue(path.isfile(self.mapped_hotspot))
        if generate_ref_files:
            copy(self.mapped_hotspot, ref_dir)
        self.assertTrue(cmp(self.mapped_hotspot, ref_mapped_hotspot))
    
    def test_summary(self):
        self.skipTest('Save time until custom test files reduce test execution time') #TODO: Remove when ready
        
        ymap.interface(self.yeastID, self.interface_acet, self.mutation_prot, self.mapped_interface_acet, self.summary)
        ymap.interface(self.yeastID, self.interface_phos, self.mutation_prot, self.mapped_interface_phos, self.summary)
        ymap.interface(self.yeastID, self.interface_ubiq, self.mutation_prot, self.mapped_interface_ubiq, self.summary)
        ymap.ppi(self.yeastID, self.interact_acet, self.mutation_prot, self.mapped_interact_acet, self.summary)
        ymap.ppi(self.yeastID, self.interact_phos, self.mutation_prot, self.mapped_interact_phos, self.summary)
        ymap.ppi(self.yeastID, self.interact_ubiq, self.mutation_prot, self.mapped_interact_ubiq, self.summary)
        ymap.withinPro(self.yeastID, self.within_prot, self.mutation_prot, self.mapped_within_prot, self.summary)
        ymap.betweenPro(self.yeastID, self.between_prot, self.mutation_prot, self.mapped_between_prot, self.summary)
        ymap.hotspot(self.yeastID, self.hotspot, self.mutation_prot, self.mapped_hotspot, self.summary)
        
        self.assertTrue(path.isfile(self.summary))
        if generate_ref_files:
            copy(self.summary, ref_func_summary)
        self.assertTrue(cmp(self.summary, ref_func_summary))



    # ygenes
@unittest.skip("showing class skipping")
class GeneToProteinTest(unittest.TestCase):
    '''
    Test conversion of gene-level mutation file to protein-level mutation file
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.mutation_gene = copy(ref_mutation_gene, ymap.mutation_gene_file) # Input file
        self.gff = copy(ref_gff, ymap.gff_file)
        self.d_id_map = copy(ref_d_id_map, ymap.d_id_map_file)

        # Setup temporary file paths
        self.mutation_prot = path.join(self.test_dir, ymap.mutation_prot_file)         
                    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
        
    def test_mutation_file(self):
        ymap.mutation_file(self.mutation_gene, self.gff, self.d_id_map, self.mutation_prot)
        
        self.assertTrue(path.isfile(self.mutation_prot))
        if generate_ref_files:
            copy(self.mutation_prot, ref_mutation_prot_converted)
        self.assertTrue(cmp(self.mutation_prot, ref_mutation_prot_converted))
    
    def test_revcomp(self):
        pass
    
    def test_translate_dna(self):
        pass


    
    # GO enrichment, Biogrid network and final report
class FinalAnalysisTest(unittest.TestCase):
    '''
    Test production of files containing p-values for GO enrichment and files containing Biogrid IDs
    Test conversion of summary file into final report file
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = mkdtemp()
        chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.uniprot_biogrid = copy(ref_uniprot_biogrid, ymap.uniprot_biogrid_file)
        #TODO: Don't think it's necessary to test each file separately. Only test the final report file?
#         self.mapped_interface_acet = copy(ref_mapped_interface_acet, ymap.mapped_interface_acet_file)
#         self.mapped_interface_phos = copy(ref_mapped_interface_phos, ymap.mapped_interface_phos_file)
#         self.mapped_interface_ubiq = copy(ref_mapped_interface_ubiq, ymap.mapped_interface_ubiq_file)
#         self.mapped_interact_acet = copy(ref_mapped_interact_acet, ymap.mapped_interact_acet_file)
#         self.mapped_interact_phos = copy(ref_mapped_interact_phos, ymap.mapped_interact_phos_file)
#         self.mapped_interact_ubiq = copy(ref_mapped_interact_ubiq, ymap.mapped_interact_ubiq_file)
#         self.mapped_hotspot = copy(ref_mapped_hotspot, ymap.mapped_hotspot_file)
#         self.mapped_within_prot = copy(ref_mapped_within_prot, ymap.mapped_within_prot_file)
#         self.mapped_between_prot = copy(ref_mapped_between_prot, ymap.mapped_between_prot_file)
#         self.mapped_sites = copy(ref_mapped_sites, ymap.mapped_sites_file)
#         self.mapped_ptms = copy(ref_mapped_ptms, ymap.mapped_ptms_file)
#         self.mapped_domains = copy(ref_mapped_domains, ymap.mapped_domains_file)
#         self.mapped_nucleotide = copy(ref_mapped_nucleotide, ymap.mapped_nucleotide_file)
#         self.mapped_struct = copy(ref_mapped_struct, ymap.mapped_struct_file)

        # Setup temporary file paths
        self.final_report = path.join(self.test_dir, ymap.final_report_file)
        self.p_value = path.join(self.test_dir, ymap.p_value_file)
        self.biog = path.join(self.test_dir, ymap.biog_file)
                    
    def tearDown(self):
        chdir(ref_dir) # navigate out of test directory before deleting
        rmtree(self.test_dir)
    
    def test_sum_file_map(self):
        self.summary = copy(ref_summary, ymap.summary_file)
        
        ymap.sum_file_map(self.summary, self.final_report)
        
        self.assertTrue(path.isfile(self.final_report))
        if generate_ref_files:
            copy(self.final_report, ref_final_report)
        self.assertTrue(cmp(self.final_report, ref_final_report))
    
    def test_enrich(self):
        self.final_report = copy(ref_final_report, ymap.final_report_file)
        
        c.enrich(self.final_report)
        
        self.assertTrue(path.isfile(self.p_value))
        if generate_ref_files:
            copy(self.p_value, ref_p_value)
        self.assertTrue(cmp(self.p_value, ref_p_value))
    
    def test_preWeb(self):
        self.final_report = copy(ref_final_report, ymap.final_report_file)
        
        c.preWeb(self.uniprot_biogrid, self.biog)
        
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
