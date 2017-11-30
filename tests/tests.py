from __future__ import print_function, absolute_import, division, unicode_literals
import unittest
from tempfile import mkdtemp
from shutil import copy, rmtree
from os import chdir, mkdir, makedirs, path
from time import clock
from filecmp import cmp
import ymap.ymap as ymap 

# Setup filenames
# Input files
mutation_prot = 'mutation.txt'
mutation_gene = 'mutated_proteins.txt'

# Intermediate processing files
uniprot = 'uniprot_mod_raw.txt' # downloaded
bact = 'bact.txt'
ptms = 'PTMs.txt'
domains = 'domains.txt'
nucleotide = 'nucleotide.txt'
pdb = 'pdb.txt'

gff = 'gff.txt' # downloaded
frmt = 'frmt.txt'

yeastID = 'yeastID.txt' # downloaded
d_id_map = 'd_id_map.txt' # Links Uniprot ID to Gene ID, ORF start, ORF stop, strand orientation 
sites_id = 'sites_id.txt'
ptm_id = 'ptm_id_file.txt'
domain_id = 'id_domain.txt'
nucleotide_id = 'id_nucleotide.txt' # Nucleotide binding sites

# Downloaded and copied PTMfunc and PTMcode data files
interface_acet = '3DID_aceksites_interfaceRes_sc.txt'
interface_phos = '3DID_phosphosites_interfaceRes_sc.txt'
interface_ubiq = '3DID_ubisites_interfaceRessc_sc.txt'
regulatory_hotspots = 'schotspot_updated.txt'
interact_acet = 'SC_acet_interactions.txt'
interact_phos = 'SC_psites_interactions_sc.txt'
interact_ubiq = 'SC_ubi_interactions_sc.txt'
within_prot = 'sc_within_proteins.txt'
between_prot = 'sc_btw_proteins.txt'

# Output files
summary = 'summary.txt'

mapped_sites = 'ab_mutation_file.txt'
mapped_ptms = 'mapped_ptms.txt'
mapped_domains = 'domains_mapped.txt'
mapped_nucleotide = 'nucleotide_map.txt'
mapped_mutation_pos = 'mutation_id.txt' 
mapped_struct = 'stru_mutation.txt' # Structural regions of protein e.g. beta sheet, alpha helix, turn 

mapped_interface_acet = 'interface_mutation.txt'
mapped_interface_phos = 'interface_mutation.txt'
mapped_interface_ubiq = 'interface_mutation.txt'
mapped_interact_acet = 'ppi_mutation.txt'
mapped_interact_phos = 'ppi_mutation.txt'
mapped_interact_ubiq = 'ppi_mutation.txt'
mapped_hotspot = 'hotspot.txt'
mapped_within_prot = 'within_protein.txt'
mapped_between_prot = 'ptm_between_proteins.txt'



# Reference directory paths
ref_dir = os.path.abspath('test_files')

ref_mutation_prot = os.path.join(ref_dir, mutation_prot)
ref_mutation_prot_converted = os.path.join(ref_dir, 'converted_' + mutation_prot)
ref_mutation_gene = os.path.join(ref_dir, mutation_gene)


ref_uniprot_file = os.path.join(ref_dir, uniprot)
ref_bact = os.path.join(ref_dir, bact)
ref_ptms = os.path.join(ref_dir, ptms)
ref_domains = os.path.join(ref_dir, domains)
ref_nucleotide = os.path.join(ref_dir, nucleotide)
ref_pdb = os.path.join(ref_dir, pdb)

ref_gff = os.path.join(ref_dir, gff)
ref_frmt = os.path.join(ref_dir, frmt)

ref_yeastID = os.path.join(ref_dir, yeastID)
ref_d_id_map = os.path.join(ref_dir, d_id_map)
ref_sites_id = os.path.join(ref_dir, sites_id)
ref_ptm_id = os.path.join(ref_dir, ptm_id)
ref_domain_id = os.path.join(ref_dir, domain_id)
ref_nucleotide_id = os.path.join(ref_dir, nucleotide_id)


ref_summary = os.path.join(ref_dir, summary)

ref_mapped_sites = os.path.join(ref_dir, mapped_sites)
ref_mapped_ptms = os.path.join(ref_dir, mapped_ptms)
ref_mapped_domains = os.path.join(ref_dir, mapped_domains)
ref_mapped_nucleotide = os.path.join(ref_dir, mapped_nucleotide)
ref_mapped_mutation_pos = os.path.join(ref_dir, mapped_mutation_pos)
ref_mapped_struct = os.path.join(ref_dir, mapped_struct)

ref_interface_acet = os.path.join(ref_dir, interface_acet)
ref_interface_phos = os.path.join(ref_dir, interface_phos)
ref_interface_ubiq = os.path.join(ref_dir, interface_ubiq)
ref_hotspot = os.path.join(ref_dir, regulatory_hotspots)
ref_interact_acet = os.path.join(ref_dir, interact_acet)
ref_interact_phos = os.path.join(ref_dir, interact_phos)
ref_interact_ubiq = os.path.join(ref_dir, interact_ubiq)
ref_within_prot = os.path.join(ref_dir, within_prot)
ref_between_prot = os.path.join(ref_dir, between_prot)

ref_mapped_interface_acet = os.path.join(ref_dir, 'acet_' + mapped_interface_acet)
ref_mapped_interface_phos = os.path.join(ref_dir, 'phos_' + mapped_interface_phos)
ref_mapped_interface_ubiq = os.path.join(ref_dir, 'ubiq_' + mapped_interface_ubiq)
ref_mapped_interact_acet = os.path.join(ref_dir, 'acet_' + mapped_interact_acet)
ref_mapped_interact_phos = os.path.join(ref_dir, 'phos_' + mapped_interact_phos)
ref_mapped_interact_ubiq = os.path.join(ref_dir, 'ubiq_' + mapped_interact_ubiq)
ref_mapped_hotspot = os.path.join(ref_dir, mapped_hotspot)
ref_mapped_within_prot = os.path.join(ref_dir, mapped_within_prot)
ref_mapped_between_prot = os.path.join(ref_dir, mapped_between_prot)

# Create object to allow function calls
#TODO: Remove the class from ymap.py - because it is pointless
c = YGtPM()


class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start
   
@unittest.skip("showing class skipping")
class DataDownloadTest(unittest.TestCase):
    '''
    Test ydata downloads - these tests are low priority
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        
    def tearDown(self):
        os.chdir(ref_dir) # navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)
        
    def test_gff(self):
        pass
    
    def test_iD(self):
        # This test may fail because yeastID.txt is downloaded and data may have changed from that in test_files folder.
#         c.iD()
#         self.assertTrue(os.path.isfile('yeastID.txt'))
#         self.assertTrue(filecmp.cmp('yeastID.txt', os.path.join(ref_dir, 'yeastID.txt')))
        pass
    
    def test_pTMdata(self):
        pass
    
    def test_bioGrid(self):
        pass
    
    def test_resc(self):
        # To test copying text and zip files - UNNECESSARY?
        pass
    
# ydata processing
class UniprotRawProcessingTest(unittest.TestCase):        
    '''
    Test uniprot_mod_raw.txt processing
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.uniprot_file = shutil.copy(ref_uniprot_file, self.test_dir)
        
        # Setup temporary file paths
        self.bact = os.path.join(self.test_dir, bact)
        self.ptms = os.path.join(self.test_dir, ptms)
        self.domains = os.path.join(self.test_dir, domains)
        self.nucleotide = os.path.join(self.test_dir, nucleotide)
        self.pdb = os.path.join(self.test_dir, pdb)
        
    def tearDown(self):
        os.chdir(ref_dir) # navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)
    
    def test_ab(self):
        c.ab(self.uniprot_file)
        
        self.assertTrue(os.path.isfile(self.bact))
        shutil.copy(self.bact, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.bact, ref_bact))
        
    def test_clean(self):
        c.clean(self.uniprot_file)
        
        self.assertTrue(os.path.isfile(self.ptms))
        shutil.copy(self.ptms, ref_dir) #TODO: Remove this (after file obtained)
        self.assertTrue(filecmp.cmp(self.ptms, ref_ptms))
    
    def test_dclean(self):
        c.dclean(self.uniprot_file)
        
        self.assertTrue(os.path.isfile(self.domains))
        shutil.copy(self.domains, ref_dir) #TODO: Remove this (after file obtained)
        self.assertTrue(filecmp.cmp(self.domains, ref_domains))
    
    def test_nucleotide(self):
        c.nucleotide(self.uniprot_file)
        
        self.assertTrue(os.path.isfile(self.nucleotide))
        shutil.copy(self.nucleotide, ref_dir) #TODO: Remove this (after file obtained)
        self.assertTrue(filecmp.cmp(self.nucleotide, ref_nucleotide))
    
    def test_pdb_c(self):
        c.pdb_c(self.uniprot_file)
        
        self.assertTrue(os.path.isfile(self.pdb))
        shutil.copy(self.pdb, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.pdb, ref_pdb))

    # gff.txt processing
class GffProcessingTest(unittest.TestCase):        
    '''
    Test gff.txt processing
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.gff = shutil.copy(ref_gff, self.test_dir)
        
        # Setup temporary file paths
        self.frmt = os.path.join(self.test_dir, frmt)
    
    def tearDown(self):
        os.chdir(ref_dir) # navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)
        
    def test_frmt(self):
        c.frmt(self.gff)
        
        self.assertTrue(os.path.isfile(self.frmt))
        shutil.copy(self.frmt, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.frmt, ref_frmt))
    
    # yeastID.txt processing
class YeastIDProcessingTest(unittest.TestCase):        
    '''
    Test yeastID.txt processing
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.yeastID = shutil.copy(ref_yeastID, yeastID)
        self.frmt = shutil.copy(ref_frmt, self.test_dir)
        self.bact = shutil.copy(ref_bact, bact)
        self.ptms = shutil.copy(ref_ptms, ptms)
        self.domains = shutil.copy(ref_domains, domains)
        self.nucleotide = shutil.copy(ref_nucleotide, nucleotide)
        
        # Setup temporary file paths                         
        self.d_id_map = os.path.join(self.test_dir, d_id_map)
        self.sites_id = os.path.join(self.test_dir, sites_id)
        self.ptm_id = os.path.join(self.test_dir, ptm_id)
        self.domain_id = os.path.join(self.test_dir, domain_id)
        self.nucleotide_id = os.path.join(self.test_dir, nucleotide_id)                        
    
    def tearDown(self):
        os.chdir(ref_dir) # navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)
    
    def test_id_map(self):
        c.id_map(self.yeastID, self.frmt)
        
        self.assertTrue(os.path.isfile(self.d_id_map))
        shutil.copy(self.d_id_map, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.d_id_map, ref_d_id_map))
    
    def test_id(self):
        c.id(self.bact, self.yeastID) #TODO: Again, poor consistency in the parameter order and the yeastID parameter is unused in ymap.py
        
        self.assertTrue(os.path.isfile(self.sites_id))
        shutil.copy(self.sites_id, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.sites_id, ref_sites_id))
    
    def test_pmap(self):
        c.pmap(self.yeastID, self.ptms)
        
        self.assertTrue(os.path.isfile(self.ptm_id))
        shutil.copy(self.ptm_id, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.ptm_id, ref_ptm_id))
    
    def test_d_map(self):
        c.d_map(self.yeastID, self.domains)
        
        self.assertTrue(os.path.isfile(self.domain_id))
        shutil.copy(self.domain_id, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.domain_id, ref_domain_id))
    
    def test_n_map(self):
        c.n_map(self.yeastID, self.nucleotide)
        
        self.assertTrue(os.path.isfile(self.nucleotide_id))
        shutil.copy(self.nucleotide_id, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.nucleotide_id, ref_nucleotide_id))
    
    
    
    # yproteins
    # uniprot_data()
@unittest.skip("showing class skipping")
class UniprotMappingTest(unittest.TestCase):
    '''
    Test mapping of mutations to UniProt annotations
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.mutation_prot = shutil.copy(ref_mutation_prot, mutation_prot) # Input file
        self.yeastID = shutil.copy(ref_yeastID, yeastID)
        self.d_id_map = shutil.copy(ref_d_id_map, d_id_map)
        self.sites_id = shutil.copy(ref_sites_id, sites_id)
        self.ptm_id = shutil.copy(ref_ptm_id, ptm_id)
        self.domain_id = shutil.copy(ref_domain_id, domain_id)
        self.nucleotide_id = shutil.copy(ref_nucleotide_id, nucleotide_id)
        self.pdb = shutil.copy(ref_pdb, pdb)
        
        # Setup temporary file paths
        self.summary = os.path.join(self.test_dir, summary)         
        self.mapped_sites = os.path.join(self.test_dir, mapped_sites)
        self.mapped_ptms = os.path.join(self.test_dir, mapped_ptms)
        self.mapped_domains = os.path.join(self.test_dir, mapped_domains)
        self.mapped_nucleotide = os.path.join(self.test_dir, mapped_nucleotide)
        self.mapped_mutation_pos = os.path.join(self.test_dir, mapped_mutation_pos)
        self.mapped_struct = os.path.join(self.test_dir, mapped_struct)
                    
    def tearDown(self):
        os.chdir(ref_dir) # navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)    
        
    def test_mu_map(self):
        c.mu_map(self.yeastID, self.mutation_prot)
        
        self.assertTrue(os.path.isfile(self.mapped_mutation_pos))
        shutil.copy(self.mapped_mutation_pos, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_mutation_pos, ref_mapped_mutation_pos))
    
    def test_pdb(self):
        # Setup temporary file (Note: generated by mu_map so shouldn't go in setup?)
        self.mapped_mutation_pos = shutil.copy(ref_mapped_mutation_pos, mapped_mutation_pos)
        
        c.pdb(self.pdb, self.mapped_mutation_pos)
        
        self.assertTrue(os.path.isfile(self.mapped_struct))
        shutil.copy(self.mapped_struct, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_struct, ref_mapped_struct))
        #TODO: Not sure how to implement check of summary...
#         self.assertTrue(os.path.isfile(self.summary)) 
#         shutil.copy(self.summary, ref_dir) #TODO: Remove this (after file initially obtained)
#         self.assertTrue(filecmp.cmp(self.summary, ref_summary))

    def test_mmap(self):
        c.mmap(self.mutation_prot, self.sites_id)
        
        self.assertTrue(os.path.isfile(self.mapped_sites))
        shutil.copy(self.mapped_sites, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_sites, ref_mapped_sites))
    
    def test_ptm_map(self):
        c.ptm_map(self.mutation_prot, self.ptm_id)
        
        self.assertTrue(os.path.isfile(self.mapped_ptms))
        shutil.copy(self.mapped_ptms, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_ptms, ref_mapped_ptms))
    
    def test_dmap(self):
        c.dmap(self.mutation_prot, self.domain_id)
        
        self.assertTrue(os.path.isfile(self.mapped_domains))
        shutil.copy(self.mapped_domains, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_domains, ref_mapped_domains))
    
    def test_nucleotide_map(self):
        c.nucleotide_map(self.mutation_prot, self.nucleotide_id)
        
        self.assertTrue(os.path.isfile(self.mapped_nucleotide))
        shutil.copy(self.mapped_nucleotide, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_nucleotide, ref_mapped_nucleotide))
    

        
    # functional_data()
@unittest.skip("showing class skipping")
class FunctionalMappingTest(unittest.TestCase):
    '''
    Test mapping of mutations to PTMFunc and PTMcode annotations
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.mutation_prot = shutil.copy(ref_mutation_prot, mutation_prot) # Input file
        self.yeastID = shutil.copy(ref_yeastID, yeastID)        
        self.interface_acet = shutil.copy(ref_interface_acet, interface_acet)
        self.interface_phos = shutil.copy(ref_interface_phos, interface_phos)
        self.interface_ubiq = shutil.copy(ref_interface_ubiq, interface_ubiq)
        self.interact_acet = shutil.copy(ref_interact_acet, interact_acet)
        self.interact_phos = shutil.copy(ref_interact_phos, interact_phos)
        self.interact_ubiq = shutil.copy(ref_interact_ubiq, interact_ubiq)
        self.hotspot = shutil.copy(ref_hotspot, regulatory_hotspots)        
        self.within_prot = shutil.copy(ref_within_prot, within_prot)
        self.between_prot = shutil.copy(ref_between_prot, between_prot)

        # Setup temporary file paths
        self.summary = os.path.join(self.test_dir, summary)         
        self.mapped_interface_acet = os.path.join(self.test_dir, mapped_interface_acet)
        self.mapped_interface_phos = os.path.join(self.test_dir, mapped_interface_phos)
        self.mapped_interface_ubiq = os.path.join(self.test_dir, mapped_interface_ubiq)
        self.mapped_interact_acet = os.path.join(self.test_dir, mapped_interact_acet)
        self.mapped_interact_phos = os.path.join(self.test_dir, mapped_interact_phos)
        self.mapped_interact_ubiq = os.path.join(self.test_dir, mapped_interact_ubiq)
        self.mapped_hotspot = os.path.join(self.test_dir, mapped_hotspot)
        self.mapped_within_prot = os.path.join(self.test_dir, mapped_within_prot)
        self.mapped_between_prot = os.path.join(self.test_dir, mapped_between_prot)
                    
    def tearDown(self):
        os.chdir(ref_dir) # navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)
        
    # intf()    
    def test_interface(self):
        # For each 3DID_<site>_interfaceRes_sc.txt file
        interface(self.interface_acet, self.mutation_prot)
        
        interface(self.interface_phos, self.mutation_prot)
        
        interface(self.interface_ubiq, self.mutation_prot)
        
        self.assertTrue(os.path.isfile(self.mapped_interface_acet))
        self.assertTrue(os.path.isfile(self.mapped_interface_phos))
        self.assertTrue(os.path.isfile(self.mapped_interface_ubiq))
        
        shutil.copy(self.mapped_interface_acet, ref_mapped_interface_acet) #TODO: Remove this (after file initially obtained)
        shutil.copy(self.mapped_interface_phos, ref_mapped_interface_phos) #TODO: Remove this (after file initially obtained)
        shutil.copy(self.mapped_interface_ubiq, ref_mapped_interface_ubiq) #TODO: Remove this (after file initially obtained)
        
        self.assertTrue(filecmp.cmp(self.mapped_interface_acet, ref_mapped_interface_acet))
        self.assertTrue(filecmp.cmp(self.mapped_interface_acet, ref_mapped_interface_phos))
        self.assertTrue(filecmp.cmp(self.mapped_interface_acet, ref_mapped_interface_ubiq))
    
    # pi()
    def test_ppi(self):
        # For each SC_<site>_interactions_sc.txt
        ppi(self.interact_acet, self.mutation_prot)
        
        ppi(self.interact_phos, self.mutation_prot)
        
        ppi(self.interact_ubiq, self.mutation_prot)
        
        self.assertTrue(os.path.isfile(self.mapped_interact_acet))
        self.assertTrue(os.path.isfile(self.mapped_interact_phos))
        self.assertTrue(os.path.isfile(self.mapped_interact_ubiq))
        
        shutil.copy(self.mapped_interact_acet, ref_mapped_interact_acet) #TODO: Remove this (after file initially obtained)
        shutil.copy(self.mapped_interact_phos, ref_mapped_interact_phos) #TODO: Remove this (after file initially obtained)
        shutil.copy(self.mapped_interact_ubiq, ref_mapped_interact_ubiq) #TODO: Remove this (after file initially obtained)
        
        self.assertTrue(filecmp.cmp(self.mapped_interact_acet, ref_mapped_interact_acet))
        self.assertTrue(filecmp.cmp(self.mapped_interact_acet, ref_mapped_interact_phos))
        self.assertTrue(filecmp.cmp(self.mapped_interact_acet, ref_mapped_interact_ubiq))
    
    # withP()    
    def test_withinPro(self):
        withinPro(self.within_prot, self.mutation_prot)
        
        self.assertTrue(os.path.isfile(self.mapped_within_prot))
        shutil.copy(self.mapped_within_prot, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_within_prot, ref_mapped_within_prot))
    
    # betweenP()
    def test_betweenPro(self):
        betweenPro(self.between_prot, self.mutation_prot)
        
        self.assertTrue(os.path.isfile(self.mapped_between_prot))
        shutil.copy(self.mapped_between_prot, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_between_prot, ref_mapped_between_prot))
    
    # hotS()
    def test_hotspot(self):
        hotspot(self.hotspot, self.mutation_prot)
        
        self.assertTrue(os.path.isfile(self.mapped_hotspot))
        shutil.copy(self.mapped_hotspot, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mapped_hotspot, ref_mapped_hotspot))
    
            
    # ygenes
@unittest.skip("showing class skipping")
class GeneToProteinTest(unittest.TestCase):
    '''
    Test conversion of gene-level mutation file to protein-level mutation file
    '''
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        
        # Copy reference files to temporary directory
        self.mutation_gene = shutil.copy(ref_mutation_gene, mutation_gene) # Input file
        self.gff = shutil.copy(ref_gff, gff)
        self.d_id_map = shutil.copy(ref_d_id_map, d_id_map)

        # Setup temporary file paths
        self.mutation_prot = os.path.join(self.test_dir, mutation_prot)         
                    
    def tearDown(self):
        os.chdir(ref_dir) # navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)
        
    def test_mutation_file(self):
        mutation_file(self.mutation_gene, self.d_id_map)
        
        self.assertTrue(os.path.isfile(self.mutation_prot))
        shutil.copy(self.mutation_prot, ref_mutation_prot_converted) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.mutation_prot, ref_mutation_prot_converted))
    
    def test_revcomp(self):
        pass
    
    def test_translate_dna(self):
        pass

    
    
    # Output functions
    def test_enrich(self):
        pass
    
    def test_preWeb(self):
        pass
    
    def test_sum_file_map(self):
        pass
    
    
    # yweb
    def test_bweb(self):
        pass


# NOTE: Run from command line. cd to root ymap directory and run 'python -m unittest tests.tests' 
if __name__ == '__main__':
    unittest.main()
