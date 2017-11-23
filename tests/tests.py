from __future__ import print_function, absolute_import, division, unicode_literals
import unittest
import tempfile
import shutil
import os
import time
import filecmp
from ymap.ymap import data, ymap_proteins, ymap_genes, web, YGtPM

ref_dir = os.path.abspath('test_files')

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
        # Setup filenames and paths
        uniprot = 'uniprot_mod_raw.txt'
        bact = 'bact.txt'
        ptms = 'PTMs.txt'
        domains = 'domains.txt'
        nucleotide = 'nucleotide.txt'
        
        self.ref_bact = os.path.join(ref_dir, bact)
        self.bact = os.path.join(self.test_dir, bact)
        self.ref_ptms = os.path.join(ref_dir, ptms)
        self.ptms = os.path.join(self.test_dir, ptms)
        self.ref_domains = os.path.join(ref_dir, domains)
        self.domains = os.path.join(self.test_dir, domains)
        self.ref_nucleotide = os.path.join(ref_dir, nucleotide)
        self.nucleotide = os.path.join(self.test_dir, nucleotide)
        
        # Copy uniprot_mod_raw.txt to temporary directory
        self.ref_uniprot_file = os.path.join(ref_dir, uniprot)
        self.uniprot_file = shutil.copy(self.ref_uniprot_file, self.test_dir)
        
    def tearDown(self):
        os.chdir(ref_dir) # navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)
    
    def test_ab(self):
        c.ab(self.uniprot_file)
        
        self.assertTrue(os.path.isfile(self.bact))
#         shutil.copy(self.bact, ref_dir) #TODO: Remove this (after file initially obtained)
        self.assertTrue(filecmp.cmp(self.bact, self.ref_bact))
        
    def test_clean(self):
        c.clean(self.uniprot_file)
        
        self.assertTrue(os.path.isfile(self.ptms))
#         shutil.copy(self.ptms, ref_dir) #TODO: Remove this (after file obtained)
        self.assertTrue(filecmp.cmp(self.ptms, self.ref_ptms))
    
    def test_dclean(self):
        c.dclean(self.uniprot_file)
        
        self.assertTrue(os.path.isfile(self.domains))
#         shutil.copy(self.domains, ref_dir) #TODO: Remove this (after file obtained)
        self.assertTrue(filecmp.cmp(self.domains, self.ref_domains))
    
    def test_nucleotide(self):
        c.nucleotide() #TODO: Modify nucleotide() in ymap.py to take uniprot_file as argument
        
        self.assertTrue(os.path.isfile(self.nucleotide))
#         shutil.copy(self.nucleotide, ref_dir) #TODO: Remove this (after file obtained)
        self.assertTrue(filecmp.cmp(self.nucleotide, self.ref_nucleotide))
    

    # gff.txt processing
    def test_frmt(self):
        pass
    
    # yeastID.txt processing
    def test_id_map(self):
        pass
    
    def test_id(self):
        pass
    
    def test_pmap(self):
#         self.c.pmap('yeastID.txt', 'PTMs.txt')
#         self.assertTrue(os.path.isfile('PTM_id_file.txt'))
#         self.assertTrue(filecmp.cmp('PTM_id_file.txt', os.path.join(ref_dir, 'PTM_id_file.txt')))
        pass
    
    def test_d_map(self):
        pass
    
    def test_n_map(self):
        pass
    
    def test_pdb_c(self):
        pass
    
    def test_mu_map(self):
        pass
    
    # ygenes
    # Conversion of gene-level mutation file
    def test_mutation_file(self):
        pass
    
    def test_revcomp(self):
        pass
    
    def test_translate_dna(self):
        pass
    
    
    # yproteins
    # uniprot_data()
    def test_mmap(self):
        pass
    
    def test_ptm_map(self):
        pass
    
    def test_dmap(self):
        pass
    
    def test_nucleotide_map(self):
        pass
    
    def test_pdb(self):
        pass
        
    # functional_data()
    # intf()
    def test_interface(self):
        # For each 3DID_<site>_interfaceRes_sc.txt file
        pass
    
    # pi()
    def test_ppi(self):
        # For each SC_<site>_interactions_sc.txt
        pass
    
    # withP()    
    def test_withinPro(self):
        pass
    
    # betweenP()
    def test_betweenPro(self):
        pass
    
    # hotS()
    def test_hotspot(self):
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
