from __future__ import print_function, absolute_import, division, unicode_literals
import unittest
import tempfile
import shutil
import os
import time
import filecmp
from ymap.ymap import data, ymap_proteins, ymap_genes, web, YGtPM

ref_dir = os.path.abspath('test_files')


class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start


class YGtPMTest(unittest.TestCase):
    """
    YMap tests
    """

    def setUp(self):
        # Create a temporary directory
        # self.ref_dir = os.path.abspath('test_files')
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        self.c = YGtPM()
        # self.store_ref = False
        # at various places output is compared to stored reference files. put this to true
        # to regenerate all reference files

    def tearDown(self):
        # Remove the directory after the test
        os.chdir(ref_dir) # need to navigate out of test directory before deleting
        shutil.rmtree(self.test_dir)

    def test_class_creation(self):
        """
        testing the creation of YGtPM
        """
        self.assertIsInstance(self.c, YGtPM)
        self.c.pTMdata()
        self.assertTrue(os.path.isfile('uniprot_mod_raw.txt'))

    def test_data(self):
        """
        testing the steps in data()
        """

        self.assertIsInstance(self.c, YGtPM)
        try:  # the actual download can only be performed a limit number of times per hour
            raise RuntimeError  # this is to prevent connecting the database during development comment for real tesing
            self.c.pTMdata()
            if self.store_ref:
                shutil.copy('uniprot_mod_raw.txt', ref_dir)  # saving the reference data
        except:
            shutil.copy(os.path.join(ref_dir, 'uniprot_mod_raw.txt'), '.')
        self.assertTrue(os.path.isfile('uniprot_mod_raw.txt'))

        self.c.clean('uniprot_mod_raw.txt')  # produces PTMs.txt
        self.assertTrue(os.path.isfile('PTMs.txt'))
        self.assertTrue(filecmp.cmp('PTMs.txt', os.path.join(ref_dir, 'PTMs.txt')))

        # This test may fail because yeastID.txt is downloaded and data may have changed from that in test_files folder.
        self.c.iD()
        self.assertTrue(os.path.isfile('yeastID.txt'))
        self.assertTrue(filecmp.cmp('yeastID.txt', os.path.join(ref_dir, 'yeastID.txt')))

        self.c.pmap('yeastID.txt', 'PTMs.txt')
        self.assertTrue(os.path.isfile('PTM_id_file.txt'))
        self.assertTrue(filecmp.cmp('PTM_id_file.txt', os.path.join(ref_dir, 'PTM_id_file.txt')))
        
        
    # 1 TEST PER FUNCTION
    # ydata downloads - these tests are low priority
    def test_gff(self):
        pass
    
    def test_iD(self):
        # See Michiel's implementation
        pass
    
    def test_pTMdata(self):
        pass
    
    def test_bioGrid(self):
        pass
    
    def test_resc(self):
        # To test copying text and zip files - UNNECESSARY?
        pass
    
    
    # ydata processing
    # uniprot_mod_raw.txt processing
    def test_ab(self):
        pass
    
    def test_clean(self):
        pass
    
    def test_dclean(self):
        pass
    
    def test_nucleotide(self):
        pass 
    
    # ydata processing
    # gff.txt processing
    def test_frmt(self):
        pass
    
    # yeastID.txt processing
    def test_id_map(self):
        pass
    
    def test_id(self):
        pass
    
    def test_pmap(self):
        # See Michiel's implementation
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
