# yMap - Yeast Genotype to Phenotype Map

yMap is a python based application to map yeast (*Saccharomyces cerevisiae*) mutations/variants, at either protein or DNA level, to 

	- protein post-translational modifications (PTMs)
	- protein domains
	- protein-nucleotide binding domains
	- protein structural regions
	- protein active and binding sites
	- protein networks visualisation


In three user friendly steps, yMap generates a folder (with sub-directories) containing output files that can be used in further data analyses. A summary file (final_report.txt) reports all the non-synonymous mutations that overlap with or fall within proteins functional regions (described above). An enrichment file (pvalue.txt) gives p-values for significantly enriched gene ontology (GO) terms among the set of mutated proteins or genes. Another file (biog.txt) contains [BioGrid](https://thebiogrid.org) IDs for visualisation of protein networks.

Data used by yMap were downloaded from [UniProt](https://www.uniprot.org/), [*Saccharomyces* Genome Database (SGD)](https://www.yeastgenome.org/), and sources with annotated PTMs like [PTMcode 2.0](http://ptmcode.embl.de/) and [PTMfunc](http://ptmfunc.com/). For more details, see [Introduction to data used by yMap](#1-\--introduction-to-data-used-by-ymap).

See also [yMoose](https://github.com/eregenyi/yWebsite), or **y**east **M**apper **O**nline **O**pen-**S**ource **E**nhanced, a project providing a web-interface for yMap.

## Dependencies

yMap depends on:

- python 2.6.x
- python 3.x
- [Orange bioinformatics](http://pythonhosted.org/Orange-Bioinformatics/#installation)


## Installation
#### 1. Via wheel (.whl) file (RECOMMENDED METHOD)

If you are using Anaconda, open Anaconda Prompt. To avoid errors during the following installation process, it is recommended that you first create an environment for ymap (using `conda create -n <env_name>` where <env_name> is an environment name of your choosing). Activate this environment (using `activate <env_name>` in Windows, or `source activate <env_name>` in Linux/MacOS). Next, install pip in this new environment (using `conda install pip`). Finally, install Orange3 as directed [here](https://orange.biolab.si/download/) (using instructions for Anaconda, if possible).

Now we create the wheel (.whl) file. Change directory to the root directory of downloaded yMap package (directory containing setup.py file). Create the .whl file with the following command:

	python setup.py bdist_wheel --universal

This makes two folders, called 'build' and 'dist' in the current directory. The .whl file to be used for installation is found in 'dist'. Navigate to the 'dist' directory and run:

	pip install wheel_filename.whl

yMap should now be installed.


#### 2. Via pip (DO NOT USE)

Until we update the PyPi repository, pip will install a very outdated version of yMap. Hence, **DO NOT USE** the following method:

	pip install ymap

### PyPi 
https://pypi.python.org/pypi/ymap

## Video demo (OUTDATED)

[![yMap tutorial video](http://img.youtube.com/vi/pcmkuWvLRzI/0.jpg)](https://www.youtube.com/watch?v=pcmkuWvLRzI)

# Usage

### Run from installed yMap:

In the prompt/terminal, navigate to a newly created directory and run the following commands:

**Step 1** - Download prerequisite data:  
`$ ydata`  
Downloads all data needed for proper execution of ymap. WARNING: May take some time (5-20 min).  
If you have done this for a previous analysis and are happy to use old data, you can skip this step to save time, as long as you manually copy the data files from an old `ymap_data` directory (used for an old analysis) to `./ymap_data` in the current working directory (for the current analysis). If this is the first time you are using ymap, you should not skip this step.

**Step 2** - Supply an input mutations file:    
Copy and paste the protein-level mutation file (mutations.txt) OR DNA-level mutation file (mutated_proteins.txt) into `./input` (a subdirectory called input in the current working directory). See `example_mutation_files/` subdirectory of yMap package for example input file formats.

**Step 3** - Map the mutations:    
If using mutations.txt (mutations at proteins level): `$ yproteins`  
If using mutated_proteins.txt (mutations at the chromosome level with genetic coordinates): `$ ygenes`
Populates `./output` with output files (mapping mutations to protein PTMs and functional regions).

**Step 4** - Visualize protein networks on BioGrid website: `$ yweb`    
You will be prompted with `Enter a path to a biog.txt file: `. Specify the absolute path to a biog.txt file (see output of Step 3). This command will then open tabs in a web browser, each corresponding to the BioGrid page (see [example](https://thebiogrid.org/31517)) of one of the mutated proteins identified in a biog.txt output file. Allows visualization of protein networks that involve the mutated protein.


### Run from source code:

In the prompt/terminal, navigate to a newly created directory. Copy ymap.py into this directory. Then run the following (corresponding to steps described above):

**Step 1**: `$ python ymap.py -d`  

**Step 3**: `$ python ymap.py -p`  

**Step 3**: `$ python ymap.py -g`  

**Step 4**: `$ python ymap.py -w <path_to_biog.txt>`    
User must specify a single argument: an appropriately quoted path to a biog.txt file e.g. "./path/to/biog.txt" in Linux/MacOS or "C:\path\to\biog.txt" in Windows).


# Detailed Guide:
1 - [Introduction to data used by yMap](#1-\--introduction-to-data-used-by-ymap)  
2 - [Interpreting output](#2-\--interpreting-output)  
3 - [Documentation for ymap.py functions](#3-\--documentation-for-ymappy-functions)  
4 - [Troubleshooting](#4-\--troubleshooting)  


## 1 - Introduction to data used by yMap

### Input files

A tab-delimited text file called "mutated_proteins.txt" or "mutation.txt".

- **mutation.txt** should contain columns (1) protein common names and (2) mutated residue positions.  
- **mutated\_proteins.txt** should contain columns (1) chromosome number e.g. chrXI, (2) index (location) of point mutation, (3) original nucleotide e.g. C, (4) mutated nucleotide e.g. T, (5) common gene name e.g. YPK1 or INTERGENIC.

*Note: Headings in these files are not currently supported.*

### Downloaded data files (during execution of `ydata`)

#### (i) Files downloaded from [UniProt](https://www.uniprot.org/)

1. **uniprot_mod_raw.txt** - Mapping all yeast proteins to protein features, in [GFF format](https://www.ensembl.org/info/website/upload/gff.html)  
2. **uniprot_bioGrid.txt** - Mapping all yeast proteins to BioGrid IDs  
3. **yeastID.txt** - Mapping all yeast proteins to ordered locus names (OLN) and primary gene names  

#### (ii) Files downloaded from [*Saccharomyces* Genome Database (SGD)](https://www.yeastgenome.org/)

1. **gff.txt** - All yeast chromosome features, in [GFF format](https://www.ensembl.org/info/website/upload/gff.html) (to map DNA level mutations to proteins)  

#### (iii) Files from [PTMcode 2.0](http://ptmcode.embl.de/)

1. **sc_btw_proteins.zip** (containing sc_btw_proteins.txt) - Pairwise interactions between different PTMs on different yeast proteins  
2. **sc_within_proteins.txt** - Pairwise interactions between PTMs within the same yeast protein  
3. **schotspot_updated.txt** - Yeast protein hotspots (small regions in the protein sequence that are enriched in PTMs).  

#### (iv) Files from [PTMfunc](http://ptmfunc.com/)

1. **3DID_aceksites_interfaceRes_sc.txt** - Acetylated residues predicted to be at an interface, based on domain-domain interactions in [3DID database](https://3did.irbbarcelona.org/)  
2. **3DID_phosphosites_interfaceRes_sc.txt** - Phosphorylated residues predicted to be at an interface, based on domain-domain interactions in [3DID database](https://3did.irbbarcelona.org/)  
3. **3DID_ubisites_interfaceRessc_sc.txt** - Ubiquintinated residues predicted to be at an interface, based on domain-domain interactions in [3DID database](https://3did.irbbarcelona.org/)  
4. **SC_acet_interactions.txt** - Acetylated residues predicted to be at an interface, based on docking and homology models or X-ray crystallography  
5. **SC_psites_interactions_sc.txt** - Phosphorylated residues predicted to be at an interface, based on docking and homology models or X-ray crystallography  
6. **SC_ubi_interactions_sc.txt** - Ubiquintinated residues predicted to be at an interface, based on docking and homology models or X-ray crystallography  

### Processed data files (Output of `ydata`)

1.  **PTM_id_file.txt**  
2.  **PTMs.txt**  
3.  **bact.txt**  
4.  **between_prot_id.txt**  
5.  **d_id_map.txt**  
6.  **domains.txt**  
7.  **frmt.txt**  
8.  **id_domain.txt**  
9.  **id_nucleotide.txt**  
10. **id_struct.txt**  
11. **interact_acet_id.txt**  
12. **interact_phos_id.txt**  
13. **interact_ubiq_id.txt**  
14. **interface_acet_id.txt**  
15. **interface_phos_id.txt**  
16. **interface_ubiq_id.txt**  
17. **nucleotide.txt**  
18. **pdb.txt**  
19. **processed_file_names.txt**  
20. **regulatory_hotspots_id.txt**  
21. **sites_id.txt**  
22. **within_prot_id.txt**  


## 2 - Interpreting Output

Within the output directory `./output/yMap-resultsXXXXX` reside a number of files and subdirectories. Each subdirectory contains three files:

1. File mapping mutations to protein features (mutated proteins, mutation positions, mutated functional region and data source)
2. File giving GO (pathways) enrichment analysis results: pvalue.txt
3. File mapping mutated proteins to BioGrid IDs: biog.txt

Additionally, in the root folder `./output/yMap-resultsXXXXX` can be found a summary file: final_report.txt
This file simply contains all the lines from the individual mapping files (found in the subdirectories of the root folder). In the future, we hope to unify its format and make it easily machine-readable.

An error file (errors.txt) containing unrecognized gene names/failed mappings is also included in the `./output` directory.  

All files have tab-separated values and no headings.

### Full description of mapping files

**/A-B-sites/ab_mutation_file.txt**  
Maps mutated positions in proteins to active sites and binding sites.  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) Binding site/Active site, (5) Mutated position, (6) Binds to/Activity, (7) Data source

**/Domains/domains_mapped.txt**  
Maps mutated positions in proteins to domains.  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) Mutated position, (5) Domain start, (6) Domain end, (7) Domain name, (8) PROSITE reference, (9) Data source

**/Nucleotide_binding/nucleotide_map.txt**  
Maps mutated positions in proteins to nucleotide binding sites.  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) 'Nucleotide binding', (5) Mutated position, (6) Nucleotide binding site start, (7) Nucleotide binding site end, (8) Bound nucleotide, (9) Data source

**/PDB/stru_mutation.txt**  
Maps mutated positions in proteins to secondary protein structures.  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) Type of secondary structure, (5) Mutated position, (6) Secondary structure start, (7) Secondary structure end, (8) PDB ID, (9) Data source

#### Output on PTMs  

**/PTMs/mapped_ptms.txt**  
Maps mutated positions in proteins to PTMs.  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) Mutated position, (5) PTM, (9) Data source

**/PPI/Acetylation/ppi_mutation.txt**  
Maps mutated positions in proteins to acetylated residues predicted to interface with a second protein.  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) Mutated position, (5) Affected residue, (6) Uniprot ID (Protein 2), (7) OLN (Protein 2), (8) Common/Primary gene name (Protein 2), (9) Prediction method, (10) Data source

**/PPI/Phosphorylation/ppi_mutation.txt**  
As above, but for phosphorylated residues.

**/PPI/Ubiquitination/ppi_mutation.txt**  
As above, but for ubiquitinated residues.

**/Interface/Acetylation/interface_mutation.txt**  
Maps mutated positions in proteins to acetylated residues predicted to interface with another protein, based on domain-domain interactions in 3DID database.  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) Mutated position, (5) Affected residue, (6) PFAM domain, (7) Data source

**/Interface/Ubiquitination/interface_mutation.txt**  
As above, but for ubiquitinated residues.

**/Interface/Phosphorylation/interface_mutation.txt**  
As above, but for phosphorylated residues.

**/PTMs_hotSpot/hotspot.txt**  
Maps mutated positions in proteins to PTM hotspots (See [Beltrao et al. 2012](https://doi.org/10.1016/j.cell.2012.05.036)).  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) Mutated position, (5) Affected residue, (6) PFAM domain, (7) PDB ID, (8) Data source

**/PTMs_between_Proteins/ptm_between_proteins.txt**  
Maps mutated positions in proteins to PTMs involved in crosstalk between proteins. (See [Minguez et al. 2015](https://dx.doi.org/10.1093%2Fnar%2Fgku1081))  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) PTM, (5) Mutated position, (6) Affected residue, (7) Uniprot ID (Protein 2), (8) OLN (Protein 2), (9) Common/Primary gene name (Protein 2), (10) PTM (Protein 2), (11) PTM position (Protein 2), (12) Affected residue (Protein 2), (13) Data source

**/PTMs_within_Proteins/ptm_within_proteins.txt**  
Maps mutated positions in proteins to PTMs involved in crosstalk within proteins. (See [Minguez et al. 2015](https://dx.doi.org/10.1093%2Fnar%2Fgku1081))  
Columns: (1) Uniprot ID, (2) OLN, (3) Common/Primary gene name, (4) PTM, (5) Mutated position, (6) Affected residue, (7) PTM 2, (8) PTM position 2, (9) Affected residue 2, (13) Data source


## 3 - Documentation for ymap.py functions

**ab**(mutations, biogrid_IDs)  
Map mutations to active/binding sites.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**betweenP**(mutations, biogrid_IDs)  
Map mutations to PTMs between interacting proteins. See PTMcode2 for further details.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**betweenPro**(mut_prot_input, between_prot_id_input, mapped_between_prot_output, summary_output)  
Map the positions of mutations in mutated proteins to pairwise combinations of PTMs known/predicted to affect protein-protein interactions (PPIs)????  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
between_prot_id_input -- file path, mapping UniProt ID, ordered locus name, common name, PTM position, residue to a PTM on another protein  
mapped_between_prot_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name and inter-protein PTM data
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**betweenPro_map**(gene_names_by_locus, between_prot_input, between_prot_id_output)  
Map UniProt IDs to protein names and pairwise combinations of PTMs known/predicted to affect protein-protein interactions (PPIs)????  

See Minguez et al. "PTMcode v2: a resource for functional associations of post-translational modifications within and between proteins." Nucleic Acids Res. 2015 Jan 28; 43 for further details.  

*Arguments*:  
gene_names_by_locus -- dictionary, key = locus_name, value = list of gene names  
between_prot_input -- file path, mapping locus name to protein-protein interaction (PPI) data  
between_prot_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and PPI data  

___

**bioGrid**(uniprot_biogrid_output)  
Download a list of UniProt IDs and their associated BioGrid IDs.  

___

**console_web**()  
Wrapper of web() function for console_script functionality.  

___

**d_map**(yeastID_input, domains_input, domain_id_output)  
Map UniProt IDs to protein names and domains and PROSITE references  

*Arguments*:  
yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names  
domains_input -- file path, mapping UniProt ID to domains, their start and end positions in the polypeptide and PROSITE references
domain_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and domain data  

___

**data**()  
Process data files into intermediate files, including only relevant data and  

___

**dmap**(mut_prot_input, domain_id_input, mapped_domains_output, summary_output)  
Map the positions of mutations in mutated proteins to domains.  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
domain_id_input -- file path, mapping UniProt ID, ordered locus name, common name, domain position (start and end) and domain name  
mapped_domains_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, domain position and domain name  
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**domain**(mutations, biogrid_IDs)  
Map mutations to domains.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**download**()  
Download/Copy all required data into a specified directory.  

___

**enrich**(mapped_mut_input)  
Perform Gene Ontology (GO) enrichment on a set of genes (obtained from the second column of an tab-delimited input file).  

Write the results of the enrichment to a file. If a GO term has been found to be overrepresented, then its ID, name, the significance of the enrichment (for further details on how p-value is calculated see ?????????????), a list of the genes in the gene set that are annotated with the enriched GO term, and a reference count (???????) are written as a line in the file.  

___

**extractZips**(dir, extract_to_dir)  
Recursively search given dir for zip files and extract them to extract_to_dir.  

___

**frmt**(gff_input, frmt_output)  
Reformat a General Feature Format (GFF) file into a tab-delimited file (tsv) with columns gene id, start and end with strand orientation.  

___

**functional_data**(mutations, biogrid_IDs)  
Call all functions that work with PTMfunc and PTMcode data.  

___

**gff**(gff_output)  
Downloads the current General Feature Format (GFF) file for the Saccharomyces cerevisiae genome  

___

**hotS**(mutations, biogrid_IDs)  
Map mutations to PTMs in PTM hotspots. See PTMfunc for further details.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**hotspot**(mut_prot_input, hotspot_id_input, mapped_hotspot_output, summary_output)  
Map the positions of mutations in mutated proteins to PTM hotspots.  

PTM-containing motifs in close proximity are named hotspots??
See Beltrao et al. Cell 2012 for further details.  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
hotspot_id_input -- file path, mapping UniProt ID, ordered locus name, common name and hotspot data, including PTM position, residue, PFAM domain and PDB ID  
mapped_hospot_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name and hotspot data
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**hotspot_map**(gene_names_by_locus, hotspot_input, hotspot_id_output)  
Map UniProt IDs to protein names and PTM hotspots.  

PTM-containing motifs in close proximity are named hotspots??
See Beltrao et al. Cell 2012 for further details.  

*Arguments*:  
gene_names_by_locus -- dictionary, key = locus_name, value = list of gene names  
hotspot_input -- file path, mapping locus name to hotspot data, including PTM position, modified residue, PFAM domain, PDB ID  
hotspot_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and hotspot data  

___

**iD**(yeastID_output)  
Download up-to-date UniProt data file using the following query settings: Species: Saccharomyces cerevisiae, ...  

___

**id**(yeastID_input, bact_input, sites_id_output)  
Map UniProt IDs to protein names and binding sites and active sites.  

*Arguments*:  
yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names  
bact_input -- file path, mapping UniProt ID to binding sites and active site, their position in the polypeptide and their binding partner/activity  
sites_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and binding/active site data  

___

**id_map**(gene_names_by_locus, frmt_input, d_id_map_output)  
Map UniProt IDs to gene names and genomic loci (start, end positions and strand orientation).  

gene_names_by_locus -- dictionary, key = locus_name, value = list of gene names  
frmt_input -- file path, mapping gene (locus) name to start and end positions of locus (on a chromosome) and strand orientation)  
d_id_map_output -- file path, mapping UniProt IDs to gene names and genomic locus data  

___

**interface**(mut_prot_input, interface_id_input, mapped_interface_output, summary_output)  
Map the positions of mutations in mutated proteins to PTMs at an interface.  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
interface_id_input -- file path, mapping UniProt ID, ordered locus name, common name, PTM position, residue and PFAM domain name  
mapped_interface_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name and interface data  
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**interface_map**(gene_names_by_locus, interface_sites_input, interface_id_output)  
Map UniProt IDs to protein names and PTMs at an interface that are known to affect interactions.  

See Beltrao et al. Cell 2012 for further details.  

*Arguments*:  
gene_names_by_locus -- dictionary, key = locus_name, value = list of gene names  
interface_sites_input -- file path, mapping locus name to PTMs, their position in a protein interface, the modified residue and the PFAM domain.  
interface_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and interface data  

___

**intf**(mutations, biogrid_IDs)  
Map mutations to PTMs (acetylation, phosphorylation, ubiquitination) at the inferfaces of mutated proteins.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**makeDirTree**(dirs)  

___

**make_bact_file**(uniprot_input, bact_output)  
Extract binding site and active site data from a downloaded UniProt file and write to tab-delimited (tsv) file.  

___

**make_domains_file**(uniprot_input, domains_output)  
Extract domain data from a downloaded UniProt file and write to tab-delimited (tsv) file.  

___

**make_nucleotide_file**(uniprot_input, nucleotide_output)  
Extract nucleotide binding site data from a downloaded UniProt file and write to tab-delimited (tsv) file.  

___

**make_pdb_file**(uniprot_input, pdb_output)  
Extract structural data from a downloaded UniProt file and write to tab-delimited (tsv) file.  

___

**make_ptms_file**(uniprot_input, ptms_output)  
Extract Post-Translational Modification (PTM) data from a downloaded Uniprot file and write to tab-delimited (tsv) file.  

___

**mmap**(mut_prot_input, sites_id_input, mapped_sites_output, summary_output)  
Map the positions of mutations in mutated proteins to binding sites and active sites.  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
sites_id_input -- file path, mapping UniProt ID, ordered locus name, common name, binding/active site, their position and binding partner/activity
mapped_sites_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, binding/active site data  
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**mu_map**(yeastID_input, struct_input, struct_id_output)  
Map UniProt IDs to protein names and structural elements (helices, beta strands, turns).  

*Arguments*:  
yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names  
struct_input -- file path, mapping UniProt ID to PDB ID, structural elements (helices, beta strands, turns) and their position (start and end) in the polypeptide  
struct_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and structural data  

___

**mutation_file**(mut_gene_input, gff_input, d_id_input, mut_prot_output, errors_output)  
Converts a file of mutations at DNA/gene-level to mutations at the protein-level.  

*Arguments*:  
mut_gene_input -- file path, mapping mutated genes (common names) to point mutations on specified chromosomes  
gff_input -- file path, General Feature Format (GFF) file for yeast genome  
d_id_input -- file path, mapping UniProt IDs and gene identifiers to chromosome loci (start, end, strand orientation)  
mut_prot_output -- file path, mapping mutated proteins (common names) to non-synonymous or stop mutations
errors_output -- file path, listing any errors that were raised during data processing  

___

**n_map**(yeastID_input, nucleotide_input, nucleotide_id_output)  
Map UniProt IDs to protein names and nucleotide binding sites.  

*Arguments*:  
yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names  
nucleotide_input -- file path, mapping UniProt ID to nucleotide binding sites, their position (start and end) in the polypeptide and their binding nucleotide  
nucleotide_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and nucleotide binding site data  

___

**nucleo**(mutations, biogrid_IDs)  
Map mutations to nucleotide binding sites.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**nucleotide_map**(mut_prot_input, nucleotide_id_input, mapped_nucleotide_output, summary_output)  
Map the positions of mutations in mutated proteins to nucleotide binding sites.  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
nucleotide_id_input -- file path, mapping UniProt ID, ordered locus name, common name, the start and end position of nucleotide binding sites and the binding nucleotide  
mapped_nucleotide_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, nucleotide binding site data  
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**pTMdata**(uniprot_output)  
Download up-to-date UniProt data file using the following query settings: Species: Saccharomyces cerevisiae, ...  

___

**parse_biogrid**(uniprot_biogrid_input)  
Return dictionary mapping UniProt IDs to BioGrid IDs.  

___

**parse_gene_loci**(d_id_input)  
Parse a file that maps UniProt IDs and gene identifiers to chromosome loci (start, end, strand orientation) into a dictionary.  

Returns the dictionary. Key = gene name. Value = list of start position, end position, strand orientation and chromosome on which the gene is located.  

*Arguments*:  
d_id_input -- file path, mapping UniProt IDs and gene identifiers to chromosome loci (start, end, strand orientation)  

___

**parse_gene_mutations**(mut_gene_input)  
Parse input DNA/gene-level mutation file into a dictionary of mutations.  

Returns the dictionary. Key = mutated gene name. Value = set of tuples, each containing data  
about a mutation: chromosome, position, reference nucleotide (as per reference genome) and
alternative nucleotide.  

*Arguments*:  
mut_gene_input -- file path, mapping mutated genes (common names) to point mutations on specified chromosomes.  

___

**parse_gene_names**(yeastID_input)  
Parse the yeastID_input file (mapping UniProt IDs to gene names) into a dictionary.  

___

**parse_genome**(gff_input)  
Parse General Feature Format (GFF) file for yeast genome into a dictionary.  

Returns the dictionary. Key = chromosome id e.g. chrI. Value = sequence of chromosome  
(from reference genome).  

*Arguments*:  
gff_input -- file path, General Feature Format (GFF) file for yeast genome.  

___

**parse_mutations**(mut_prot_input)  
Parse the mut_prot_input file (mapping mutated proteins to mutation positions) into a dictionary.  

___

**parse_names_by_locus**(gene_names)  
Parse the gene_names into a dictionary (key = locus_name, value = list of gene names).  

*Arguments*:  
gene_names -- dictionary, key = UniProt ID; value = list of (sgd_name, common_name) pairs  

___

**pdb**(mut_prot_input, struct_id_input, mapped_struct_output, summary_output)  
Map the positions of mutations in mutated proteins to structural elements (helices, beta strands, turns).  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
struct_id_input -- file path, mapping UniProt ID, ordered locus name, common name, PDB ID, and the start and end position of each structural element  
mapped_struct_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, nucleotide binding site data  
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**pi**(mutations, biogrid_IDs)  
Map mutations to PTMs known to affect protein-protein interactions (PPIs).  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**pmap**(yeastID_input, ptms_input, ptm_id_output)  
Map UniProt IDs to protein names and post-translational modifications (PTMs).  

*Arguments*:  
yeastID_input -- dictionary, mapping UniProt ID to ordered locus and common (gene) names  
ptms_input -- file path, mapping UniProt ID to a PTM and its position in polypeptide  
ptm_id_output -- file path, mapping UniProt ID, ordered locus name, common name, PTM position and PTM  

___

**ppi**(mut_prot_input, ppi_id_input, mapped_ppi_output, summary_output)  
Map the positions of mutations in mutated proteins to PTMs known to affect protein protein interactions (PPIs).  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
ppi_id_input -- file path, mapping UniProt ID, ordered locus name, common name, PTM position, residue and interacting protein partner  
mapped_interface_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name and PPI data
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**ppi_map**(gene_names_by_locus, ppi_sites_input, ppi_id_output)  
Map UniProt IDs to protein names and PTMs known to affect protein-protein interactions (PPIs).  

See Beltrao et al. Cell 2012 for further details.  

*Arguments*:  
gene_names_by_locus -- dictionary, key = locus_name, value = list of gene names  
ppi_sites_input -- file path, mapping locus name to protein-protein interaction (PPI) data, including the interacting partner (name and locus name), position of PTM, residue affected, and method by which interaction was determined  
ppi_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and PPI data  

___

**preWeb**(uniprot_biogrid_input, mapped_mut_input)  
Map mutations to BioGrid IDs and write to file.  

*Arguments*:  
uniprot_biogrid_input -- dictionary, mapping UniProt IDs to BioGrid IDs  
mapped_mut_input -- file path, mutations mapped to features  

___

**ptm**(mutations, biogrid_IDs)  
Map mutations to PTMs.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**ptm_map**(mut_prot_input, ptm_id_input, mapped_ptms_output, summary_output)  
Map the positions of mutations in mutated proteins to post-translational modifications (PTMs).  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
ptm_id_input -- file path, mapping UniProt ID, ordered locus name, common name, PTM position and PTM  
mapped_ptms_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name, PTM position and PTM  
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**resc**(output_dir)  
Copy pre-downloaded data files (PTMfunc and PTMcode database files) from the ymap package to output dir.  

___

**revcomp**(dna, reverse=True, complement=True)  
Return the reverse complement of a given DNA sequence.  

___

**struc_map**(mutations, biogrid_IDs)  
Map mutations to active/binding sites.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**sum_file_map**(summary_input, final_report_output)  
Generate a final report file.  

___

**translate_dna**(dna)  
Return the protein sequence encoded by the given cDNA sequence.  

___

**uniprot_data**(mutations, biogrid_IDs)  
Call all functions that work with UniProt data.  

___

**web**(biog_input)  
For each BioGrid ID in biog_input, open the corresponding BioGrid database entry in web browser (one tab per entry).  

___

**withP**(mutations, biogrid_IDs)  
Map mutations to PTMs within interacting proteins. See PTMcode2 for further details.  

Perform mapping, GO enrichment on the genes/proteins containing mapped mutations, generate file of associated BioGrid IDs.  
Move files to respective output folders.  

___

**withinPro**(mut_prot_input, within_prot_id_input, mapped_within_prot_output, summary_output)  
Map the positions of mutations in mutated proteins to pairwise combinations of PTMs known/predicted to affect protein-protein interactions (PPIs)????  

*Arguments*:  
mut_prot_input -- dictionary, mapping mutated proteins (common names) to mutated positions in each protein  
ppi_id_input -- file path, mapping UniProt ID, ordered locus name, common name, PTM position, residue and interacting protein partner  
mapped_interface_output -- file path, mapping each mutation in mut_prot_input to UniProt ID, ordered locus name, common name and PPI data
summary_output -- file path, recording a summary of all features to which mutations have been mapped  

___

**withinPro_map**(gene_names_by_locus, within_prot_input, within_prot_id_output)  
Map UniProt IDs to protein names and pairwise combinations of PTMs known/predicted to affect protein-protein interactions (PPIs)????  

See Minguez et al. "PTMcode v2: a resource for functional associations of post-translational modifications within and between proteins." Nucleic Acids Res. 2015 Jan 28; 43 for further details.  

*Arguments*:  
gene_names_by_locus -- dictionary, key = locus_name, value = list of gene names  
within_prot_input -- file path, mapping locus name to protein-protein interaction (PPI) data  
within_prot_id_output -- file path, mapping UniProt ID, ordered locus name, common name, and PPI data  

___

**ymap_genes**()  
Call all functions to analyse DNA-level mutation input file.  

Convert the mutation input file at DNA-level to a mutations file at protein level.  
Call all the functions that map mutations from the protein-level mutation file to features.  

___

**ymap_proteins**()  
returns all the results of all the codes of yMap; starting from proteins level mutation positions  


## 4 - Troubleshooting

**Problem**: During installation, an error occurred stating that "Microsoft Visual C++ 14.0 is required".  
**Solution**: You may need to first install the latest Visual Studio Build Tools from [https://visualstudio.microsoft.com/downloads/](https://visualstudio.microsoft.com/downloads/). Ensure that you are using an up-to-date version of setuptools. However, performing these steps may not resolve the error at all.

**Problem**: yweb fails to locate the directory.  
**Solution**: In python 2.x, the path should be given as "path/to/biog.txt". In python 3.x the path should be given without quotation marks path/to/biog.txt

**Problem**: Only one tab is opened when I use yweb, but the biog.txt maps more than one protein to a BioGrid ID.  
**Solution**: Try first opening your web browser before running yweb.


# Reference
**Ahmed Arslan and Vera van Noort**, *yMap: An automated method to map yeast variants to protein modifications and functional regions*  
Bioinformatics October 22, 2016 doi:10.1093/bioinformatics/btw658

# Contributors

[http://www.biw.kuleuven.be/CSB/](http://www.biw.kuleuven.be/CSB/)

This work is supported by KU Leuven research fund. 
