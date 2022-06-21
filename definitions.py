aspect_list = {
    'F': "Molecular function",
    'P': "Biological process",
    'C': "Cellular component"
}
aspect_order = ["Molecular function", "Biological process", "Cellular component"]

experimental_evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
# computational_analysis_evidence_codes = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA']

annot_type_all = ['Experimental evidence', 'Predicted', 'Unknown']
annot_type_no_exp = ['Predicted', 'Unknown']
# colors = ['#769f4f', '#2c7fb8', '#000000']
colors_all = ['#70ad47', '#3365ff', '#000000']
colors_no_exp = ['#3365ff', '#000000']

# Date
date = 'Generated on 06/20/22\n© 2022 Rhee-Lab'


# Species Class
class Species(object):
    """Object for running all classifiers
    """
    def __init__(self, species_id, species_name, species_abbr, species_taxon, species_type,
                 annot_type, annot_file, annot_source, annot_version, annot_link,
                 gene_type, gene_file, gene_source, gene_version, gene_link,
                 mapping_file=None, mapping_version=None, mapping_link=None):
        self.species_id = species_id
        self.species_name = species_name
        self.species_abbr = species_abbr
        self.species_taxon = species_taxon
        # Model, JGI, Uniprot
        self.species_type = species_type
        # GAF, Phytozome
        self.annot_type = annot_type
        self.annot_file = annot_file
        self.annot_source = annot_source
        self.annot_version = annot_version
        self.annot_link = annot_link
        # GFF, Gene, UniProt
        self.gene_type = gene_type
        self.gene_file = gene_file
        self.gene_source = gene_source
        self.gene_version = gene_version
        self.gene_link = gene_link
        # Mapping between annot & gene
        self.mapping_file = mapping_file
        self.mapping_version = mapping_version
        self.mapping_link = mapping_link

        # Sequences
        self.num_of_seqs = 0
        self.num_of_exp = 0
        self.num_of_comp = 0
        self.num_of_unknown = 0

        # Domain
        self.num_of_mf_seqs = 0
        self.num_of_mf_exp = 0
        self.num_of_mf_comp = 0
        self.num_of_mf_unknown = 0
        self.num_of_cc_seqs = 0
        self.num_of_cc_exp = 0
        self.num_of_cc_comp = 0
        self.num_of_cc_unknown = 0
        self.num_of_bp_seqs = 0
        self.num_of_bp_exp = 0
        self.num_of_bp_comp = 0
        self.num_of_bp_unknown = 0
