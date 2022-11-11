import os.path
import re

from functions import delim_parser, read_gff, read_tair2uniprot, read_tab_file_to_kv, read_gaf, \
    read_phytozome_annotation_info, flatten_list


# Species Class
class Species(object):
    """
    Class for a Species
    """
    def __init__(self, species_id, species_name, species_abbr, species_taxon, species_type,
                 annot_type, annot_file, annot_source, annot_version, annot_link,
                 gene_type, gene_file, gene_source, gene_version, gene_link,
                 mapping_file=None, mapping_version=None, mapping_link=None, gff_bio_type=None):
        self.date = ''
        self.go_dag = None

        self.species_id = species_id
        self.species_name = species_name
        self.species_abbr = species_abbr
        self.species_taxon = species_taxon
        # Model, JGI, Uniprot
        self.species_type = species_type
        # Annotations: GAF, Phytozome
        self.annot_type = annot_type
        if not os.path.isfile(annot_file):
            print("Annotation File not found", annot_file)
            raise SystemError
        self.annot_file = annot_file
        self.annot_source = annot_source
        self.annot_version = annot_version
        self.annot_link = annot_link

        # Genes: gff, gene, uniprot
        self.gene_type = gene_type
        self.gene_file = gene_file
        if not os.path.isfile(gene_file):
            print("Gene File not found", gene_file)
            raise SystemError
        self.gene_source = gene_source
        self.gene_version = gene_version
        self.gene_link = gene_link
        # GFF settings
        self.feature_type = None
        self.attribute_type = None

        # Mapping between annot & gene
        self.mapping_file = mapping_file
        if mapping_file:
            if not os.path.isfile(mapping_file):
                print("Mapping File not found", mapping_file)
                raise SystemError
        self.mapping_version = mapping_version
        self.mapping_link = mapping_link
        self.protein_to_gene_map = None
        self.gff_bio_type = gff_bio_type

        # Sequences
        self.list_of_seqs = []
        self.list_of_exp = []
        self.list_of_comp = []
        self.list_of_unknown = []

        # Domain
        self.list_of_mf_seqs = []
        self.list_of_mf_exp = []
        self.list_of_mf_comp = []
        self.list_of_mf_unknown = []
        self.list_of_cc_seqs = []
        self.list_of_cc_exp = []
        self.list_of_cc_comp = []
        self.list_of_cc_unknown = []
        self.list_of_bp_seqs = []
        self.list_of_bp_exp = []
        self.list_of_bp_comp = []
        self.list_of_bp_unknown = []

    def read_gene_file(self):
        if self.gene_type.lower() == 'gff':
            seqs_from_gff = read_gff(self.gene_file, feature_type=self.feature_type, attribute_type=self.attribute_type,
                                     bio_type=self.gff_bio_type)
            self.list_of_seqs = sorted(set(flatten_list(seqs_from_gff)))
        elif self.gene_type.lower() == 'uniprot-proteome.tab':
            seqs_from_tab = set()
            protein_to_genes = {}
            with open(self.gene_file, 'r') as fp:
                for line_idx, line in enumerate(fp):
                    if line_idx == 0 and 'Entry\t' in line:
                        info = [i.strip() for i in line.rstrip('\n').split('\t')]
                        try:
                            entry_idx = info.index('Entry')
                        except ValueError:
                            entry_idx = 0
                        try:
                            gene_names_idx = info.index('Gene names')
                        except ValueError:
                            gene_names_idx = 9
            with open(self.gene_file, 'r') as fp:
                for line_num, info in delim_parser(fp, ignore=['Entry'], indices=[entry_idx, gene_names_idx]):
                    if info is None:
                        continue
                    entry, gene_names = info
                    if gene_names == '':
                        genes = {entry}
                    else:
                        genes = set(re.split(r'[;\s]+', gene_names))
                    # Get first of sorted all possible gene names
                    gene = sorted(genes)[0]
                    seqs_from_tab.add(gene)
                    protein_to_genes.setdefault(entry, [gene])
            if len(protein_to_genes) == 0:
                print('No proteins found', self.gene_file)
                raise SystemError
            self.list_of_seqs = sorted(seqs_from_tab)
            self.protein_to_gene_map = protein_to_genes
        else:
            print('Gene File type unknown', self.gene_type)
            raise SystemError

    def read_mapping_file(self):
        if self.mapping_file is not None and self.protein_to_gene_map is None:
            # Ara specific
            if self.species_taxon == 3702:
                if len(self.list_of_seqs) == 0:
                    print('No genes found, run \'read_gene_file()\' first.')
                else:
                    self.protein_to_gene_map = read_tair2uniprot(self.mapping_file, self.list_of_seqs)
            else:
                self.protein_to_gene_map = read_tab_file_to_kv(self.mapping_file, key_idx=1, val_idx=0, val_sep=',')

    def read_annot_file(self):
        if len(self.list_of_seqs) == 0:
            print('No genes found, run \'read_gene_file()\' first.')
            raise SystemError
        if self.mapping_file is not None and len(self.protein_to_gene_map) == 0:
            print('Mapping file exists, run \'read_mapping_file()\' first.')
            raise SystemError
        if self.annot_type == 'gaf':
            self.list_of_exp, self.list_of_comp, self.list_of_mf_seqs, self.list_of_mf_exp, self.list_of_mf_comp, \
                self.list_of_cc_seqs, self.list_of_cc_exp, self.list_of_cc_comp, \
                self.list_of_bp_seqs, self.list_of_bp_exp, self.list_of_bp_comp =\
                read_gaf(self.annot_file, protein_to_gene=self.protein_to_gene_map, id_list=self.list_of_seqs)

        elif self.annot_type == 'annotation_info':
            if self.go_dag is None:
                print("GO dag not generated, run 'generate_go_dag() first.'")
            self.list_of_comp, self.list_of_mf_seqs, self.list_of_mf_comp, self.list_of_cc_seqs, self.list_of_cc_comp, \
                self.list_of_bp_seqs, self.list_of_bp_comp = \
                read_phytozome_annotation_info(self.annot_file, self.go_dag)
        else:
            print("Annotation File type unknown", self.annot_type)
            raise SystemError

        self.list_of_unknown = sorted(set(self.list_of_seqs) - set(self.list_of_exp) - set(self.list_of_comp))
        self.list_of_mf_unknown = \
            sorted(set(self.list_of_seqs) - set(self.list_of_mf_exp) - set(self.list_of_mf_comp))
        self.list_of_cc_unknown = \
            sorted(set(self.list_of_seqs) - set(self.list_of_cc_exp) - set(self.list_of_cc_comp))
        self.list_of_bp_unknown = \
            sorted(set(self.list_of_seqs) - set(self.list_of_bp_exp) - set(self.list_of_bp_comp))



