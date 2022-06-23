from goatools.obo_parser import GODag

from util.species import Species

output_folder = './output'

# GO link: http://release.geneontology.org/2022-05-16/ontology/go-basic.obo
go_obo_path = '~/Downloads/go-basic.obo'
g = GODag(go_obo_path, prt=None)

# Date
date = '06/13/22'
lab_copyright = 'Generated on ' + date + '\n© 2022 Rhee-Lab'

# Model
model = []
ara = Species("ara", "Arabidopsis thaliana", "A.thaliana", 3702, "model", "gaf",
              "~/Downloads/tair.gaf", "Gene Ontology",
              "2022-05-16", "http://release.geneontology.org/2022-05-16/annotations/tair.gaf.gz", "gff",
              "~/Downloads/Araport11_GFF3_genes_transposons.May2022.gff",
              "The Arabidopsis Information Resource (TAIR)", "2022-05-12",
              "https://www-arabidopsis-org.stanford.idm.oclc.org/download_files/Genes/Araport11_genome_release/"
              "Araport11_GFF3_genes_transposons.May2022.gff.gz",
              "~/Downloads/TAIR2UniprotMapping.txt",
              "2021-07-05",
              "https://www-arabidopsis-org.stanford.idm.oclc.org/download_files/Proteins/Id_conversions/"
              "TAIR2UniprotMapping.txt.gz"
              )
ara.__setattr__('date', date)
model.append(ara)
# JGI
jgi = []

# Uniprot
uniprot = []