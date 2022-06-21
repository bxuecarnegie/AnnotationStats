from definitions import Species

# Model
# yeast Map: in gff ID <-> dbxref
yeast = Species("yeast", "Saccharomyces cerevisiae", "S.cerevisiae", 4932, "model", "gaf",
                "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/S.cerevisiae/sgd.gaf", "Gene Ontology",
                "2022-05-16", "http://release.geneontology.org/2022-05-16/annotations/sgd.gaf.gz", "gff",
                "/Users/bxue/Documents/Carnegie/SourceData/GFFs/S.cerevisiae/S288C_reference_genome_R64-3-1_20210421/"
                "saccharomyces_cerevisiae_R64-3-1_20210421.gff",
                "Saccharomyces Genome Database (sgd)", "R64-3-1",
                "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/"
                "S288C_reference_genome_R64-3-1_20210421.tgz"
                )
# human: gff: ID=gene: <-> map: Entry name
human = Species("human", "Homo sapiens", "H.sapiens", 9606, "model", "gaf",
                "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/H.sapiens/goa_human.gaf", "Gene Ontology",
                "2022-05-16", "http://release.geneontology.org/2022-05-16/annotations/goa_human.gaf.gz", "gff",
                "/Users/bxue/Documents/Carnegie/SourceData/GFFs/H.sapiens/Ensembl_release_106-Apr_2022/"
                "Homo_sapiens.GRCh38.106.genes.gff3",
                "EBI Gene Ontology Annotation Database (goa)", "Ensembl Release 106",
                "http://ftp.ensembl.org/pub/release-106/gff3/homo_sapiens/",
                "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/H.sapiens/"
                "uniprot-yourlist_M2022061492C7BAECDB1C5C413EE0E0348724B68250950FU.tab",
                "UniProt 2022_02", "https://www.uniprot.org/uploadlists/"
                )
# ara: gff: ID=
ara = Species("ara", "Arabidopsis thaliana", "A.thaliana", 3702, "model", "gaf",
              "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/A.thaliana/tair.gaf", "Gene Ontology",
              "2022-05-16", "http://release.geneontology.org/2022-05-16/annotations/tair.gaf.gz", "gff",
              "/Users/bxue/Documents/Carnegie/SourceData/GFFs/A.thaliana/TAIR_2022-05-12/"
              "Araport11_GFF3_genes_transposons.May2022.gff",
              "The Arabidopsis Information Resource (TAIR)", "2022-05-12",
              "https://www-arabidopsis-org.stanford.idm.oclc.org/download_files/Genes/Araport11_genome_release/"
              "Araport11_GFF3_genes_transposons.May2022.gff.gz",
              "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/A.thaliana/"
              "Araport11_TAIRlocusaccessionID_AGI_mapping.txt",
              "2021-02-24",
              "https://www-arabidopsis-org.stanford.idm.oclc.org/download_files/Genes/Araport11_genome_release/"
              "Araport11_TAIRlocusaccessionID_AGI_mapping.txt"
              )
# fly: gff: ID=
fly = Species("fly", "Drosophila melanogaster", "D.melanogaster", 7227, "model", "gaf",
              "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/D.melanogaster/fb.gaf", "Gene Ontology",
              "2022-05-16", "http://release.geneontology.org/2022-05-16/annotations/fb.gaf.gz", "gff",
              "/Users/bxue/Documents/Carnegie/SourceData/GFFs/D.melanogaster/FB2022_02/dmel-genes-r6.45.gff",
              "FlyBase (fb)", "FB2022_02",
              "http://ftp.flybase.net/releases/FB2022_02/dmel_r6.45/gff/dmel-all-r6.45.gff.gz"
              )
# mouse map: in gff ID <-> gene_id
mouse = Species("mouse", "Mus musculus", "M.musculus", 10090, "model", "gaf",
                "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/M.musculus/mgi.gaf", "Gene Ontology",
                "2022-05-16", "http://release.geneontology.org/2022-05-16/annotations/mgi.gaf.gz", "gff",
                "/Users/bxue/Documents/Carnegie/SourceData/GFFs/M.musculus/GRCm39/MGI.gff3",
                "Mouse Genome Informatics (mgi)", "6.19",
                "http://www.informatics.jax.org/downloads/mgigff3/archive/monthly/MGI.202206.gff3.gz"
                )
# zebrafish: gff: ID=
zebrafish = Species("zebrafish", "Danio rerio", "D.rerio", 7955, "model", "gaf",
                    "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/D.rerio/zfin.gaf", "Gene Ontology",
                    "2022-05-16", "http://release.geneontology.org/2022-05-16/annotations/zfin.gaf.gz", "gff",
                    "/Users/bxue/Documents/Carnegie/SourceData/GFFs/D.rerio/ZFIN_13-Jun-2022/zfin_genes.gff3",
                    "Zebrafish Information Network (zfin)", "13 Jun 2022",
                    "https://zfin.org/downloads/archive/2022.06.13/zfin_genes.gff3"
                    )
# roundworm: gff: ID=Gene:
roundworm = Species("roundworm", "Caenorhabditis elegans", "C.elegans", 6239, "model", "gaf",
                    "/Users/bxue/Documents/Carnegie/SourceData/GO/2022-05-16/C.elegans/wb.gaf", "Gene Ontology",
                    "2022-05-16", "http://release.geneontology.org/2022-05-16/annotations/wb.gaf.gz", "gff",
                    "/Users/bxue/Documents/Carnegie/SourceData/GFFs/C.elegans/WormBase_WS283/"
                    "c_elegans.PRJNA13758.WS283.genes.gff3",
                    "WormBase database of nematode biology (wb)", "WS283",
                    "https://downloads.wormbase.org/releases/WS283/species/c_elegans/PRJNA13758/"
                    "c_elegans.PRJNA13758.WS283.annotations.gff3.gz",
                    )

# JGI
poplar = Species("poplar", "Populus trichocarpa", "P.trichocarpa", 3694, "doe", "annotation_info",
                 "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Ptrichocarpa_v4.1/annotation/"
                 "Ptrichocarpa_533_v4.1.annotation_info.txt", "Phytozome", "13",
                 "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
                 "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Ptrichocarpa_v4.1/annotation/"
                 "Ptrichocarpa_533_v4.1.gene.gff3", "Phytozome", "13",
                 "https://data.jgi.doe.gov/refine-download/phytozome"
                 )
brachypodium = Species("brachypodium", "Brachypodium distachyon", "B.distachyon", 15368, "doe", "annotation_info",
                       "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Bdistachyon_v3.2/annotation/"
                       "Bdistachyon_556_v3.2.annotation_info.txt", "Phytozome", "13",
                       "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
                       "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Bdistachyon_v3.2/annotation/"
                       "Bdistachyon_556_v3.2.gene.gff3", "Phytozome", "13",
                       "https://data.jgi.doe.gov/refine-download/phytozome"
                       )
phallii = Species("phallii", "Panicum hallii", "P.hallii", 1504633, "doe", "annotation_info",
                  "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Phallii_v3.2/annotation/"
                  "Phallii_590_v3.2.annotation_info.txt", "Phytozome", "13",
                  "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
                  "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Phallii_v3.2/annotation/"
                  "Phallii_590_v3.2.gene.gff3", "Phytozome", "13",
                  "https://data.jgi.doe.gov/refine-download/phytozome"
                  )
soy = Species("soy", "Glycine max", "G.max", 3847, "doe", "annotation_info",
              "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Gmax_Wm82.a4.v1/annotation/"
              "Gmax_508_Wm82.a4.v1.annotation_info.txt", "Phytozome", "13",
              "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
              "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Gmax_Wm82.a4.v1/annotation/"
              "Gmax_508_Wm82.a4.v1.gene.gff3", "Phytozome", "13",
              "https://data.jgi.doe.gov/refine-download/phytozome"
              )
chlamy = Species("chlamy", "Chlamydomonas reinhardtii", "C.reinhardtii", 3055, "doe", "annotation_info",
                 "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Creinhardtii_v5.6/annotation/"
                 "Creinhardtii_281_v5.6.annotation_info.txt", "Phytozome", "13",
                 "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
                 "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Creinhardtii_v5.6/annotation/"
                 "Creinhardtii_281_v5.6.gene.gff3", "Phytozome", "13",
                 "https://data.jgi.doe.gov/refine-download/phytozome"
                 )
setaria = Species("setaria", "Setaria italica", "S.italica", 4555, "doe", "annotation_info",
                  "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV12/Sitalica_312_v2.2/annotation/"
                  "Sitalica_312_v2.2.annotation_info.txt", "Phytozome", "13",
                  "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
                  "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV12/Sitalica_312_v2.2/annotation/"
                  "Sitalica_312_v2.2.gene.gff3", "Phytozome", "13",
                  "https://data.jgi.doe.gov/refine-download/phytozome"
                  )
sorghumbicolor = Species("sorghumbicolor", "Sorghum bicolor", "S.bicolor", 4558, "doe", "annotation_info",
                         "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV12/Sbicolor_454_v3.1.1/"
                         "annotation/Sbicolor_454_v3.1.1.annotation_info.txt", "Phytozome", "13",
                         "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
                         "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV12/Sbicolor_454_v3.1.1/"
                         "annotation/Sbicolor_454_v3.1.1.gene.gff3", "Phytozome", "13",
                         "https://data.jgi.doe.gov/refine-download/phytozome"
                         )
switchgrass = Species("switchgrass", "Panicum virgatum", "P.virgatum", 38727, "doe", "annotation_info",
                      "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Pvirgatum_v5.1/annotation/"
                      "Pvirgatum_516_v5.1.annotation_info.txt", "Phytozome", "13",
                      "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
                      "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV13/Pvirgatum_v5.1/annotation/"
                      "Pvirgatum_516_v5.1.gene.gff3", "Phytozome", "13",
                      "https://data.jgi.doe.gov/refine-download/phytozome"
                      )
moss = Species("moss", "Physcomitrium patens", "P.patens", 3218, "doe", "annotation_info",
               "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV12/Ppatens_318_v3.3/annotation/"
               "Ppatens_318_v3.3.annotation_info.txt", "Phytozome", "13",
               "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
               "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV12/Ppatens_318_v3.3/annotation/"
               "Ppatens_318_v3.3.gene.gff3", "Phytozome", "13",
               "https://data.jgi.doe.gov/refine-download/phytozome"
               )
msinensis = Species("msinensis", "Miscanthus sinensis", "M.sinensis", 62337, "doe", "annotation_info",
                    "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV12/Msinensis_497_v7.1/annotation/"
                    "Msinensis_497_v7.1.annotation_info.txt", "Phytozome", "13",
                    "https://data.jgi.doe.gov/refine-download/phytozome", "gff",
                    "/Users/bxue/Documents/Carnegie/SourceData/Phytozome/PhytozomeV12/Msinensis_497_v7.1/annotation/"
                    "Msinensis_497_v7.1.gene.gff3", "Phytozome", "13",
                    "https://data.jgi.doe.gov/refine-download/phytozome"
                    )

# Uniprot
breadwheat = Species("breadwheat", "Triticum aestivum", "T.aestivum", 4565, "uniprot", "GAF",
                     "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/T.aestivum/"
                     "QuickGO-annotations-1655794633092-20220621.gaf", "Uniprot GOA", "2022_02",
                     "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                     "uniprot-proteome.tab",
                     "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/T.aestivum/"
                     "uniprot-proteome_UP000019116.tab", "UniProt", "2022_02",
                     "https://www.uniprot.org/proteomes/?query=taxonomy:4565"
                     )
tax_3635 = Species("TAX-3635", "Gossypium hirsutum", "G.hirsutum", 3635, "uniprot", "GAF",
                   "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/G.hirsutum/"
                   "QuickGO-annotations-1655795547755-20220621.gaf", "Uniprot GOA", "2022_02",
                   "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                   "uniprot-proteome.tab",
                   "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/G.hirsutum/"
                   "uniprot-proteome_UP000189702.tab", "UniProt", "2022_02",
                   "https://www.uniprot.org/proteomes/?query=taxonomy:3635"
                   )
corn = Species("corn", "Zea mays", "Z.mays", 4577, "uniprot", "GAF",
               "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/Z.mays/"
               "QuickGO-annotations-1655797575418-20220621.gaf", "Uniprot GOA", "2022_02",
               "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
               "uniprot-proteome.tab",
               "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/Z.mays/"
               "uniprot-proteome_UP000007305.tab", "UniProt", "2022_02",
               "https://www.uniprot.org/proteomes/?query=taxonomy:4577"
               )
ntabacum_tn90 = Species("ntabacum_tn90", "Nicotiana tabacum", "N.tabacum", 4097, "uniprot", "GAF",
                        "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/N.tabacum/"
                        "QuickGO-annotations-1655796483482-20220621.gaf", "Uniprot GOA", "2022_02",
                        "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                        "uniprot-proteome.tab",
                        "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/N.tabacum/"
                        "uniprot-proteome_UP000084051.tab", "UniProt", "2022_02",
                        "https://www.uniprot.org/proteomes/?query=taxonomy:4097"
                        )
soleracea = Species("soleracea", "Spinacia oleracea", "S.oleracea", 3562, "uniprot", "GAF",
                    "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/S.oleracea/"
                    "QuickGO-annotations-1655797028570-20220621.gaf", "Uniprot GOA", "2022_02",
                    "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                    "uniprot-proteome.tab",
                    "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/S.oleracea/"
                    "uniprot-proteome_UP000054095.tab", "UniProt", "2022_02",
                    "https://www.uniprot.org/proteomes/?query=taxonomy:3562"
                    )
mtruncatula = Species("mtruncatula", "Medicago truncatula", "M.truncatula", 3880, "uniprot", "GAF",
                      "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/M.truncatula/"
                      "QuickGO-annotations-1655795989299-20220621.gaf", "Uniprot GOA", "2022_02",
                      "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                      "uniprot-proteome.tab",
                      "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/M.truncatula/"
                      "uniprot-proteome_UP000002051.tab", "UniProt", "2022_02",
                      "https://www.uniprot.org/proteomes/?query=taxonomy:3880"
                      )
tax_3469 = Species("TAX-3469", "Papaver somniferum", "P.somniferum", 3469, "uniprot", "GAF",
                   "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/P.somniferum/"
                   "QuickGO-annotations-1655793536290-20220621.gaf", "Uniprot GOA", "2022_02",
                   "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                   "uniprot-proteome.tab",
                   "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/P.somniferum/"
                   "uniprot-proteome_UP000316621.tab", "UniProt", "2022_02",
                   "https://www.uniprot.org/proteomes/?query=taxonomy:3469"
                   )
castorbean = Species("castorbean", "Ricinus communis", "R.communis", 3988, "uniprot", "GAF",
                     "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/R.communis/"
                     "QuickGO-annotations-1655796118330-20220621.gaf", "Uniprot GOA", "2022_02",
                     "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                     "uniprot-proteome.tab",
                     "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/R.communis/"
                     "uniprot-proteome_UP000008311.tab", "UniProt", "2022_02",
                     "https://www.uniprot.org/proteomes/?query=taxonomy:3988"
                     )
potato = Species("potato", "Solanum tuberosum", "S.tuberosum", 4113, "uniprot", "GAF",
                 "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/S.tuberosum/"
                 "QuickGO-annotations-1655796691613-20220621.gaf", "Uniprot GOA", "2022_02",
                 "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                 "uniprot-proteome.tab",
                 "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/S.tuberosum/"
                 "uniprot-proteome_UP000011115.tab", "UniProt", "2022_02",
                 "https://www.uniprot.org/proteomes/?query=taxonomy:4113"
                 )
oryza = Species("oryza", "Oryza sativa Japonica Group", "O.sativa", 39947, "uniprot", "GAF",
                "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/O.sativa/"
                "QuickGO-annotations-1655798173033-20220621.gaf", "Uniprot GOA", "2022_02",
                "http://release.geneontology.org/2022-05-16/annotations/goa_uniprot_all.gaf.gz",
                "uniprot-proteome.tab",
                "/Users/bxue/Documents/Carnegie/SourceData/UniProt/2022_02/O.sativa/"
                "uniprot-proteome_UP000059680.tab", "UniProt", "2022_02",
                "https://www.uniprot.org/proteomes/?query=taxonomy:39947"
                )
