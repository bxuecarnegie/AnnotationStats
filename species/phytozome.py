import os
import re

from Bio.SeqIO import parse
from goatools.obo_parser import GODag, GraphEngines

from definitions import aspect_list, date
from util.draw import draw_pie_from_input, go_domain_by_species, completeness_by_species
from util.functions import write_tsv


go_obo_path = '/Users/bxue/Documents/Carnegie/SourceData/GO/go-basic_21-11-08.obo'
g = GODag(go_obo_path, prt=None)

namespace_to_aspect = {
    "molecular_function": "F",
    "biological_process": "P",
    "cellular_component": "C"
}


def count_num_of_seq(input_fasta):
    count = 0
    for _ in parse(input_fasta, "fasta"):
        count += 1
    return count


def read_phytozome_annotation_info(annotation_info_path, with_exp=True):
    annotation_dict = {}
    for aspect in sorted(aspect_list.keys()):
        aspect_dict = {
            "pred": set(),
            "exp": set()
        }
        if with_exp is True:
            aspect_dict.setdefault('exp', set())
        annotation_dict.setdefault(aspect, aspect_dict)
    unfound_go_terms = set()
    with open(annotation_info_path, 'r') as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n')
            info = line.split('\t')
            try:
                protein_id = info[3]
                go_id_list = [go.strip() for go in re.split(r'[, ]', info[9]) if go.strip() != '']
                if len(go_id_list) > 0:
                    for go_id in go_id_list:
                        res = g.query_term(go_id)
                        if res is not None:
                            try:
                                aspect = namespace_to_aspect[res.namespace]
                                annotation_dict[aspect]['pred'].add(protein_id)
                            except KeyError:
                                print("Aspect not found", res)
                                continue
                        else:
                            unfound_go_terms.add(go_id)
            except IndexError:
                print('Index Error', info)
                continue
    print("Unfound GO", unfound_go_terms)
    return annotation_dict


def find_files_of_phytozome(input_folder):
    annotatino_info_file = None
    fasta_file = None
    for f in os.listdir(input_folder):
        file_path = os.path.join(input_folder, f)
        if f.endswith('.annotation_info.txt') and os.path.isfile(file_path):
            annotatino_info_file = file_path
        elif (f.endswith('fa') or f.endswith('fasta')) and os.path.isfile(file_path):
            fasta_file = file_path
    return annotatino_info_file, fasta_file


def process_phytozome_folder(input_folder, species_name, plot_path=None, with_exp=False, plot_venn=False, date=date):
    annotatino_info_file, fasta_file = find_files_of_phytozome(input_folder)
    if None not in (annotatino_info_file, fasta_file):
        annotation_dict = read_phytozome_annotation_info(annotatino_info_file)
        total_protein_num = count_num_of_seq(fasta_file)
        if with_exp is False:
            total_protein_num, aspect_pie, species_gene_count = \
                draw_pie_from_input(total_protein_num, annotation_dict, plot_path=plot_path, title=species_name,
                                    with_exp=with_exp, plot_venn=plot_venn, extra_txt=date)
        else:
            total_protein_num, aspect_pie, species_gene_count, _, _ = \
                draw_pie_from_input(total_protein_num, annotation_dict, plot_path=plot_path, title=species_name,
                                    with_exp=with_exp, plot_venn=plot_venn, extra_txt=date)
        return total_protein_num, aspect_pie, species_gene_count


def plot_phytozome_list(list_of_species_name, root_folder_path, plot_folder):
    domain_plot_path = os.path.join(plot_folder, "domain.plot.png")
    completeness_plot_path = os.path.join(plot_folder, "completeness.plot.png")
    tsv_path = os.path.join(plot_folder, "doe_stats.tsv")

    list_of_total_gene = []
    list_of_aspect_pie = []
    list_of_species_gene_count = []
    for species_name in list_of_species_name:
        print(species_name)
        species_folder = os.path.join(root_folder_path, species_name)
        species_plot_path = os.path.join(plot_folder, species_name + '.plots.png')
        if os.path.isdir(species_folder):
            out = process_phytozome_folder(species_folder, species_name, plot_path=species_plot_path, with_exp=True)
            if out is not None and len(out) > 1:
                total_gene_num, aspect_pie, species_gene_count = out
                list_of_total_gene.append(total_gene_num)
                list_of_aspect_pie.append(aspect_pie)
                list_of_species_gene_count.append(species_gene_count)

    go_domain_by_species(list_of_aspect_pie, list_of_species_name, domain_plot_path, extra_txt=date)
    completeness_by_species(list_of_species_gene_count, list_of_total_gene, list_of_species_name,
                            completeness_plot_path, extra_txt=date)
    print(list_of_aspect_pie)
    write_tsv(list_of_species_name, list_of_total_gene, list_of_species_gene_count, tsv_path, with_exp=False,
              list_of_aspect_pie=list_of_aspect_pie)


if __name__ == '__main__':
    list_of_species_name = ["P.trichocarpa", "B.distachyon", "P.hallii", "G.max", "C.reinhardtii", "S.italica", "S.bicolor", "P.virgatum", "P.patens", "M.sinensis"]
    # list_of_species_name = ["P.trichocarpa", "B.distachyon"]
    root_folder_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/DOE_Flagship/'
    plot_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/DOE_Flagship/plots_protein'
    plot_phytozome_list(list_of_species_name, root_folder_path, plot_folder)