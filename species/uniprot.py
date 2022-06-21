import os

from definitions import date
from util.draw import draw_pie_from_input, go_domain_by_species, completeness_by_species, experimental_chart_by_species
from util.functions import read_tab_file_to_kv, read_gaf, write_tsv


def find_files_of_uniprot(input_folder):
    gaf_file = None
    uniprot_tab = None
    for f in os.listdir(input_folder):
        file_path = os.path.join(input_folder, f)
        if f.endswith('gaf') and os.path.isfile(file_path):
            gaf_file = file_path
        elif f.startswith('uniprot-') and f.endswith('.tab') and os.path.isfile(file_path) and 'taxonomy' in f:
            uniprot_tab = file_path
    return gaf_file, uniprot_tab


def read_uniprot_taxon_tab(uniprot_taxon_tab_path):
    uniprot_entry_to_name = read_tab_file_to_kv(uniprot_taxon_tab_path, key_idx=0, val_idx=1, one_val=True)
    return uniprot_entry_to_name


def process_uniprot_folder(input_folder, species_name, plot_path=None):
    gaf_file, uniprot_tab = find_files_of_uniprot(input_folder)

    if None not in (gaf_file, uniprot_tab):
        gaf_dict = read_gaf(gaf_file, db_obj_id_only=True)
        uniprot_entry_to_name = read_uniprot_taxon_tab(uniprot_tab)
        total_protein_num = len(uniprot_entry_to_name)
        total_protein_num, aspect_pie, species_gene_count, exp_aspect_sets, exp_aspect_labels = \
            draw_pie_from_input(total_protein_num, gaf_dict, plot_path=plot_path, title=species_name, extra_txt=date)
        return total_protein_num, aspect_pie, species_gene_count, exp_aspect_sets, exp_aspect_labels


def plot_uniprot_list(list_of_species_name, root_folder_path, plot_folder):
    domain_plot_path = os.path.join(plot_folder, "domain.plot.png")
    completeness_plot_path = os.path.join(plot_folder, "completeness.plot.png")
    experimental_plot_path = os.path.join(plot_folder, "experimental.plot.png")
    tsv_path = os.path.join(plot_folder, "uniprot_stats.tsv")

    list_of_total_protein = []
    list_of_aspect_pie = []
    list_of_species_protein_count = []
    list_of_exp_aspect_sets = []
    list_of_exp_aspect_labels = []
    for species_name in list_of_species_name:
        print(species_name)
        species_folder = os.path.join(root_folder_path, species_name)
        species_plot_path = os.path.join(plot_folder, species_name + '.plots.png')
        if os.path.isdir(species_folder):
            out = process_uniprot_folder(species_folder, species_name, plot_path=species_plot_path)
            if out is not None and len(out) > 1:
                total_protein_num, aspect_pie, species_protein_count, exp_aspect_sets, exp_aspect_labels = out
                list_of_total_protein.append(total_protein_num)
                list_of_aspect_pie.append(aspect_pie)
                list_of_species_protein_count.append(species_protein_count)
                list_of_exp_aspect_sets.append(exp_aspect_sets)
                list_of_exp_aspect_labels.append(exp_aspect_labels)
            else:
                print('error', species_name)
    go_domain_by_species(list_of_aspect_pie, list_of_species_name, domain_plot_path, extra_txt=date)
    completeness_by_species(list_of_species_protein_count, list_of_total_protein, list_of_species_name,
                            completeness_plot_path, extra_txt=date)
    experimental_chart_by_species(list_of_exp_aspect_sets, list_of_exp_aspect_labels, list_of_species_name,
                                  experimental_plot_path, extra_txt=date)
    write_tsv(list_of_species_name, list_of_total_protein, list_of_species_protein_count, tsv_path, list_of_aspect_pie=list_of_aspect_pie)


def plot_main():
    list_of_species_name = ["G.hirsutum", "M.truncatula", "N.tabacum", "O.sativa", "S.oleracea", "S.tuberosum", "T.aestivum", "Z.mays"]
    root_folder_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/NCBI_Viridiplantae/top_species'
    plot_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/NCBI_Viridiplantae/top_species/plots_protein'
    plot_uniprot_list(list_of_species_name, root_folder_path, plot_folder)


if __name__ == '__main__':
    list_of_species_name = ["T.aestivum", "G.hirsutum", "Z.mays", "N.tabacum", "S.oleracea", "M.truncatula",
                            "P.somniferum", "R.communis", "S.tuberosum", "O.sativa"]
    root_folder_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/NCBI_Viridiplantae/top_species'
    plot_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/NCBI_Viridiplantae/top_species/plots_protein'
    plot_uniprot_list(list_of_species_name, root_folder_path, plot_folder)