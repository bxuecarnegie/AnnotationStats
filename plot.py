import os
import sys

from util.draw import go_domain_by_species, completeness_by_species, experimental_chart_by_species, \
    draw_pie_from_species

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'util'))

from util.functions import output_list_of_species_class, delim_parser
from config import output_folder, model, doe, uniprot, lab_copyright


def add_stats_to_species(species_class, stat_table):
    with open(stat_table, 'r') as fp:
        for line_num, info in \
                delim_parser(fp, indices=[0, 9, 10, 11, 12, 13, 15, 17, 19, 21, 23, 25, 27, 29], ignore=['ID']):
            if info is not None:
                if getattr(species_class, 'species_id') == info[0]:
                    setattr(species_class, 'list_of_seqs', [''] * int(info[1]))
                    setattr(species_class, 'list_of_exp', [''] * int(info[2]))
                    setattr(species_class, 'list_of_comp', [''] * int(info[3]))
                    setattr(species_class, 'list_of_unknown', [''] * int(info[4]))
                    setattr(species_class, 'list_of_mf_exp', [''] * int(info[5]))
                    setattr(species_class, 'list_of_mf_comp', [''] * int(info[6]))
                    setattr(species_class, 'list_of_mf_unknown', [''] * int(info[7]))
                    setattr(species_class, 'list_of_bp_exp', [''] * int(info[8]))
                    setattr(species_class, 'list_of_bp_comp', [''] * int(info[9]))
                    setattr(species_class, 'list_of_bp_unknown', [''] * int(info[10]))
                    setattr(species_class, 'list_of_cc_exp', [''] * int(info[11]))
                    setattr(species_class, 'list_of_cc_comp', [''] * int(info[12]))
                    setattr(species_class, 'list_of_cc_unknown', [''] * int(info[13]))


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def plot_type(species_tsv, species_type, species_class_list, output, max_col=8):
    for species_class in species_class_list:
        species_name = getattr(species_class, 'species_abbr')
        print(species_name)
        species_plot = os.path.join(output, species_type, species_name + '.plots.png')
        add_stats_to_species(species_class, species_tsv)
        # draw_pie_from_species(species_class, species_plot, title=species_name, extra_txt=lab_copyright)
    for idx, species_class_chunks in enumerate(divide_chunks(species_class_list, max_col)):
        if len(species_class_chunks) == 0:
            continue
        domain_plot_path = os.path.join(output, species_type, species_type + '.' + str(idx) + ".domain.plot.png")
        completeness_plot_path = os.path.join(output, species_type, species_type + '.' + str(idx) + ".completeness.plot.png")
        experimental_plot_path = os.path.join(output, species_type, species_type + '.' + str(idx) + ".experimental.plot.png")
        go_domain_by_species(species_class_chunks, domain_plot_path, extra_txt=lab_copyright)
        completeness_by_species(species_class_chunks, completeness_plot_path, extra_txt=lab_copyright)
        experimental_chart_by_species(species_class_chunks, experimental_plot_path, extra_txt=lab_copyright)


if __name__ == '__main__':
    model_tsv = 'output/model.tsv'
    doe_tsv = 'output/doe.tsv'
    uniprot_tsv = 'output/uniprot.tsv'

    plot_type(model_tsv, 'model', model, output_folder)
    plot_type(doe_tsv, 'doe', doe, output_folder, max_col=10)
    plot_type(uniprot_tsv, 'uniprot', uniprot, output_folder, max_col=10)

