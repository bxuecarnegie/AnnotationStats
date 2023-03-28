import os
import sys
from joblib import Parallel, delayed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'util'))

from util.functions import output_list_of_species_stats
from util.draw import draw_pie_from_species, go_domain_by_species, completeness_by_species, \
    experimental_chart_by_species
from config import output_folder, model, doe, uniprot, lab_copyright


def write_tsv(gene_list, list_of_mf_seqs, list_of_cc_seqs, list_of_bp_seqs, output_path):
    with open(output_path, 'w') as op:
        op.write('Gene\tMF?\tCC?\tBP?\n')
        for gene in sorted(gene_list):
            mf = 'No'
            if gene in list_of_mf_seqs:
                mf = 'Yes'
            cc = 'No'
            if gene in list_of_cc_seqs:
                cc = 'Yes'
            bp = 'No'
            if gene in list_of_bp_seqs:
                bp = 'Yes'
            op.write('\t'.join([gene, mf, cc, bp]) + '\n')


def run_on_species_class(species_class, annot_ext='tsv'):
    print('Processing Species:', getattr(species_class, 'species_name'))
    print('Reading Gene File...')
    species_class.read_gene_file()
    print('Reading Mapping File...')
    species_class.read_mapping_file()
    print('Reading Annotation File...')
    species_class.read_annot_file()
    species_name = getattr(species_class, 'species_abbr')
    try:
        os.makedirs(os.path.join(output_folder, getattr(species_class, 'species_type')))
    except OSError:
        pass
    species_plot = os.path.join(output_folder, getattr(species_class, 'species_type'), species_name + '.plots.png')
    draw_pie_from_species(species_class, species_plot, title=species_name, extra_txt=lab_copyright)
    exp_gaf = os.path.join(output_folder, getattr(species_class, 'species_type'), species_name + '.exp.' + annot_ext)
    with open(exp_gaf, 'w') as op:
        for gene_id in sorted(getattr(species_class, 'exp_gaf').keys()):
            op.write('\n'.join(sorted(getattr(species_class, 'exp_gaf')[gene_id])) + '\n')
    comp_gaf = os.path.join(output_folder, getattr(species_class, 'species_type'), species_name + '.comp.' + annot_ext)
    with open(comp_gaf, 'w') as op:
        for gene_id in sorted(getattr(species_class, 'comp_gaf').keys()):
            op.write('\n'.join(sorted(getattr(species_class, 'comp_gaf')[gene_id])) + '\n')
    # exp_tsv = os.path.join(output_folder, getattr(species_class, 'species_type'), species_name + '.exp.tsv')
    # write_tsv(getattr(species_class, 'list_of_exp'), getattr(species_class, 'list_of_mf_seqs'),
    #           getattr(species_class, 'list_of_cc_seqs'), getattr(species_class, 'list_of_bp_seqs'), exp_tsv)
    # comp_tsv = os.path.join(output_folder, getattr(species_class, 'species_type'), species_name + '.comp.tsv')
    # write_tsv(getattr(species_class, 'list_of_comp'), getattr(species_class, 'list_of_mf_seqs'),
    #           getattr(species_class, 'list_of_cc_seqs'), getattr(species_class, 'list_of_bp_seqs'), comp_tsv)
    unknown_tsv = os.path.join(output_folder, getattr(species_class, 'species_type'), species_name + '.unknown.tsv')
    with open(unknown_tsv, 'w') as op:
        op.write('Gene\n')
        for gene in sorted(getattr(species_class, 'list_of_unknown')):
            op.write(gene + '\n')


def run_list_of_species_classes(species_class_list, species_type):
    output_file = os.path.join(output_folder, species_type + '.tsv')
    # if species_type == 'doe':
    if species_type:
        for species_class in species_class_list:
            # print(len(species_class.list_of_exp))
            run_on_species_class(species_class)
            # print(len(species_class.list_of_exp))
    else:
        Parallel(n_jobs=5)(delayed(run_on_species_class)(species_class) for species_class in species_class_list)
    output_list_of_species_stats(species_class_list, output_file)
    # domain_plot_path = os.path.join(output_folder, species_type, '.'.join([species_type, "domain.plot.png"]))
    # completeness_plot_path = os.path.join(output_folder, species_type,
    #                                       '.'.join([species_type, "completeness.plot.png"]))
    # experimental_plot_path = os.path.join(output_folder, species_type,
    #                                       '.'.join([species_type, "experimental.plot.png"]))
    # go_domain_by_species(species_class_list, domain_plot_path, extra_txt=lab_copyright)
    # completeness_by_species(species_class_list, completeness_plot_path, extra_txt=lab_copyright)
    # experimental_chart_by_species(species_class_list, experimental_plot_path, extra_txt=lab_copyright)


if __name__ == '__main__':
    # run_list_of_species_classes([model[0]], 'model')
    # run_list_of_species_classes(model, 'model')
    run_list_of_species_classes(doe, 'doe')
    # run_list_of_species_classes(uniprot, 'uniprot')

