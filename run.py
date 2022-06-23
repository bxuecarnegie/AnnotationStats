import os
import sys
from joblib import Parallel, delayed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'util'))

from util.functions import output_list_of_species_class
from util.draw import draw_pie_from_species, go_domain_by_species, completeness_by_species, \
    experimental_chart_by_species
from config import output_folder, model, jgi, uniprot, lab_copyright


def run_on_species_class(species_class):
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


def run_list_of_species_classes(species_class_list, species_type):
    output_file = os.path.join(output_folder, species_type + '.tsv')

    # Parallel(n_jobs=5)(delayed(run_on_species_class)(species_class) for species_class in species_class_list)
    for species_class in species_class_list:
        run_on_species_class(species_class)
    output_list_of_species_class(species_class_list, output_file)
    domain_plot_path = os.path.join(output_folder, species_type, '.'.join([species_type, "domain.plot.png"]))
    completeness_plot_path = os.path.join(output_folder, species_type,
                                          '.'.join([species_type, "completeness.plot.png"]))
    experimental_plot_path = os.path.join(output_folder, species_type,
                                          '.'.join([species_type, "experimental.plot.png"]))
    go_domain_by_species(species_class_list, domain_plot_path, extra_txt=lab_copyright)
    completeness_by_species(species_class_list, completeness_plot_path, extra_txt=lab_copyright)
    experimental_chart_by_species(species_class_list, experimental_plot_path, extra_txt=lab_copyright)


if __name__ == '__main__':
    run_list_of_species_classes(model, 'model')
    run_list_of_species_classes(jgi, 'jgi')
    run_list_of_species_classes(uniprot, 'uniprot')

