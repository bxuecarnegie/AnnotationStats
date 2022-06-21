import os
import re

from Bio.SeqIO import parse

from definitions import date
from util.draw import go_domain_by_species, completeness_by_species, experimental_chart_by_species, draw_pie_from_input
from util.functions import read_gaf, read_gff, read_tab_file_to_kv, flatten_list, write_tsv


def get_seq_from_fasta(input_fasta, rm_version=False, id_split=None, id_idx=0,
                       in_description=None, description_split=None):
    if description_split is None:
        description_split = [" "]
    if id_split is None:
        id_split = [" "]
    seq_set = set()
    for record in parse(input_fasta, "fasta"):
        if in_description is not None:
            description = re.split('|'.join(description_split), record.description)
            if len([re.sub(in_description, '', i) for i in description if re.search(in_description, i)]) > 0:
                seq_id = [re.sub(in_description, '', i) for i in description if re.search(in_description, i)][0]
            else:
                continue
        else:
            seq_id = record.id
        if rm_version is True:
            seq_id = re.sub(r'\.\d+$', '', seq_id)
        if len([i.strip() for i in re.split('|'.join(id_split), seq_id) if i.strip() != ""]) > id_idx + 1:
            seq_id = [i.strip() for i in re.split('|'.join(id_split), seq_id) if i.strip() != ""][id_idx]
        seq_set.add(seq_id)
    return seq_set


def find_files_of_goa(input_folder):
    gff_file = None
    gaf_file = None
    fasta_file = None
    for f in os.listdir(input_folder):
        file_path = os.path.join(input_folder, f)
        if f.endswith('gaf') and os.path.isfile(file_path):
            gaf_file = file_path
        elif (f.endswith('gff3') or f.endswith('gff')) and os.path.isfile(file_path):
            gff_file = file_path
        elif (f.endswith('fa') or f.endswith('fasta')) and os.path.isfile(file_path):
            fasta_file = file_path
    return gff_file, gaf_file, fasta_file


def process_fasta_folder(input_folder, species_name, map_file="", key_prefix="locus:", gaf_id_idx=1, plot_path=None,
                         rm_version=True, id_split=None, id_idx=0, in_description=None, description_split=None):
    _, gaf_file, fasta_file = find_files_of_goa(input_folder)
    if None not in (gaf_file, fasta_file):
        seq_set = get_seq_from_fasta(fasta_file, rm_version=rm_version, id_split=id_split, id_idx=id_idx,
                                     in_description=in_description, description_split=description_split)
        total_num_of_seq = len(seq_set)
        if os.path.isfile(map_file):
            id_mapping = read_tab_file_to_kv(map_file, key_idx=0, val_idx=1, key_prefix=key_prefix)
            gaf_dict = read_gaf(gaf_file, id_idx=gaf_id_idx, protein_to_gene=id_mapping, id_list=seq_set)
        else:
            gaf_dict = read_gaf(gaf_file, id_idx=gaf_id_idx, id_list=seq_set)
        total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
            draw_pie_from_input(total_num_of_seq, gaf_dict, plot_path=plot_path, title=species_name, extra_txt=date)
        return total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels


def process_gff_folder(input_folder, species_name, gaf_id_idx=1, feature_type=None, plot_path=None,
                       attribute_type=None, attribute_must_include=None):
    if attribute_type is None:
        attribute_type = ["gene_id="]
    gff_file, gaf_file, _ = find_files_of_goa(input_folder)
    if None not in (gff_file, gaf_file):
        gff_list = read_gff(gff_file, feature_type=feature_type,
                            attribute_type=attribute_type, attribute_must_include=attribute_must_include)
        seq_set = set(flatten_list(gff_list))
        total_num_of_seq = len(seq_set)
        gaf_dict = read_gaf(gaf_file, id_idx=gaf_id_idx, id_list=seq_set)
        total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
            draw_pie_from_input(total_num_of_seq, gaf_dict, plot_path=plot_path, title=species_name, extra_txt=date)
        return total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels


def ara_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count, list_of_exp_aspect_sets, list_of_exp_aspect_labels):
    input_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/A.thaliana'
    species_name = 'A.thaliana'
    mapping_file = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/A.thaliana/Araport11_TAIRlocusaccessionID_AGI_mapping.txt'
    plot_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/plots_goa/A.thaliana.plots.png'
    total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
        process_fasta_folder(input_folder, species_name, map_file=mapping_file, plot_path=plot_path, rm_version=True)
    list_of_species_name.append(species_name)
    list_of_total_seq.append(total_num_of_seq)
    list_of_aspect_pie.append(aspect_pie)
    list_of_species_seq_count.append(species_seq_count)
    list_of_exp_aspect_sets.append(exp_aspect_sets)
    list_of_exp_aspect_labels.append(exp_aspect_labels)


def roundworm_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count, list_of_exp_aspect_sets, list_of_exp_aspect_labels):
    input_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/C.elegans'
    species_name = 'C.elegans'
    plot_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/plots_goa/C.elegans.plots.png'
    total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
        process_fasta_folder(input_folder, species_name, plot_path=plot_path, in_description="gene=")
    list_of_species_name.append(species_name)
    list_of_total_seq.append(total_num_of_seq)
    list_of_aspect_pie.append(aspect_pie)
    list_of_species_seq_count.append(species_seq_count)
    list_of_exp_aspect_sets.append(exp_aspect_sets)
    list_of_exp_aspect_labels.append(exp_aspect_labels)


def zebrafish_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count, list_of_exp_aspect_sets, list_of_exp_aspect_labels):
    input_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/D.rerio'
    species_name = 'D.rerio'
    plot_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/plots_goa/D.rerio.plots.png'
    total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
        process_gff_folder(input_folder, species_name, plot_path=plot_path)
    list_of_species_name.append(species_name)
    list_of_total_seq.append(total_num_of_seq)
    list_of_aspect_pie.append(aspect_pie)
    list_of_species_seq_count.append(species_seq_count)
    list_of_exp_aspect_sets.append(exp_aspect_sets)
    list_of_exp_aspect_labels.append(exp_aspect_labels)


def fly_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count, list_of_exp_aspect_sets, list_of_exp_aspect_labels):
    input_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/D.melanogaster'
    species_name = 'D.melanogaster'
    plot_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/plots_goa/D.melanogaster.plots.png'
    total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
        process_fasta_folder(input_folder, species_name, plot_path=plot_path, in_description="parent=",
                             description_split=["; ", ","])
    list_of_species_name.append(species_name)
    list_of_total_seq.append(total_num_of_seq)
    list_of_aspect_pie.append(aspect_pie)
    list_of_species_seq_count.append(species_seq_count)
    list_of_exp_aspect_sets.append(exp_aspect_sets)
    list_of_exp_aspect_labels.append(exp_aspect_labels)


def human_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count, list_of_exp_aspect_sets, list_of_exp_aspect_labels):
    input_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/H.sapiens'
    species_name = 'H.sapiens'
    plot_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/plots_goa/H.sapiens.plots.png'
    total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
        process_fasta_folder(input_folder, species_name, plot_path=plot_path, id_split=["\|"], id_idx=1)
    list_of_species_name.append(species_name)
    list_of_total_seq.append(total_num_of_seq)
    list_of_aspect_pie.append(aspect_pie)
    list_of_species_seq_count.append(species_seq_count)
    list_of_exp_aspect_sets.append(exp_aspect_sets)
    list_of_exp_aspect_labels.append(exp_aspect_labels)


def mouse_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count, list_of_exp_aspect_sets, list_of_exp_aspect_labels):
    input_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/M.musculus'
    species_name = 'M.musculus'
    plot_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/plots_goa/M.musculus.plots.png'
    total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
        process_gff_folder(input_folder, species_name, plot_path=plot_path)
    list_of_species_name.append(species_name)
    list_of_total_seq.append(total_num_of_seq)
    list_of_aspect_pie.append(aspect_pie)
    list_of_species_seq_count.append(species_seq_count)
    list_of_exp_aspect_sets.append(exp_aspect_sets)
    list_of_exp_aspect_labels.append(exp_aspect_labels)


def yeast_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count, list_of_exp_aspect_sets, list_of_exp_aspect_labels):
    input_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/S.cerevisiae'
    species_name = 'S.cerevisiae'
    plot_path = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/plots_goa/S.cerevisiae.plots.png'
    total_num_of_seq, aspect_pie, species_seq_count, exp_aspect_sets, exp_aspect_labels = \
        process_fasta_folder(input_folder, species_name, plot_path=plot_path, in_description="SGDID:",
                             description_split=[" ", ","])
    list_of_species_name.append(species_name)
    list_of_total_seq.append(total_num_of_seq)
    list_of_aspect_pie.append(aspect_pie)
    list_of_species_seq_count.append(species_seq_count)
    list_of_exp_aspect_sets.append(exp_aspect_sets)
    list_of_exp_aspect_labels.append(exp_aspect_labels)


def plot_all(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count, list_of_exp_aspect_sets, list_of_exp_aspect_labels):
    plot_folder = '/Users/bxue/Documents/Carnegie/GoPieCharts/Selected_Species/plots_goa'
    domain_plot_path = os.path.join(plot_folder, "domain.plot.png")
    completeness_plot_path = os.path.join(plot_folder, "completeness.plot.png")
    experimental_plot_path = os.path.join(plot_folder, "experimental.plot.png")
    tsv_path = os.path.join(plot_folder, "goa_stats.tsv")

    go_domain_by_species(list_of_aspect_pie, list_of_species_name, domain_plot_path, extra_txt=date)
    completeness_by_species(list_of_species_seq_count, list_of_total_seq, list_of_species_name,
                            completeness_plot_path, extra_txt=date)
    experimental_chart_by_species(list_of_exp_aspect_sets, list_of_exp_aspect_labels, list_of_species_name,
                                  experimental_plot_path, extra_txt=date)
    write_tsv(list_of_species_name, list_of_total_seq, list_of_species_seq_count, tsv_path, list_of_aspect_pie=list_of_aspect_pie)


if __name__ == "__main__":
    fasta_protein = ["A.thaliana", "H.sapiens", "S.cerevisiae"]
    fasta_gene = ["C.elegans", "D.melanogaster"]
    #gene_id=,
    gff_gene = ["D.rerio", "M.musculus"]

    list_of_species_name = []
    list_of_total_seq = []
    list_of_aspect_pie = []
    list_of_species_seq_count = []
    list_of_exp_aspect_sets = []
    list_of_exp_aspect_labels = []

    yeast_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count,
               list_of_exp_aspect_sets, list_of_exp_aspect_labels)
    human_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count,
               list_of_exp_aspect_sets, list_of_exp_aspect_labels)
    ara_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count,
             list_of_exp_aspect_sets, list_of_exp_aspect_labels)
    fly_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count,
             list_of_exp_aspect_sets, list_of_exp_aspect_labels)
    mouse_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count,
               list_of_exp_aspect_sets, list_of_exp_aspect_labels)
    zebrafish_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count,
                   list_of_exp_aspect_sets, list_of_exp_aspect_labels)
    roundworm_help(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count,
                   list_of_exp_aspect_sets, list_of_exp_aspect_labels)
    plot_all(list_of_species_name, list_of_total_seq, list_of_aspect_pie, list_of_species_seq_count,
             list_of_exp_aspect_sets, list_of_exp_aspect_labels)



