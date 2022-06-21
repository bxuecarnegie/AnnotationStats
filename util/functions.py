import copy
import math
import os
import re

import networkx
from tqdm import tqdm

from definitions import aspect_list, experimental_evidence_codes, annot_type_all, annot_type_no_exp, aspect_order


def flatten_list(input_list):
    return [item for sublist in input_list for item in sublist]


def tab_itr(fp, skip="!", indices=None):
    for line in fp:
        line = line.rstrip('\n')
        if line.startswith(skip):
            continue
        info = [i.strip() for i in line.split('\t')]
        try:
            if indices is None:
                yield info
            else:
                yield [info[idx] for idx in indices]
        except IndexError:
            continue


def read_gff(gff_path, attribute_split=';', feature_type=None, attribute_type=None, attribute_must_include=None):
    if feature_type is None:
        feature_type = ['gene']
    if attribute_type is None:
        attribute_type = ['ID=']
    if attribute_must_include is None:
        attribute_must_include = ''
    attr_regex = re.compile('|'.join(['^' + a for a in attribute_type]))
    res = []
    with open(gff_path, 'r') as gp:
        for feature, attribute in tab_itr(gp, indices=[2, -1]):
            if feature in feature_type and attribute_must_include in attribute:
                info = [re.sub(attr_regex, '', attr) for attr in attribute.split(attribute_split)
                        if re.search(attr_regex, attr)]
                res.append(info)
    return res


def read_gaf(gaf_path, id_idx=1, taxon_list=None, db_obj_id_only=True, protein_to_gene=None, id_list=None):
    gaf_dict = {}
    for aspect in sorted(aspect_list.keys()):
        aspect_dict = {
            "exp": set(),
            "pred": set()
        }
        gaf_dict.setdefault(aspect, aspect_dict)
    with open(gaf_path, 'r') as fp:
        for key_id, db_obj_symbol, go_id, ev_code, aspect, db_obj_name, db_obj_syn, db_obj_type, taxon_ids \
                in tqdm(tab_itr(fp, indices=[id_idx, 2, 4, 6, 8, 9, 10, 11, 12])):
            db_obj_syn = set([syn for syn in db_obj_syn.split('|') if syn != ''])
            db_potential_names = set(db_obj_syn | {key_id, db_obj_symbol, db_obj_name})
            if db_obj_id_only is True:
                db_potential_names = {key_id}
            if protein_to_gene is None:
                # potential_gene_ids = set([g.upper() for g in sorted(db_potential_names) if g != ''])
                potential_gene_ids = [g for g in sorted(db_potential_names) if g != '']
            else:
                # potential_gene_ids = set()
                potential_gene_ids = []
                for db_name in sorted(db_potential_names):
                    try:
                        # potential_gene_ids.update([g.lower() for g in sorted(protein_to_gene[db_name])])
                        potential_gene_ids += [g for g in protein_to_gene[db_name]]
                    except KeyError:
                        continue
            if taxon_list is not None:
                taxon_ids = set([re.sub("^taxon:", "", taxon) for taxon in taxon_ids.split('|')])
                if len(taxon_ids & set(taxon_list)) == 0:
                    continue
            if len(potential_gene_ids) > 0:
                # gene_id = sorted(potential_gene_ids)[0]
                if id_list is not None:
                    gene_id = sorted(set(id_list) & set(potential_gene_ids))
                    if len(gene_id) > 0:
                        gene_id = sorted(set(id_list) & set(potential_gene_ids))[0]
                    else:
                        continue
                else:
                    gene_id = potential_gene_ids[0]
                if aspect in gaf_dict and ev_code in experimental_evidence_codes:
                    gaf_dict[aspect]['exp'].add(gene_id)
                elif aspect in gaf_dict:
                    gaf_dict[aspect]['pred'].add(gene_id)
    return gaf_dict


def read_tab_file_to_kv(tab_file_path, key_idx=1, val_idx=0, val_sep=None, one_val=False, key_prefix=""):
    if val_sep is None:
        val_sep = [' ']
    key_val_map = {}
    with open(tab_file_path, 'r') as fp:
        for info in tab_itr(fp, skip='yourlist', indices=[key_idx, val_idx]):
            key = key_prefix + info[0]
            val_list = [re.sub(r'\W+$', '', re.sub(r'^\W+', '', i)) for i in re.split('|'.join(val_sep), info[1])
                        if re.sub(r'\W+$', '', re.sub(r'^\W+', '', i)) != '']
            if len(val_list) == 0:
                continue
            if one_val is True:
                val_list = [val_list[0]]
            try:
                key_val_map[key] += val_list
            except KeyError:
                key_val_map.setdefault(key, val_list)
    return key_val_map


def read_names_dmp(names_dmp_path):
    taxon_to_scientific_name = {}
    with open(names_dmp_path, 'r') as fp:
        for line in fp:
            info = re.split(r'\t+\|\t+', line)
            scientific = False
            for i in info:
                if 'scientific name' in i:
                    scientific = True
            if scientific is True:
                taxon = info[0].strip()
                name = info[1].strip()
                taxon_to_scientific_name.setdefault(taxon, name)
    return taxon_to_scientific_name


# Reverse the labeling
def subset_dict_by_keys(key_list, input_dict):
    subset_dict = {}
    not_found_dict = {}
    input_copy = dict(input_dict)
    for key in input_copy.keys():
        val = copy.deepcopy(input_copy[key])
        if key in key_list:
            subset_dict.setdefault(key, val)
        else:
            not_found_dict.setdefault(key, val)
    return subset_dict, not_found_dict


def highest_power_of_10(n):
    return 10 ** (int(math.log10(n)))


def combine_non_list_val_dict_by_keys(list_of_dicts):
    combined_dict = {}
    for d in tqdm(list_of_dicts):
        d_copy = dict(d)
        for key in d_copy.keys():
            val = d_copy[key]
            if type(d_copy[key]) is not int and type(d_copy[key]) is not str:
                print('Val Type Error', type(d_copy[key]), val)
            else:
                combined_dict.setdefault(key, val)
    return combined_dict


def combine_list_val_dict_by_keys(list_of_dicts):
    combined_dict = {}
    for d in tqdm(list_of_dicts):
        d_copy = dict(d)
        for key in d_copy.keys():
            val = d_copy[key]
            try:
                if type(d_copy[key]) is set:
                    combined_dict[key] += list(sorted(val))
                elif type(d_copy[key]) is not list:
                    combined_dict[key].append(val)
                else:
                    combined_dict[key] += list(val)
            except KeyError:
                if type(d_copy[key]) is set:
                    combined_dict.setdefault(key, list(sorted(val)))
                elif type(d_copy[key]) is not list:
                    combined_dict.setdefault(key, [val])
                else:
                    combined_dict.setdefault(key, list(val))
    return combined_dict


def avg_list_of_lists(list_of_lists):
    avg_list = []
    for idx, l in enumerate(list_of_lists):
        try:
            if idx == 0:
                avg_list = [float(num) for num in l]
            else:
                for i in range(0, len(l)):
                    avg_list[i] += float(l[i])
        except ValueError:
            continue
    for idx, n in enumerate(avg_list):
        avg_list[idx] = n / len(list_of_lists)
    return avg_list


def to_graph(l):
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G


def to_edges(l):
    """
        treat `l` as a Graph and returns it's edges
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current


def write_tsv(list_of_species_name, list_of_total_sequence_num, list_of_species_sequence_count, output_path,
              seq_type='Proteins', with_exp=True, list_of_aspect_pie=None):
    with open(output_path, 'w') as op:
        if with_exp is True:
            if list_of_aspect_pie is None:
                op.write('Species\tNum. of ' + seq_type + '\tNum. of ' + seq_type +
                         ' w/ Experimental evidence\tNum. of ' + seq_type +
                         ' w/ Computational evidence\tNum. of ' + seq_type + ' with no GO annotations\n')
            else:
                annot_type = annot_type_all
                pie_str = '\t'.join(['\t'.join(
                    ['Num. of ' + aspect + ' (' + annot + ')\tPct. of ' + aspect + ' (' + annot + ')' for annot in
                     annot_type]) for aspect in aspect_order])
                op.write('Species\tNum. of ' + seq_type + '\tNum. of ' + seq_type +
                         ' w/ Experimental evidence\tNum. of ' + seq_type +
                         ' w/ Computational evidence\tNum. of ' + seq_type + ' with no GO annotations\t' +
                         pie_str + '\n')
        else:
            if list_of_aspect_pie is None:
                op.write('Species\tNum. of ' + seq_type + '\tNum. of ' + seq_type +
                         ' w/ Computational evidence\tNum. of ' + seq_type + ' with no GO annotations\n')
            else:
                annot_type = annot_type_no_exp
                pie_str = '\t'.join(['\t'.join(
                    ['Num. of ' + aspect + ' (' + annot + ')\tPct. of ' + aspect + ' (' + annot + ')' for annot in
                     annot_type]) for aspect in aspect_order])
                op.write('Species\tNum. of ' + seq_type + '\tNum. of ' + seq_type +
                         ' w/ Computational evidence\tNum. of ' + seq_type + ' with no GO annotations\t' +
                         pie_str + '\n')

        for idx, species_name in enumerate(list_of_species_name):
            total_seq_num = list_of_total_sequence_num[idx]
            species_seq_count = list_of_species_sequence_count[idx]

            if with_exp is True:
                if list_of_aspect_pie is None:
                    op.write('\t'.join([species_name, str(total_seq_num), str(species_seq_count[0]),
                                        str(species_seq_count[1]), str(species_seq_count[2])]) + '\n')
                else:
                    aspect_pie = list_of_aspect_pie[idx]
                    output_str = aspect_pie_str_helper(total_seq_num, aspect_pie, with_exp)
                    op.write('\t'.join([species_name, str(total_seq_num), str(species_seq_count[0]),
                                        str(species_seq_count[1]), str(species_seq_count[2])]) + '\t' +
                             '\t'.join(output_str) + '\n')
            else:
                if list_of_aspect_pie is None:
                    try:
                        op.write('\t'.join([species_name, str(total_seq_num), str(species_seq_count[1]),
                                            str(species_seq_count[2])]) + '\n')
                    except IndexError:
                        op.write('\t'.join([species_name, str(total_seq_num), str(species_seq_count[0]),
                                            str(species_seq_count[1])]) + '\n')
                else:
                    aspect_pie = list_of_aspect_pie[idx]
                    output_str = aspect_pie_str_helper(total_seq_num, aspect_pie, with_exp)
                    try:
                        op.write('\t'.join([species_name, str(total_seq_num), str(species_seq_count[1]),
                                            str(species_seq_count[2])]) + '\t' + '\t'.join(output_str) + '\n')
                    except IndexError:
                        op.write('\t'.join([species_name, str(total_seq_num), str(species_seq_count[0]),
                                            str(species_seq_count[1])]) + '\t' + '\t'.join(output_str) + '\n')
            print(species_name, total_seq_num, sum(species_seq_count), species_seq_count)


def aspect_pie_str_helper(total_seq_num, aspect_pie, with_exp=True):
    aspect_output = [[(str(num), str('{:.2%}'.format(num / total_seq_num))) for num in aspect_pie[aspect]]
                     for aspect in aspect_order]
    if with_exp is True:
        output_str = ['\t'.join(('\t'.join(t[0]), '\t'.join(t[1]), '\t'.join(t[2]))) for t in aspect_output]
    else:
        try:
            output_str = ['\t'.join(('\t'.join(t[1]), '\t'.join(t[2]))) for t in aspect_output]
        except IndexError:
            output_str = ['\t'.join(('\t'.join(t[0]), '\t'.join(t[1]))) for t in aspect_output]
    return output_str

