import re

from tqdm import tqdm

from util.definitions import aspect_list, experimental_evidence_codes, annot_type_all, annot_type_no_exp, aspect_order


def flatten_list(input_list):
    return [item for sublist in input_list for item in sorted(sublist)]


def delim_parser(fp, delim='\t', indices=None, ignore=None, return_line=False):
    if indices is None:
        indices = [0, 1]
    if ignore is None:
        ignore = ['!', '#']
    for idx, line in enumerate(fp):
        if len(line.lstrip()) == 0:
            if return_line is True:
                yield idx, None, line
            else:
                yield idx, None
        elif True in [line.lstrip().startswith(i) for i in ignore]:
            if return_line is True:
                yield idx, None, line
            else:
                yield idx, None
        else:
            line = line.rstrip('\n')
            info = re.split(delim, line)
            try:
                if return_line is True:
                    yield idx, [info[idx] for idx in indices], line
                else:
                    yield idx, [info[idx] for idx in indices]
            except IndexError:
                if return_line is True:
                    yield idx, None, line
                else:
                    yield idx, None


def dict_set_helper(input_dict, key, val):
    try:
        if type(val) is list or type(val) is set:
            input_dict[key].update(val)
        else:
            input_dict[key].add(val)
    except KeyError:
        if type(val) is list or type(val) is set:
            input_dict.setdefault(key, set(val))
        else:
            input_dict.setdefault(key, {val})
    return input_dict


def read_tab_file_to_kv(tab_file_path, key_idx=1, val_idx=0, val_sep=None, one_val=False, key_prefix="", skip=None):
    if val_sep is None:
        val_sep = [' ']
    if skip is None:
        skip = ['yourlist']
    key_val_map = {}
    with open(tab_file_path, 'r') as fp:
        for line_num, info in tqdm(delim_parser(fp, ignore=skip, indices=[key_idx, val_idx])):
            if info is None:
                continue
            key = key_prefix + info[0]
            val_list = [re.sub(r'\W+$', '', re.sub(r'^\W+', '', i)) for i in re.split('|'.join(val_sep), info[1])
                        if re.sub(r'\W+$', '', re.sub(r'^\W+', '', i)) != '']
            if len(val_list) == 0:
                continue
            if one_val is True:
                val_list = [val_list[0]]
            key_val_map = dict_set_helper(key_val_map, key, val_list)
    return key_val_map


def read_tair2uniprot(file_path, agi_list):
    ids_to_agi = {}
    with open(file_path, 'r') as fp:
        for line_num, info in tqdm(delim_parser(fp, indices=[0, 1, 2])):
            if info is not None:
                # remove isoforms
                uniprot_id = re.sub('-\d$', '', info[0])
                tair_id = info[1]
                agi_id = info[2]
                ids_to_agi = dict_set_helper(ids_to_agi, tair_id, agi_id)
                ids_to_agi = dict_set_helper(ids_to_agi, uniprot_id, agi_id)
        for agi in sorted(set(agi_list)):
            ids_to_agi = dict_set_helper(ids_to_agi, agi, agi)
    return ids_to_agi


# gff bio_type work around
bio_type_dict = {
    "biotype": "protein_coding",
    "locus_type": "protein_coding",
    "so_term_name": "protein_coding_gene"
}


def read_gff(gff_path, attribute_split=';', feature_type=None, attribute_type=None, bio_type=None):
    if feature_type is None:
        feature_type = ['gene']
    if attribute_type is None:
        attribute_type = ['ID=']
    attr_regex = re.compile('|'.join(['^' + a for a in attribute_type]))
    res = []
    with open(gff_path, 'r', errors='replace') as gp:
        for line_num, info in tqdm(delim_parser(gp, indices=[2, -1])):
            if info is None:
                continue
            feature, attribute = info
            if feature in feature_type:
                info = [re.sub(attr_regex, '', attr) for attr in attribute.split(attribute_split)
                        if re.search(attr_regex, attr)]
                if bio_type is not None:
                    if type(bio_type) is str:
                        bio_type = [bio_type]
                    bio_type_regex = re.compile('|'.join([b for b in bio_type]))
                    if re.search(bio_type_regex, attribute):
                        res.append(info)
                else:
                    res.append(info)
    return res


def lists_of_seq_helper(gene_id, aspect, list_of_ev, list_of_mf, list_of_mf_of_ev, list_of_cc, list_of_cc_of_ev,
                        list_of_bp, list_of_bp_of_ev):
    list_of_ev.add(gene_id)
    if aspect == 'F' or aspect == 'molecular_function':
        list_of_mf.add(gene_id)
        list_of_mf_of_ev.add(gene_id)
    elif aspect == 'C' or aspect == 'cellular_component':
        list_of_cc.add(gene_id)
        list_of_cc_of_ev.add(gene_id)
    elif aspect == 'P' or aspect == 'biological_process':
        list_of_bp.add(gene_id)
        list_of_bp_of_ev.add(gene_id)
    return list_of_ev, list_of_mf, list_of_mf_of_ev, list_of_cc, list_of_cc_of_ev, list_of_bp, list_of_bp_of_ev


def read_gaf(gaf_path, id_idx=1, taxon_list=None, db_obj_id_only=False, protein_to_gene=None, id_list=None,
             gaf_db="TRL"):
    exp_gaf = {}
    comp_gaf = {}
    # ev_code in experimental_evidence_codes
    list_of_exp = set()
    list_of_comp = set()

    # aspect == 'F' (molecular function)
    list_of_mf_seqs = set()
    list_of_mf_exp = set()
    list_of_mf_comp = set()

    # aspect == 'C' (cellular component).
    list_of_cc_seqs = set()
    list_of_cc_exp = set()
    list_of_cc_comp = set()

    # aspect == 'P' (biological process)
    list_of_bp_seqs = set()
    list_of_bp_exp = set()
    list_of_bp_comp = set()

    with open(gaf_path, 'r') as fp:
        for line_num, info, line in \
                tqdm(delim_parser(fp, indices=[id_idx, 2, 4, 6, 8, 9, 10, 11, 12], return_line=True)):
            if info is None:
                continue
            key_id, db_obj_symbol, go_id, ev_code, aspect, db_obj_name, db_obj_syn, db_obj_type, taxon_ids = info
            if taxon_list is not None:
                taxon_ids = set([re.sub("^taxon:", "", taxon) for taxon in taxon_ids.split('|')])
                if len(taxon_ids & set(taxon_list)) == 0:
                    continue
            db_obj_syn = set([syn for syn in db_obj_syn.split('|') if syn != ''])
            db_potential_names = [key_id]
            if db_obj_id_only is not True:
                db_potential_names += sorted(set(db_obj_syn | {db_obj_symbol, db_obj_name}) - {key_id})
            if protein_to_gene is None:
                potential_gene_ids = [g for g in db_potential_names if g != '']
            else:
                potential_gene_ids = []
                for db_name in sorted(db_potential_names):
                    try:
                        potential_gene_ids += [g for g in sorted(protein_to_gene[db_name])]
                    except KeyError:
                        continue
            if len(potential_gene_ids) > 0:
                if id_list is not None:
                    gene_id_list = sorted(set(id_list) & set(potential_gene_ids))
                    if len(gene_id_list) > 0:
                        gene_id = gene_id_list[0]
                    else:
                        gene_id = None
                else:
                    gene_id = potential_gene_ids[0]
            else:
                gene_id = None
            # Replace the DB and DB identifier to the one we used.
            if type(gaf_db) is str and gene_id is not None:
                gaf_line_info = line.split('\t')
                gaf_line_info[0] = gaf_db
                gaf_line_info[1] = gene_id
                line = '\t'.join(gaf_line_info)
            if ev_code in experimental_evidence_codes and gene_id is not None:
                list_of_exp, list_of_mf_seqs, list_of_mf_exp, list_of_cc_seqs, list_of_cc_exp, \
                list_of_bp_seqs, list_of_bp_exp = lists_of_seq_helper(
                    gene_id, aspect, list_of_exp, list_of_mf_seqs, list_of_mf_exp, list_of_cc_seqs, list_of_cc_exp,
                    list_of_bp_seqs, list_of_bp_exp)
                dict_set_helper(exp_gaf, gene_id, line)
            elif gene_id is not None and ev_code != "ND":
                # "ND" removed
                list_of_comp, list_of_mf_seqs, list_of_mf_comp, list_of_cc_seqs, list_of_cc_comp, \
                list_of_bp_seqs, list_of_bp_comp = lists_of_seq_helper(
                    gene_id, aspect, list_of_comp, list_of_mf_seqs, list_of_mf_comp,
                    list_of_cc_seqs, list_of_cc_comp, list_of_bp_seqs, list_of_bp_comp)
                dict_set_helper(comp_gaf, gene_id, line)
    exp_gaf_only = {}
    comp_gaf_only = {}
    for gene_id in sorted(exp_gaf.keys()):
        dict_set_helper(exp_gaf_only, gene_id, exp_gaf[gene_id])
        if gene_id in comp_gaf_only:
            dict_set_helper(exp_gaf_only, gene_id, comp_gaf[gene_id])
    for gene_id in sorted(comp_gaf.keys()):
        if gene_id not in exp_gaf_only:
            dict_set_helper(comp_gaf_only, gene_id, comp_gaf[gene_id])
    return \
        sorted(list_of_exp), sorted(list_of_comp - list_of_exp), \
        sorted(list_of_mf_seqs), sorted(list_of_mf_exp), sorted(list_of_mf_comp - list_of_mf_exp), \
        sorted(list_of_cc_seqs), sorted(list_of_cc_exp), sorted(list_of_cc_comp - list_of_cc_exp), \
        sorted(list_of_bp_seqs), sorted(list_of_bp_exp), sorted(list_of_bp_comp - list_of_bp_exp), \
        exp_gaf_only, comp_gaf_only


def read_phytozome_annotation_info(annotation_info_path, go_dag):
    comp_annotation_info = {}

    # ev_code in experimental_evidence_codes
    list_of_comp = set()

    # res.namespace == 'molecular_function'
    list_of_mf_seqs = set()
    list_of_mf_comp = set()

    # res.namespace == 'cellular_component'
    list_of_cc_seqs = set()
    list_of_cc_comp = set()

    # res.namespace == 'biological_process'
    list_of_bp_seqs = set()
    list_of_bp_comp = set()

    unfound_go_terms = set()
    with open(annotation_info_path, 'r') as fp:
        for line_num, info, line in tqdm(delim_parser(fp, indices=[1, 9], return_line=True)):
            if info is None:
                continue
            locus_id, go_field = info
            go_id_list = [go.strip() for go in re.split(r'[, ]', go_field) if go.strip() != '']
            if len(go_id_list) > 0:
                for go_id in go_id_list:
                    res = go_dag.query_term(go_id)
                    if res is not None:
                        lists_of_seq_helper(locus_id, res.namespace, list_of_comp, list_of_mf_seqs, list_of_mf_comp,
                                            list_of_cc_seqs, list_of_cc_comp, list_of_bp_seqs, list_of_bp_comp)
                        dict_set_helper(comp_annotation_info, locus_id, line)
                    else:
                        unfound_go_terms.add(go_id)
    print("Unfound GO", unfound_go_terms)
    return \
        sorted(list_of_comp), sorted(list_of_mf_seqs), sorted(list_of_mf_comp), \
        sorted(list_of_cc_seqs), sorted(list_of_cc_comp), \
        sorted(list_of_bp_seqs), sorted(list_of_bp_comp), {}, comp_annotation_info


def get_num_and_pct(list_of_seqs, list_of_exp, list_of_comp, list_of_unknown):
    num_of_genes = len(list_of_seqs)
    num_of_genes_exp = len(list_of_exp)
    num_of_genes_comp = len(list_of_comp)
    num_of_genes_unknown = len(list_of_unknown)
    if num_of_genes == 0:
        pct_of_genes_exp = str(0)
        pct_of_genes_comp = str(0)
        pct_of_genes_unknown = str(0)
    else:
        pct_of_genes_exp = '{percent:.2%}'.format(percent=float(num_of_genes_exp) / float(num_of_genes))
        pct_of_genes_comp = '{percent:.2%}'.format(percent=float(num_of_genes_comp) / float(num_of_genes))
        pct_of_genes_unknown = '{percent:.2%}'.format(percent=float(num_of_genes_unknown) / float(num_of_genes))
    return \
        str(num_of_genes), str(num_of_genes_exp), pct_of_genes_exp, str(num_of_genes_comp), pct_of_genes_comp, \
        str(num_of_genes_unknown), pct_of_genes_unknown


def output_list_of_species_class(list_of_species_class, output_path):
    with open(output_path, 'w') as op:
        op.write(
            '\t'.join(['ID', 'Taxon', 'Species', 'Type', 'Sequence Source', 'Sequence URL', 'GAF Source', 'GAF URL',
                       'Name', 'Num. of Genes', 'Num. of Genes w/ Experimental evidence',
                       'Num. of Genes w/ Computational evidence', 'Num. of Genes with no GO annotations',
                       'Num. of Molecular function (Experimental evidence)',
                       'Pct. of Molecular function (Experimental evidence)',
                       'Num. of Molecular function (Predicted)', 'Pct. of Molecular function (Predicted)',
                       'Num. of Molecular function (Unknown)', 'Pct. of Molecular function (Unknown)',
                       'Num. of Biological process (Experimental evidence)',
                       'Pct. of Biological process (Experimental evidence)',
                       'Num. of Biological process (Predicted)', 'Pct. of Biological process (Predicted)',
                       'Num. of Biological process (Unknown)', 'Pct. of Biological process (Unknown)',
                       'Num. of Cellular component (Experimental evidence)',
                       'Pct. of Cellular component (Experimental evidence)',
                       'Num. of Cellular component (Predicted)', 'Pct. of Cellular component (Predicted)',
                       'Num. of Cellular component (Unknown)', 'Pct. of Cellular component (Unknown)',
                       'Date', ]) + '\n')
        for species_class in list_of_species_class:
            date = getattr(species_class, 'date')
            species_id = getattr(species_class, 'species_id')
            taxon = str(getattr(species_class, 'species_taxon'))
            species = getattr(species_class, 'species_name')
            name = getattr(species_class, 'species_abbr')
            species_type = getattr(species_class, 'species_type')

            sequence_source = getattr(species_class, 'gene_source')
            # gene_version = getattr(species_class, 'gene_version')
            sequence_url = getattr(species_class, 'gene_link')

            gaf_source = getattr(species_class, 'annot_source')
            # annot_version = getattr(species_class, 'annot_version')
            gaf_url = getattr(species_class, 'annot_link')

            # Sequences
            list_of_seqs = getattr(species_class, 'list_of_seqs')
            list_of_exp = getattr(species_class, 'list_of_exp')
            list_of_comp = getattr(species_class, 'list_of_comp')
            list_of_unknown = getattr(species_class, 'list_of_unknown')

            num_of_genes, num_of_genes_exp, pct_of_genes_exp, num_of_genes_comp, pct_of_genes_comp, \
            num_of_genes_unknown, pct_of_genes_unknown = \
                get_num_and_pct(list_of_seqs, list_of_exp, list_of_comp, list_of_unknown)

            # Domain
            # list_of_mf_seqs = getattr(species_class, 'list_of_mf_seqs')
            list_of_mf_exp = getattr(species_class, 'list_of_mf_exp')
            list_of_mf_comp = getattr(species_class, 'list_of_mf_comp')
            list_of_mf_unknown = getattr(species_class, 'list_of_mf_unknown')

            _, num_of_mfs_exp, pct_of_mfs_exp, num_of_mfs_comp, pct_of_mfs_comp, \
            num_of_mfs_unknown, pct_of_mfs_unknown = \
                get_num_and_pct(list_of_seqs, list_of_mf_exp, list_of_mf_comp, list_of_mf_unknown)

            # list_of_cc_seqs = getattr(species_class, 'list_of_cc_seqs')
            list_of_cc_exp = getattr(species_class, 'list_of_cc_exp')
            list_of_cc_comp = getattr(species_class, 'list_of_cc_comp')
            list_of_cc_unknown = getattr(species_class, 'list_of_cc_unknown')

            _, num_of_ccs_exp, pct_of_ccs_exp, num_of_ccs_comp, pct_of_ccs_comp, \
            num_of_ccs_unknown, pct_of_ccs_unknown = \
                get_num_and_pct(list_of_seqs, list_of_cc_exp, list_of_cc_comp, list_of_cc_unknown)

            # list_of_bp_seqs = getattr(species_class, 'list_of_bp_seqs')
            list_of_bp_exp = getattr(species_class, 'list_of_bp_exp')
            list_of_bp_comp = getattr(species_class, 'list_of_bp_comp')
            list_of_bp_unknown = getattr(species_class, 'list_of_bp_unknown')

            _, num_of_bps_exp, pct_of_bps_exp, num_of_bps_comp, pct_of_bps_comp, \
            num_of_bps_unknown, pct_of_bps_unknown = \
                get_num_and_pct(list_of_seqs, list_of_bp_exp, list_of_bp_comp, list_of_bp_unknown)

            op.write(
                '\t'.join([species_id, taxon, species, species_type, sequence_source, sequence_url, gaf_source, gaf_url,
                           name, num_of_genes, num_of_genes_exp, num_of_genes_comp, num_of_genes_unknown,
                           num_of_mfs_exp, pct_of_mfs_exp, num_of_mfs_comp, pct_of_mfs_comp,
                           num_of_mfs_unknown, pct_of_mfs_unknown, num_of_bps_exp, pct_of_bps_exp,
                           num_of_bps_comp, pct_of_bps_comp, num_of_bps_unknown, pct_of_bps_unknown,
                           num_of_ccs_exp, pct_of_ccs_exp, num_of_ccs_comp, pct_of_ccs_comp,
                           num_of_ccs_unknown, pct_of_ccs_unknown, date]) + '\n')


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


def read_gaf_og(gaf_path, id_idx=1, taxon_list=None, db_obj_id_only=True, protein_to_gene=None, id_list=None):
    gaf_dict = {}
    for aspect in sorted(aspect_list.keys()):
        aspect_dict = {
            "exp": set(),
            "pred": set()
        }
        gaf_dict.setdefault(aspect, aspect_dict)
    with open(gaf_path, 'r') as fp:
        for line_num, info in delim_parser(fp, indices=[id_idx, 2, 4, 6, 8, 9, 10, 11, 12]):
            if info is None:
                continue
            key_id, db_obj_symbol, go_id, ev_code, aspect, db_obj_name, db_obj_syn, db_obj_type, taxon_ids = info
            db_obj_syn = set([syn for syn in db_obj_syn.split('|') if syn != ''])
            db_potential_names = set(db_obj_syn | {key_id, db_obj_symbol, db_obj_name})
            if db_obj_id_only is True:
                db_potential_names = {key_id}
            if protein_to_gene is None:
                potential_gene_ids = [g for g in sorted(db_potential_names) if g != '']
            else:
                potential_gene_ids = []
                for db_name in sorted(db_potential_names):
                    try:
                        potential_gene_ids += [g for g in protein_to_gene[db_name]]
                    except KeyError:
                        continue
            if taxon_list is not None:
                taxon_ids = set([re.sub("^taxon:", "", taxon) for taxon in taxon_ids.split('|')])
                if len(taxon_ids & set(taxon_list)) == 0:
                    continue
            if len(potential_gene_ids) > 0:
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
