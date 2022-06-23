import os

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib_venn import venn3_unweighted, venn2_unweighted

from util.definitions import annot_type_all, colors_all, aspect_order, annot_type_no_exp, colors_no_exp


def add_pie_to_plt_axes_helper(percentage, total, autotexts_idx, text_size, text_color, pie_count):
    if percentage is True and total is not None:
        for i, a in enumerate(autotexts_idx):
            if text_size is not None:
                a.set_size(text_size)
            if text_color is not None:
                a.set_color(text_color)
            a.set_text("{}\n({:.2f}%)".format(pie_count[i], pie_count[i] * 100 / total))


def add_pie_to_plt_axes(axs, idx, pie_count, pie_labels=None, pie_title=None, title_size=20, text_size=None,
                        text_color=None, pie_color=colors_all, percentage=False, radius=1, total=None):
    if pie_labels is None:
        pie_labels = annot_type_all
    if type(idx) is int:
        try:
            p_idx, tx_idx, autotexts_idx = axs[idx].pie(pie_count, radius=radius, labels=pie_labels, colors=pie_color,
                                                        autopct="")
        except ValueError:
            print(pie_count)
            raise SystemError
        if pie_title is not None:
            if title_size is not None:
                axs[idx].set_title(pie_title, fontsize=title_size)
            else:
                axs[idx].set_title(pie_title)
        add_pie_to_plt_axes_helper(percentage, total, autotexts_idx, text_size, text_color, pie_count)
    elif type(idx) is tuple:
        p_idx, tx_idx, autotexts_idx = axs[idx[0]][idx[1]].pie(pie_count, radius=radius, labels=pie_labels,
                                                               colors=pie_color, autopct="")
        if text_size is not None:
            [_.set_fontsize(text_size) for _ in tx_idx]
        if pie_title is not None:
            if title_size is not None:
                axs[idx[0]][idx[1]].set_title(pie_title, fontsize=title_size)
            else:
                axs[idx[0]][idx[1]].set_title(pie_title)
        add_pie_to_plt_axes_helper(percentage, total, autotexts_idx, text_size, text_color, pie_count)


def add_venn_to_plt_axes(axs, idx, venn_sets, venn_labels, three_part=True, font_size=None, title=None, ylabel=None,
                         xlabel=None):
    if type(idx) is int:
        if three_part is True:
            v = venn3_unweighted([set(v) for v in venn_sets], set_labels=venn_labels, ax=axs[idx])
        else:
            v = venn2_unweighted([set(v) for v in venn_sets], set_labels=venn_labels, ax=axs[idx])
        if title is not None:
            axs[idx].set_title(title)
        if ylabel is not None:
            axs[idx].text(-1, 0, ylabel, size=30, verticalalignment='center')
        if xlabel is not None:
            axs[idx].text(0, -1.4, xlabel, size=30)
        if font_size is not None:
            for text in v.set_labels:
                text.set_fontsize(font_size)
            for x in range(len(v.subset_labels)):
                if v.subset_labels[x] is not None:
                    v.subset_labels[x].set_fontsize(font_size)
    elif type(idx) is tuple:
        if three_part is True:
            v = venn3_unweighted([set(v) for v in venn_sets], set_labels=venn_labels, ax=axs[idx[0]][idx[1]])
        else:
            v = venn2_unweighted([set(v) for v in venn_sets], set_labels=venn_labels, ax=axs[idx[0]][idx[1]])
        if title is not None:
            axs[idx[0]][idx[1]].set_title(title)
        if ylabel is not None:
            axs[idx[0]][idx[1]].text(-1, 0, ylabel, size=30, verticalalignment='center')
        if xlabel is not None:
            axs[idx[0]][idx[1]].text(0, -1.4, xlabel, size=30)
        if font_size is not None:
            for text in v.set_labels:
                text.set_fontsize(font_size)
            for x in range(len(v.subset_labels)):
                if v.subset_labels[x] is not None:
                    v.subset_labels[x].set_fontsize(font_size)


def draw_pie_from_input_helper(axs, annot_type_all):
    axs[0][0].legend(labels=annot_type_all, prop={'size': 20}, bbox_to_anchor=(-0.2, -0.5), loc="lower left")


def draw_pie_from_species(species_class, plot_path=None, title="", plot_venn=True, extra_txt=None):
    if len(getattr(species_class, 'list_of_exp')) == 0:
        with_exp = False
    else:
        with_exp = True
    if title != "" and plot_path is not None and os.path.isdir(os.path.dirname(plot_path)):
        fig = plt.figure(figsize=(20, 12))

        gs = GridSpec(2, 6, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0:2])
        ax2 = fig.add_subplot(gs[0, 2:4])
        ax3 = fig.add_subplot(gs[0, 4:6])
        ax4 = fig.add_subplot(gs[-1, 0:3])
        ax5 = fig.add_subplot(gs[-1, 3:6])
        axs = [[ax1, ax2, ax3], [ax4, ax5]]
    else:
        fig, axs = None, None
    for idx, aspect_name in enumerate(aspect_order):
        if aspect_name == 'Molecular function':
            annot_count = [len(getattr(species_class, 'list_of_mf_exp')),
                           len(getattr(species_class, 'list_of_mf_comp')),
                           len(getattr(species_class, 'list_of_mf_unknown'))]
        elif aspect_name == 'Biological process':
            annot_count = [len(getattr(species_class, 'list_of_bp_exp')),
                           len(getattr(species_class, 'list_of_bp_comp')),
                           len(getattr(species_class, 'list_of_bp_unknown'))]
        elif aspect_name == 'Cellular component':
            annot_count = [len(getattr(species_class, 'list_of_cc_exp')),
                           len(getattr(species_class, 'list_of_cc_comp')),
                           len(getattr(species_class, 'list_of_cc_unknown'))]
        else:
            print('Aspect Error', aspect_name)
            raise SystemError
        if title != "" and plot_path is not None and os.path.isdir(os.path.dirname(plot_path)):
            if with_exp is True:
                add_pie_to_plt_axes(axs, (0, idx), annot_count, pie_title=aspect_name,
                                    pie_labels=[''] * len(annot_count), text_size=15)
            else:
                add_pie_to_plt_axes(axs, (0, idx), annot_count[1:], pie_title=aspect_name, pie_labels=annot_type_no_exp,
                                    pie_color=colors_no_exp, text_size=15)
    if title != "" and plot_path is not None and os.path.isdir(os.path.dirname(plot_path)):
        axs[0][0].set_ylabel("A", rotation=0, fontsize=30, labelpad=60)
    if with_exp is True:
        species_gene_count = [len(getattr(species_class, 'list_of_exp')), len(getattr(species_class, 'list_of_comp')),
                              len(getattr(species_class, 'list_of_unknown'))]
    else:
        species_gene_count = \
            [len(getattr(species_class, 'list_of_comp')), len(getattr(species_class, 'list_of_unknown'))]

    if title != "" and plot_path is not None and os.path.isdir(os.path.dirname(plot_path)):
        if with_exp is True:
            add_pie_to_plt_axes(axs, (1, 0), species_gene_count, pie_title="$\it{" + title + "}$", percentage=True,
                                total=len(getattr(species_class, 'list_of_seqs')),
                                pie_labels=[''] * len(species_gene_count), text_color='white', text_size=15)
        else:
            add_pie_to_plt_axes(axs, (1, 0), species_gene_count, pie_title="$\it{" + title + "}$", percentage=True,
                                total=len(getattr(species_class, 'list_of_seqs')), text_color='white',
                                pie_labels=annot_type_no_exp, pie_color=colors_no_exp, text_size=15)
        axs[1][0].set_ylabel("B", rotation=0, fontsize=30, labelpad=60)
    if plot_venn is True and with_exp is True:
        exp_aspect_sets = []
        exp_aspect_labels = []
        for aspect_name in aspect_order:
            if aspect_name == 'Molecular function':
                exp_aspect_sets.append(getattr(species_class, 'list_of_mf_exp'))
            elif aspect_name == 'Biological process':
                exp_aspect_sets.append(getattr(species_class, 'list_of_bp_exp'))
            elif aspect_name == 'Cellular component':
                exp_aspect_sets.append(getattr(species_class, 'list_of_cc_exp'))
            else:
                print('Aspect Error', aspect_name)
                raise SystemError
            exp_aspect_labels.append(aspect_name)
        if title != "" and plot_path is not None and os.path.isdir(os.path.dirname(plot_path)):
            add_venn_to_plt_axes(axs, (1, 1), exp_aspect_sets, exp_aspect_labels, ylabel='C', font_size=15)
            draw_pie_from_input_helper(axs, annot_type_all)
            if extra_txt is not None:
                axs[0][2].text(0, -2, extra_txt, size=20)
            fig.subplots_adjust(bottom=0.1)
            fig.tight_layout()
            plt.savefig(plot_path)
            plt.close()
    else:
        draw_pie_from_input_helper(axs, annot_type_all)
        if extra_txt is not None:
            axs[0][2].text(0, -2, extra_txt, size=20)
        fig.delaxes(axs[1][1])
        fig.subplots_adjust(bottom=0.1)
        fig.tight_layout()
        plt.savefig(plot_path)
        plt.close()


def go_domain_by_species(list_of_species, plot_path, num_of_rows=3, extra_txt=None):
    annot_labels = annot_type_all
    num_of_species = len(list_of_species)
    plt.figure(1)
    fig, axs = plt.subplots(nrows=num_of_rows, ncols=num_of_species, figsize=(num_of_species * 5, 15))
    for idx, species_class in enumerate(list_of_species):
        species_name = getattr(species_class, 'species_abbr')
        for r, aspect_name in enumerate(aspect_order):
            if aspect_name == 'Molecular function':
                annot_count = [len(getattr(species_class, 'list_of_mf_exp')),
                               len(getattr(species_class, 'list_of_mf_comp')),
                               len(getattr(species_class, 'list_of_mf_unknown'))]
            elif aspect_name == 'Biological process':
                annot_count = [len(getattr(species_class, 'list_of_bp_exp')),
                               len(getattr(species_class, 'list_of_bp_comp')),
                               len(getattr(species_class, 'list_of_bp_unknown'))]
            elif aspect_name == 'Cellular component':
                annot_count = [len(getattr(species_class, 'list_of_cc_exp')),
                               len(getattr(species_class, 'list_of_cc_comp')),
                               len(getattr(species_class, 'list_of_cc_unknown'))]
            axs[0][idx].set_title("$\it{" + species_name + "}$", fontsize=40)
            if len(getattr(species_class, 'list_of_exp')) == 0:
                pie_color = colors_no_exp
                annot_count = annot_count[1:]
            else:
                pie_color = colors_all
            add_pie_to_plt_axes(axs, (r, idx), annot_count, pie_labels=[''] * len(annot_count), pie_color=pie_color)
    for ax, row in zip(axs[:, 0], aspect_order):
        ax.set_ylabel(row, rotation=0, fontsize=30, labelpad=150)
    axs[-1, 0].legend(labels=annot_labels, prop={'size': 30}, loc='lower left', bbox_to_anchor=(0, -1))
    if extra_txt is not None:
        axs[2][len(list_of_species) - 1].text(-0.5, -2, extra_txt, size=20)
    fig.subplots_adjust(bottom=0.25)
    plt.savefig(plot_path, bbox_inches='tight')


def completeness_by_species(list_of_species, plot_path, extra_txt=None):
    annot_labels = annot_type_all
    num_of_species = len(list_of_species)
    plt.figure(1)
    fig, axs = plt.subplots(nrows=1, ncols=num_of_species, figsize=(num_of_species * 5, 10))
    plt.subplots_adjust(wspace=0.5)
    for idx, species_class in enumerate(list_of_species):
        total_gene = len(getattr(species_class, 'list_of_seqs'))
        species_name = getattr(species_class, 'species_abbr')
        species_gene_count = [len(getattr(species_class, 'list_of_exp')), len(getattr(species_class, 'list_of_comp')),
                              len(getattr(species_class, 'list_of_unknown'))]
        if len(getattr(species_class, 'list_of_exp')) == 0:
            pie_color = colors_no_exp
            annot_labels = annot_type_no_exp
            species_gene_count = species_gene_count[1:]
        else:
            pie_color = colors_all
        add_pie_to_plt_axes(axs, idx, species_gene_count, pie_labels=[''] * len(species_gene_count), text_size=15,
                            radius=1, percentage=True, total=total_gene, text_color='white', pie_color=pie_color)
        axs[idx].set_title("$\it{" + species_name + "}$", fontsize=30, pad=50)

    axs[0].legend(labels=annot_labels, prop={'size': 20}, loc='lower left', bbox_to_anchor=(0, -1))
    if extra_txt is not None:
        axs[len(list_of_species) - 1].text(-0.5, -1.7, extra_txt, size=20)
    fig.subplots_adjust(bottom=0.25)
    plt.savefig(plot_path, bbox_inches='tight')


def experimental_chart_by_species(list_of_species, plot_path, extra_txt=None):
    num_of_species = len(list_of_species)
    plt.figure(1)
    fig, axs = plt.subplots(nrows=1, ncols=num_of_species, figsize=(num_of_species * 8, 20))
    plt.subplots_adjust(wspace=0.5)
    for idx, species_class in enumerate(list_of_species):
        species_name = getattr(species_class, 'species_abbr')
        exp_aspect_sets = []
        exp_aspect_labels = []
        for aspect_name in aspect_order:
            if aspect_name == 'Molecular function':
                exp_aspect_sets.append(getattr(species_class, 'list_of_mf_exp'))
            elif aspect_name == 'Biological process':
                exp_aspect_sets.append(getattr(species_class, 'list_of_bp_exp'))
            elif aspect_name == 'Cellular component':
                exp_aspect_sets.append(getattr(species_class, 'list_of_cc_exp'))
            else:
                print('Aspect Error', aspect_name)
                raise SystemError
            exp_aspect_labels.append(aspect_name)
        add_venn_to_plt_axes(axs, idx, exp_aspect_sets, exp_aspect_labels, font_size=20)
        axs[idx].set_title("$\it{" + species_name + "}$", fontsize=30, pad=50)
    if extra_txt is not None:
        axs[len(list_of_species) - 1].text(0, -1.1, extra_txt, size=20)
    plt.savefig(plot_path, bbox_inches='tight')
