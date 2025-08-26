import glob
import os
import json
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import random
from matplotlib.patches import FancyArrowPatch, Patch, Polygon
from matplotlib.colors import LinearSegmentedColormap

def parse_AMR_fasta(AMR_fasta):
    with open(AMR_fasta) as i:
        rows = i.read().split(">")[1:]
    amr_genes = set()
    for r in rows:
        amr_genes.add(r.split(";")[0])
    return amr_genes

def parse_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def rv_read(read):
    rv = list(reversed(read))
    return ["+" + g[1:] if g[0] == "-" else "-" + g[1:] for g in rv]

def get_multicopy_AMR_contexts(truth_reads, amr_genes):
    gene_copies = {}
    for contig in truth_reads:
        contig_genes = truth_reads[contig][:]
        for i, g in enumerate(contig_genes):
            if g[1:] in amr_genes:
                if g[1:] not in gene_copies:
                    gene_copies[g[1:]] = []
                prefix = contig_genes[i-15:i]
                suffix = contig_genes[i+1: i+16]
                if len(prefix) < 15:
                    prefix = contig_genes[-15:] + prefix
                if len(suffix) < 15:
                    suffix = suffix + contig_genes[:16]
                path = prefix + [g] + suffix
                if g[0] == "-":
                    path = rv_read(path)
                gene_copies[g[1:]].append(path)
    for g in gene_copies:
        if len(gene_copies[g]) > 1:
            print(g, gene_copies[g], len(gene_copies[g]), "\n")

def assign_colors(shared_genes, amr_genes):
    shared_colors = {}

    # Non-AMR shared genes
    non_amr_genes = sorted(shared_genes - amr_genes)
    num_colors = len(non_amr_genes)

    anchor_colors = [
        "#332288",  # dark blue
        "#117733",  # green
        "#DDCC77",  # sand
        "#88CCEE",  # cyan
        "#CC6677",  # muted red
        "#AA4499",  # purple
        "#44AA99",  # teal
        "#999933",  # olive
        "#882255",  # wine
        "#DDDDDD",  # grey (midpoint)
        "#EE7733",  # orange
        "#000000",  # black (for fallback)
    ]
    cmap = LinearSegmentedColormap.from_list("colorblind_diverging", anchor_colors, N=max(num_colors, 3))

    # Create interpolated colors for the number of required genes
    palette = [cmap(i / (num_colors - 1)) for i in range(num_colors)] if num_colors > 1 else [cmap(0.5)]
    random.shuffle(palette)
    for i, gene in enumerate(non_amr_genes):
        shared_colors[gene] = palette[i]

    # Explicit red for AMR genes
    for gene in shared_genes:
        if gene in amr_genes:
            shared_colors[gene] = "red"  # RGB red

    return shared_colors

import matplotlib.pyplot as plt
from collections import Counter

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from collections import Counter

def plot_contexts(plot_data):
    for sample, genes in plot_data.items():
        all_genes = []
        for context_lists in genes.values():
            for context in context_lists:
                all_genes.extend([g[1:].replace("_", "") for g in context])

        gene_counts = Counter(all_genes)
        shared_genes = {g for g, count in gene_counts.items() if count > 1}

        amr_gene = list(genes.keys())[0].split(" copy")[0].replace("_", "")
        shared_colors = assign_colors(shared_genes, {amr_gene})
        num_genes = len(genes)

        fig, axes = plt.subplots(num_genes, 1, figsize=(15, 3 * num_genes), squeeze=False)
        fig.suptitle(f"{sample}", fontsize=14, x=0.55, y=1.02)

        for ax_idx, (gene, contexts) in enumerate(genes.items()):
            ax = axes[ax_idx][0]
            gene_pos = {}  # Track positions of each shared gene

            for context_idx, context in enumerate(contexts):
                y_offset = -context_idx * 0.5
                for i, g in enumerate(context):
                    g_name = g[1:].replace("_", "")
                    direction = g[0]
                    color = shared_colors.get(g_name, 'white')
                    start_x = i if direction == '+' else i + 1
                    end_x = i + 1 if direction == '+' else i

                    arrow = FancyArrowPatch(
                        (start_x, y_offset), (end_x, y_offset),
                        arrowstyle='simple,head_width=1,head_length=1,tail_width=0.5',
                        mutation_scale=10,
                        edgecolor='black',
                        facecolor=color,
                        linewidth=0.5
                    )
                    ax.add_patch(arrow)

                    # Record center position of gene arrow for shared genes
                    if g_name in shared_genes:
                        if g_name not in gene_pos:
                            gene_pos[g_name] = {}
                        if y_offset not in gene_pos[g_name]:
                            gene_pos[g_name][y_offset] = []
                        gene_pos[g_name][y_offset].append(i + 0.5)

                # Add y-label for the row
                label = "truth" if context_idx == 0 else f"gene DBG unitig {context_idx}"
                ax.text(-0.5, y_offset, label, ha='right', va='center', fontsize=14)
            # Draw connection polygons between matching shared genes
            for g_name, y_offsets in gene_pos.items():
                for y1 in y_offsets:
                    for y2 in y_offsets:
                        if y1 < y2:
                            for x1 in y_offsets[y1]:
                                for x2 in y_offsets[y2]:
                                    height = -0.2
                                    polygon = Polygon([
                                        (x1 - 0.4, y1 - height / 2),
                                        (x1 + 0.4, y1 - height / 2),
                                        (x2 + 0.4, y2 + height / 2),
                                        (x2 - 0.4, y2 + height / 2)
                                    ], closed=True, facecolor="grey", alpha=0.3, edgecolor=None)
                                    ax.add_patch(polygon)

            ax.set_ylim(-len(contexts) * 0.5, 0.5)
            ax.set_xlim(-1, max(len(c) for c in contexts) + 1)
            title_name = gene.split(" ")[0].replace("_", "")
            copy = " ".join(gene.split(" ")[1:])
            title_name = f"$\\mathit{{{title_name}}}$ {copy}"
            ax.set_title(title_name)
            ax.axis('off')

        # Add legend for shared genes
        legend_elements = []
        for gene, color in shared_colors.items():
            escaped_gene = gene.replace('_', '')
            label = f"$\\mathit{{{escaped_gene}}}$"
            legend_elements.append(Patch(facecolor=color, edgecolor='black', label=label))

        legend_elements.append(Patch(facecolor="white", edgecolor='black', label="Non-shared gene"))
        if legend_elements:
            fig.legend(
                handles=legend_elements,
                loc='center left',
                bbox_to_anchor=(1.02, 0.5),
                fontsize=14,
                frameon=False,
                ncol=2
            )

        plt.tight_layout()
        os.makedirs("multicopy_gene_plots", exist_ok=True)
        plt.savefig(f"multicopy_gene_plots/{sample}.png", dpi=900, bbox_inches='tight')
        plt.savefig(f"multicopy_gene_plots/{sample}.pdf", dpi=900, bbox_inches='tight')
        plt.close()


AMR_fasta = "AMR_alleles_unified.fa"
amr_genes = parse_AMR_fasta(AMR_fasta)
truth_files = glob.glob("assessment_results/truth_gff_jsons/*.json")
for t in truth_files:
    truth_reads = parse_json(t)
    sample = os.path.basename(t).replace(".json", "")
    print("\n", sample)
    context = get_multicopy_AMR_contexts(truth_reads, amr_genes)

plot_data = {
    "Escherichia_coli_MINF_7C": {
        "aac3IId copy 1": [
            ['+group_2075', '+group_3184', '+group_3419', '+group_3420', '+group_3421', '+group_3183', '+group_2578', '+group_3529', '+group_3422', '+group_3423', '+group_3424', '+finO', '+group_3527', '+group_3325', '+group_3537', '+aac3IId', '+group_3548', '+group_2358'],
            ["+group_2075", "+group_3184", "+group_3419", "+group_3420", "+group_3421", "+group_3183", "+group_2578", "+group_3529", "+group_3422", "+group_3423", "+group_3424", "+finO", "+group_3527", "+group_2222", "+group_3325", "+group_3537", "+aac3IId", "+group_3548", "+group_2358"]
        ],
        "aac3IId copy 2": [
            ['-yubB', '+impA', '+impB', '-group_2958', '-parM', '-group_2053', '-resD', '-group_2970', '+group_3167', '-cia', '+aac3IId', '+group_2450', '+group_1955', '+group_2383', '+group_2210', '+pilL', '+pilM', '+pilN', '+pilO2', '+group_2775', '+group_395', '+pilR', '+pilS'],
            rv_read(["-pilS", "-pilR", "-group_395", "-group_2775", "-pilO2", "-pilN", "-pilM", "-pilL", "-group_2210", "-group_2383", "-group_1955", "-group_2450", "+group_3328", "-CP042618_00163rep_cluster_2371", "-group_2599", "-group_3180", "+yagA", "-group_3548", "-aac3IId", "+cia", "-group_3167", "+group_2970", "+resD", "+group_2053", "+parM", "+group_2958", "-impB", "-impA", "+yubB"])
        ]
    },
    "Escherichia_coli_MSB1_1A": {
        "blaCTXM211NG_0574771 copy 1": [
            ['-group_602', '-cmoM', '+elyC', '-ycbJ', '-kdsB', '-ycaQ', '-lpxK', '-msbA', '-ycaI', '-rpsA', '-cmk', '-ycaL', '-aroA', '-serC', '-ycaP', '+blaCTXM211NG_0574771', '+ycaO', '+focA', '+pflB', '+group_1655', '-ycaK', '+ycaN', '-ycaM', '-ycaD', '+ycaC', '-dmsC', '-dmsB', '-dmsA', '-serS', '-rarA', '-lolA'],
            rv_read(["+lolA", "+rarA", "+serS", "+dmsA", "+dmsB", "+dmsC", "-ycaC", "+ycaD", "+ycaM", "-ycaN", "+ycaK", "-group_1655", "-pflB", "-focA", "-ycaO", "-blaCTXM211NG_0574771", "+ycaP", "+serC", "+aroA", "+ycaL", "+cmk", "+rpsA", "+ycaI", "+msbA", "+lpxK", "+ycaQ", "+kdsB", "+ycbJ", "-elyC", "+cmoM", "+group_602"])
        ],
        "blaCTXM211NG_0574771 copy 2": [
            ['-tetRB', '+tetBNG_0481711', '-tetC', '+group_3258', '+nqrC', '+group_3547', '+group_3546', '+group_3545', '+group_3544', '+group_3543', '+group_3542', '+group_3541', '+group_3540', '+group_3552', '-group_126', '+blaCTXM211NG_0574771', '+sitA', '+sitB', '+sitC', '+sitD', '+group_2871', '-group_3310', '+iucA', '+iucB', '+iucC', '+iucD', '+iutA', '-group_1092', '-group_387', '-vapC', '+resD'],
            ["-tetRB", "+tetBNG_0481711", "-tetC", "+group_3258", "+nqrC", "+group_3547", "+group_3546", "+group_3545", "+group_3544", "+group_3543", "+group_3542", "+group_3541", "+group_3540", "+group_3552", "-group_126", "+blaCTXM211NG_0574771", "+ant3Ihaac6IIdNG_0472181", "+blaOXA", "+catB3", "-group_1661", "-group_233", "-group_3538", "+group_3251", "-000101__NZ_CP012345_00042IncFIB", "+sitA", "+sitB", "+sitC", "+sitD", "+group_2871", "-group_3310", "+iucA", "+iucB", "+iucC", "+iucD", "+iutA", "-group_1092", "-group_387", "-vapC", "+resD"]
        ]
    },
    "Escherichia_coli_MSB1_3B": {
        "ermBNG_0561911 copy 1": [
            ['+group_3419', '+group_3420', '+group_3421', '+group_3183', '+traS', '+group_3529', '+group_3528', '+group_3422', '+group_3423', '+group_3424', '+finO', '+group_3325', '+group_3537', '+ermBNG_0561911', '+group_2795', '+blaCTXM110NG_0489052', '+iroN', '+group_2358', '+blaTEM239NG_0766451', '-group_622', '-tetA', '+tetRA', '-aph6Id', '-aph3Ib', '-sul2NG_0481161', '+mphANG_0479861', '+group_3255', '+group_3256', '-chrA'],
            ["+group_3422", "+group_3423", "+group_3424", "+finO", "+group_3325", "+group_3537", "+ermBNG_0561911", "+group_2795", "+blaCTXM110NG_0489052", "+iroN", "+group_2358", "+blaTEM239NG_0766451", "-group_622", "+group_3548", "-tetA", "+tetRA", "-aph6Id", "-aph3Ib", "-sul2NG_0481161", "+mphANG_0479861", "+group_3255", "+group_3256", "-chrA"],
        ],
        "ermBNG_0561911 copy 2": [
            ['+group_3419', '+group_3420', '+group_3421', '+group_3183', '+group_3085', '+group_3529', '+group_3528', '+group_3422', '+group_3423', '+group_3424', '+group_1497', '+finO', '+group_3527', '+group_2222', '+group_2698', '+ermBNG_0561911', '-stbA', '+group_3536', '+group_3535', '+ychA'],
            ["+group_3422", "+group_3423", "+group_3424", "+group_1497", "+finO", "+group_3527", "+group_2222", "+group_2698", "+group_3325", "-group_1721", "-group_2190", "+group_1843", "+group_3070", "+ermBNG_0561911", "-stbA", "+group_3536", "+group_3535", "+ychA"],
        ]
    },
    "Escherichia_coli_MSB1_6C": {
        "blaCTXM110NG0489052 copy 1": [
            ['-dacD', '+sbcB', '-tsuA', '-plaP', '-yeeY', '-yeeZ', '+hisG', '+hisD', '+group_987', '+hisB', '+hisH', '+hisA', '+hisF', '+hisIE', '-wzzB', '+blaCTXM110NG_0489052', '-ugd', '-gndA', '-group_2840', '-group_795', '-group_1089', '-group_794', '-group_1463', '-group_292', '-group_1352', '-galF', '-gnu', '-wcaM', '-wcaL', '-wcaK', '-wzxC'],
            ["-dacD", "+sbcB", "-tsuA", "-plaP", "-yeeY", "-yeeZ", "+hisG", "+hisD", "+group_987", "+hisB", "+hisH", "+hisA", "+hisF", "+hisIE", "-wzzB", "+blaCTXM110NG_0489052", "-ugd", "-gndA", "-group_2840", "-group_795", "-group_1089", "-group_794", "-group_1463", "-group_292", "-group_1352", "-galF", "-gnu", "-wcaM", "-wcaL", "-wcaK", "-wzxC"]
        ],
        "blaCTXM110NG0489052 copy 2": [
            ['+idnK', '-ahr', '+fimB.intB', '-group_3447', '-group_3447', '-group_1790', '-group_1505', '-group_2888', '-group_2494', '-group_1926', '-pefC', '-group_2474', '-group_2386', '-group_162', '-group_32', '+blaCTXM110NG_0489052', '-group_2691', '+group_3464', '-group_1239', '-group_160', '-group_419', '-group_1651', '-group_808', '+group_2076', '-group_2310', '-group_3467', '-group_2189', '+group_1662', '+group_3468', '-group_2615', '-nepI'],
            ["+idnK", "-ahr", "+fimB.intB", "-group_3447", "-group_1790", "-group_1505", "-group_2888", "-group_2494", "-group_1926", "-pefC", "-group_2474", "-group_2386", "-group_162", "-group_32", "+blaCTXM110NG_0489052", "-group_2691", "+group_3464", "-group_1239", "-group_160", "-group_419", "-group_1651", "-group_808", "+group_2076", "-group_2310", "-group_3467", "-group_3442", "-group_2189", "+group_1662", "+group_3468", "-group_2615", "-nepI"]
        ],
        "blaCTXM110NG0489052 copy 3": [
            ['-group_3183', '-group_3421', '-group_3420', '-group_3419', '-group_3184', '-psiB', '-yubM', '-group_3530', '-group_3531', '-group_3532', '-group_3533', '-group_3534', '-ychA', '-group_3535', '-group_3536', '+stbA', '+group_1154', '+group_3582', '+group_3581', '+pTR2', '+blaCTXM110NG_0489052', '-group_2222', '-group_3527', '-finO', '-group_1497', '-group_3424', '-group_3423'],
            rv_read(["+group_3423", "+group_3424", "+group_1497", "+finO", "+group_3527", "+group_2222", "+group_3325", "+group_3577", "-blaCTXM110NG_0489052", "-pTR2", "-group_3581", "-group_3582", "-group_1154", "-stbA", "+group_3536", "+group_3535", "+ychA", "+group_3534", "+group_3533", "+group_3532", "+group_3531", "+group_3530", "+yubM", "+psiB", "+psiA", "+yubP", "-yubQ", "+group_1980", "+group_3411", "+group_3412", "+group_3413", "+group_3414", "+group_3415", "+group_3416", "+group_3012", "+group_3417", "+group_3418", "+group_2075", "+group_3184", "+group_3419", "+group_3420", "+group_3421", "+group_3183"])
        ],
    },
    "Escherichia_coli_MSB1_7C": {
        "aac3IId copy 1": [
            ['+aac3IId', '+group_3548', '+group_3538', '+yaiT', '+iprA', '-ampH', '+sbmA', '+yaiW', '-ddlA', '+phoA', '+adrA', '-proC', '+yaiI', '+aroL', '+aroM', '-rdgC'],
            ["+aac3IId", "+group_3548", "+group_3538", "+yaiT", "+iprA", "-ampH", "+sbmA", "+yaiW", "-ddlA", "+phoA", "+adrA", "-proC", "+yaiI", "+aroL", "+aroM", "-rdgC"]
        ],
        "aac3IId copy 2": [
            ['-group_3543', '-group_3544', '-group_3545', '-group_3546', '-group_3547', '-nqrC', '-group_3258', '-intI1', '+dfrA17', '+aadA5', '+sul1NG_0480981', '+chrA', '-group_3256', '-group_3255', '-mphANG_0479861', '+aac3IId', '+group_3548', '+group_3538', '+group_233', '+group_1661', '+group_2943', '-group_3035', '-senB', '-group_155', '-group_1710', '-chaN', '-CP019693_00001Col156', '-cdr', '-group_312', '-group_1957', '-group_1811'],
            ["+group_233", "+group_1661", "+group_2943", "-group_3035", "-senB", "-group_155", "-group_1710", "-chaN", "-CP019693_00001Col156", "-cdr", "-group_312", "-group_1957", "-group_1811", "-group_3167", "+resD", "-group_1825", "+group_759", "+group_1219", "-group_905", "+group_1155", "-tetRA", "+tetA", "-group_1426", "-tniA", "-group_1882", "-merA", "-merC", "+merR", "+group_2358", "+blaTEM239NG_0766451", "+group_1832", "+group_3251", "-000101__NZ_CP012345_00042IncFIB", "-group_3250", "-group_3552", "-group_3540", "-group_3541", "-group_3542", "-group_3543", "-group_3544", "-group_3545", "-group_3546", "-group_3547", "-nqrC", "-group_3258", "-intI1", "+dfrA17", "+aadA5", "+sul1NG_0480981", "+chrA", "-group_3256", "-group_3255", "-mphANG_0479861", "+aac3IId", "+group_3548", "+group_3538"],
            ["+aac3IId", "+group_3548", "+group_3538", "+group_233", "+group_1661", "+group_2943"]
        ],
        "mphANG_0479861 copy 1": [
            ['-aroM', '-aroL', '-yaiI', '+proC', '-adrA', '-phoA', '+ddlA', '-yaiW', '-sbmA', '+ampH', '-iprA', '-yaiT', '-group_3538', '-group_3548', '-aac3IId', '+mphANG_0479861', '+group_3255', '+group_3256', '-chrA', '-sul1NG_0480981', '-aadA5', '-dfrA17', '+intI1', '+group_3258', '+nqrC', '+group_3547', '+group_3546', '+group_3545', '+group_3544', '+group_3543', '+group_3542'],
            rv_read(["-group_3542", "-group_3543", "-group_3544", "-group_3545", "-group_3546", "-group_3547", "-nqrC", "-group_3258", "-intI1", "+dfrA17", "+aadA5", "+sul1NG_0480981", "+chrA", "-group_3256", "-group_3255", "-mphANG_0479861", "+aac3IId", "+group_3548", "+group_3538"])
        ],
        "mphANG_0479861 copy 2": [
            ['-group_233', '-group_3538', '-group_3548', '-aac3IId', '+mphANG_0479861', '+group_3255', '+group_3256', '-chrA', '-sul1NG_0480981', '-aadA5', '-dfrA17', '+intI1', '+group_3258', '+nqrC', '+group_3547', '+group_3546', '+group_3545', '+group_3544', '+group_3543', '+group_3542'],
            []
        ],
        "sul1NG_0480981 copy 1": [
            ['-group_3250', '-group_3552', '-group_3540', '-group_3541', '-group_3542', '-group_3543', '-group_3544', '-group_3545', '-group_3546', '-group_3547', '-nqrC', '-group_3258', '-intI1', '+dfrA17', '+aadA5', '+sul1NG_0480981', '+chrA', '-group_3256', '-group_3255', '-mphANG_0479861', '+aac3IId', '+group_3548', '+group_3538', '+yaiT', '+iprA', '-ampH', '+sbmA', '+yaiW', '-ddlA', '+phoA', '+adrA'],
            ["-000101__NZ_CP012345_00042IncFIB", "-group_3250", "-group_3552", "-group_3540", "-group_3541", "-group_3542", "-group_3543", "-group_3544", "-group_3545", "-group_3546", "-group_3547", "-nqrC", "-group_3258", "-intI1", "+dfrA17", "+aadA5", "+sul1NG_0480981", "+chrA", "-group_3256", "-group_3255", "-mphANG_0479861", "+aac3IId", "+group_3548", "+group_3538"]
        ],
        "sul1NG_0480981 copy 2": [
            ['-group_3250', '-group_3552', '-group_3540', '-group_3541', '-group_3542', '-group_3543', '-group_3544', '-group_3545', '-group_3546', '-group_3547', '-nqrC', '-group_3258', '-intI1', '+dfrA17', '+aadA5', '+sul1NG_0480981', '+chrA', '-group_3256', '-group_3255', '-mphANG_0479861', '+aac3IId', '+group_3548', '+group_3538', '+group_233', '+group_1661', '+group_2943', '-group_3035', '-senB', '-group_155', '-group_1710', '-chaN'],
            []
        ],
        "aadA5 copy 1": [
            ['-000101__NZ_CP012345_00042IncFIB', '-group_3250', '-group_3552', '-group_3540', '-group_3541', '-group_3542', '-group_3543', '-group_3544', '-group_3545', '-group_3546', '-group_3547', '-nqrC', '-group_3258', '-intI1', '+dfrA17', '+aadA5', '+sul1NG_0480981', '+chrA', '-group_3256', '-group_3255', '-mphANG_0479861', '+aac3IId', '+group_3548', '+group_3538', '+yaiT', '+iprA', '-ampH', '+sbmA', '+yaiW', '-ddlA', '+phoA'],
            ["-000101__NZ_CP012345_00042IncFIB", "-group_3250", "-group_3552", "-group_3540", "-group_3541", "-group_3542", "-group_3543", "-group_3544", "-group_3545", "-group_3546", "-group_3547", "-nqrC", "-group_3258", "-intI1", "+dfrA17", "+aadA5", "+sul1NG_0480981", "+chrA", "-group_3256", "-group_3255", "-mphANG_0479861", "+aac3IId", "+group_3548", "+group_3538"]
        ],
        "aadA5 copy 2": [
            ['-000101__NZ_CP012345_00042IncFIB', '-group_3250', '-group_3552', '-group_3540', '-group_3541', '-group_3542', '-group_3543', '-group_3544', '-group_3545', '-group_3546', '-group_3547', '-nqrC', '-group_3258', '-intI1', '+dfrA17', '+aadA5', '+sul1NG_0480981', '+chrA', '-group_3256', '-group_3255', '-mphANG_0479861', '+aac3IId', '+group_3548', '+group_3538', '+group_233', '+group_1661', '+group_2943', '-group_3035', '-senB', '-group_155'],
            []
        ],
        "dfrA17 copy 1": [
            ['+group_3251', '-000101__NZ_CP012345_00042IncFIB', '-group_3250', '-group_3552', '-group_3540', '-group_3541', '-group_3542', '-group_3543', '-group_3544', '-group_3545', '-group_3546', '-group_3547', '-nqrC', '-group_3258', '-intI1', '+dfrA17', '+aadA5', '+sul1NG_0480981', '+chrA', '-group_3256', '-group_3255', '-mphANG_0479861', '+aac3IId', '+group_3548', '+group_3538', '+yaiT', '+iprA', '-ampH', '+sbmA', '+yaiW', '-ddlA'],
            ["+group_3251", "-000101__NZ_CP012345_00042IncFIB", "-group_3250", "-group_3552", "-group_3540", "-group_3541", "-group_3542", "-group_3543", "-group_3544", "-group_3545", "-group_3546", "-group_3547", "-nqrC", "-group_3258", "-intI1", "+dfrA17", "+aadA5", "+sul1NG_0480981", "+chrA", "-group_3256", "-group_3255", "-mphANG_0479861", "+aac3IId", "+group_3548", "+group_3538"]
        ],
        "dfrA17 copy 2": [
            ['+group_3251', '-000101__NZ_CP012345_00042IncFIB', '-group_3250', '-group_3552', '-group_3540', '-group_3541', '-group_3542', '-group_3543', '-group_3544', '-group_3545', '-group_3546', '-group_3547', '-nqrC', '-group_3258', '-intI1', '+dfrA17', '+aadA5', '+sul1NG_0480981', '+chrA', '-group_3256', '-group_3255', '-mphANG_0479861', '+aac3IId', '+group_3548', '+group_3538', '+group_233', '+group_1661', '+group_2943', '-group_3035', '-senB', '-group_155'],
            []
        ]
    },
    "Escherichia_coli_MSB1_9D": {
        "sul1NG_0480981 copy 1 & 2": [
            ['+group_2476', '-fipA', '-group_63', '-traJ', '-traK', '+group_2892', '+stbB', '-mucR', '-group_1671', '-ardR', '-group_2936', '-blaCTXM211NG_0574771', '+group_2072', '+group_1686', '-ssb', '-group_1434', '+group_2105', '-intI1', '+dfrA12', '+aadA8NG_0524081', '+sul1NG_0480981', '-group_2042', '+group_2782', '-qnrB41NG_0505042', '+sul1NG_0480981', '+chrA', '-group_3256', '-group_3255', '-mphANG_0479861', '-tniA', '+group_2161', '+traL', '+traB', '+group_1898', '+group_1057', '+group_1916', '+group_1422', '+group_885', '+group_1149', '+group_2476'],
            rv_read(["-group_2476", "-group_1149", "-group_885", "-group_1422", "-group_1916", "-group_1057", "-group_1898", "-traB", "-traL", "-group_2161", "+tniA", "+mphANG_0479861", "+group_3255", "+group_3256", "-chrA", "-sul1NG_0480981", "+qnrB41NG_0505042", "-group_2782", "+group_2042", "-sul1NG_0480981", "-aadA8NG_0524081", "-dfrA12", "+intI1", "-group_2105", "+group_1434", "+ssb", "-group_1686", "-group_2072", "+blaCTXM211NG_0574771", "+group_2936", "+ardR", "+group_1671", "+mucR", "-stbB", "-group_2892", "+traK", "+traJ", "+group_63", "+fipA", "-group_2476"])
        ],
    },
}

plot_contexts(plot_data)