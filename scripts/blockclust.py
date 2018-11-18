#!/usr/bin/env python
import argparse
import os
import pysam
import sys
from collections import defaultdict
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score


def add_chr(old_annotations, new_annotations):
    with open(new_annotations, 'w') as fh_new_annotations:
        with open(old_annotations) as fh_old_annotations:
            for line in fh_old_annotations:
                fh_new_annotations.write('chr' + line)
    return


def extract_d_mcl_clusters(outdir):
    d_mcl_clusters = {}
    print("Ignoring clusters with less than 3 blockgroups")
    fh_mcl = open(os.path.join(outdir, "mcl.out"))
    for n, line in enumerate(fh_mcl, 1):
        f = line.rstrip('\n').split('\t')
        if len(f) < 3:
            continue
        for i in f:
            d_mcl_clusters[i] = str(n)
    fh_mcl.close()

    d_blockgoups = defaultdict(list)
    d_bedlines = defaultdict(list)

    cluster_nr = ""
    bg_count = 0
    fh_annot_bbo = open(os.path.join(outdir, "annotated.bbo"))
    for line in fh_annot_bbo:
        if line.startswith('>'):
            f = line.rstrip('\n').split('\t')
            bg_id = f[0].replace('>cluster_', 'blockgroup_', 1)
            bg_class = f[9]
            bg_count += 1
            if str(bg_count) + ':' + bg_class + ':' + bg_id in d_mcl_clusters:
                cluster_nr = d_mcl_clusters[str(bg_count) + ':' + bg_class + ':' + bg_id]
                d_blockgoups[cluster_nr].append(line)
                bed_entry = '\t'.join([f[1], f[2], f[3],
                                       str(bg_count) + ':' + bg_class + ':' + bg_id + ':' + 'cluster_' + cluster_nr,
                                       f[5], f[4]])
                # if args.no_chr:
                #     bed_entry = "chr" + bed_entry
                d_bedlines[cluster_nr].append(bed_entry)
            else:
                cluster_nr = ""
        elif cluster_nr != "":
            d_blockgoups[cluster_nr].append(line)
    fh_annot_bbo.close()

    mcl_dir = os.path.join(outdir, "mcl_clusters")
    if not os.path.exists(mcl_dir):
        os.makedirs(mcl_dir)

    fh_all_clusters_bed = open(os.path.join(mcl_dir, "all_clusters.bed"), 'w')
    for blockgroup in sorted(d_blockgoups.keys()):
        fh_cluster_bbo = open(os.path.join(mcl_dir, 'cluster_' + blockgroup + '.bbo'), 'w')
        for bg in d_blockgoups[blockgroup]:
            fh_cluster_bbo.write(bg + '\n')
        fh_cluster_bbo.close()
        fh_cluster_bed = open(os.path.join(mcl_dir, 'cluster_' + blockgroup + '.bed'), 'w')
        for bed in d_bedlines[blockgroup]:
            fh_cluster_bed.write(bed + '\n')
            fh_all_clusters_bed.write(bed + '\n')
        fh_cluster_bed.close()
    fh_all_clusters_bed.close()
    return


def mcl_clustering(outdir):
    fh_annot = open(os.path.join(outdir, "blockgroup_annotations.txt"))
    all_bga = {}
    annotations_hash = {}
    for c, line in enumerate(fh_annot, 1):
        f = line.rstrip('\n').split('\t')
        all_bga[c] = f[1] + ':' + f[0]
        annotations_hash[c] = f[0]
    fh_annot.close()

    fh_mtx = open(os.path.join(outdir, "matrix"))
    fh_hclust_in_mtx = open(os.path.join(outdir, "hclust_input.mtx"), 'w')
    fh_tab = open(os.path.join(outdir, "discretized.gspan.tab"), 'w')
    for c, line in enumerate(fh_mtx, 1):
        f = line.rstrip('\n').rstrip().split()
        # trim the white space at the end of the line
        new_line = annotations_hash[c] + " " + line.rstrip('\n').rstrip()
        fh_hclust_in_mtx.write(new_line + '\n')
        for n, i in enumerate(f, 1):
            fh_tab.write(str(c) + ":" + all_bga[c] + "\t" + str(n) + ":" + all_bga[n] + "\t" + i + '\n')
    fh_mtx.close()
    fh_hclust_in_mtx.close()
    fh_tab.close()
    os.system(" ".join(["blockclust_plot.r clust",
                        os.path.join(outdir, "hclust_input.mtx"),
                        os.path.join(outdir, "blockgroup_annotations.txt"),
                        os.path.join(outdir, "hclust_tree.pdf")]))
    # MCL clustering
    os.system(" ".join(["mcl",
                        os.path.join(outdir, "discretized.gspan.tab"),
                        "--abc",
                        "-pi", "20",
                        "-I", "20",
                        "-o", os.path.join(outdir, "mcl.out"),
                        "2>", os.path.join(outdir, "mcl.log")]))
    extract_d_mcl_clusters(outdir)
    return


def write_train_test_targets(outdir, model_dir):
    os.system(' '.join(["cat",
                        os.path.join(model_dir, 'discretized.gspan'),
                        os.path.join(outdir, 'discretized.gspan'),
                        '>', os.path.join(outdir, 'input.gspan')]))
    os.system(' '.join(["cat",
                        os.path.join(model_dir, 'blockgroup_annotations.txt'),
                        os.path.join(outdir, 'blockgroup_annotations.txt'),
                        '>', os.path.join(outdir, 'input.blockgroup_annotations.txt')]))

    train_size = len(open(os.path.join(model_dir, 'blockgroup_annotations.txt')).readlines())
    fh_all_bga = open(os.path.join(outdir, 'input.blockgroup_annotations.txt'))
    fh_all_target = open(os.path.join(outdir, 'input.target'), 'w')
    fh_train_target = open(os.path.join(outdir, 'train.target'), 'w')
    fh_train_ids = open(os.path.join(outdir, 'train.ids'), 'w')
    fh_test_target = open(os.path.join(outdir, 'test.target'), 'w')
    fh_test_ids = open(os.path.join(outdir, 'test.ids'), 'w')

    tgt_id_count = 1
    d_train_class_targets = {}
    d_test_class_targets = {}
    d_all_class_targets = {}
    for c, line in enumerate(fh_all_bga, 1):
        bg_class = line.rstrip('\n').split('\t')[1]
        if bg_class not in d_all_class_targets:
            d_all_class_targets[bg_class] = str(tgt_id_count)
            tgt_id_count += 1
        fh_all_target.write(d_all_class_targets[bg_class] + '\n')
        if c <= train_size:
            fh_train_target.write(d_all_class_targets[bg_class] + '\n')
            fh_train_ids.write(str(c-1) + '\n')
            d_train_class_targets[bg_class] = d_all_class_targets[bg_class]
        else:
            fh_test_target.write(d_all_class_targets[bg_class] + '\n')
            fh_test_ids.write(str(c-1) + '\n')
            d_test_class_targets[bg_class] = d_all_class_targets[bg_class]
    fh_all_bga.close()
    fh_all_target.close()
    fh_train_target.close()
    fh_train_ids.close()
    fh_test_target.close()
    fh_test_ids.close()
    return d_test_class_targets, d_train_class_targets


def extract_eden_parameters(config_file):
    lines = open(config_file).readlines()
    radius = 0
    sequence_degree = 0
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith("Radius:"):
            radius = line.split()[1]
        elif line.startswith("Feature_combination:"):
            sequence_degree = len(line.split()[2].split(',')) + len(line.split()[3].split(','))
    distance = 2 * int(radius) + 1
    return str(sequence_degree), str(radius), str(distance)


def write_predictions(outdir, file, d_predictions):
    fh_bbo = open(os.path.join(outdir, "annotated.bbo"))
    fh_pred = open(file, 'w')
    c = 1
    for line in fh_bbo:
        if not line.startswith('>'):
            continue
        f = line.rstrip('\n').split('\t')
        prediction = d_predictions[c]
        # TODO: find a way to fix nochr issue before the analysis
        # if no_chr:
        #     fh_pred.write("chr")
        fh_pred.write('\t'.join([f[1], f[2], f[3], prediction, f[5], f[4]]) + '\n')
        c += 1
    fh_bbo.close()
    fh_pred.close()
    return


def nearest_neighbour_predictions(outdir, d_class_targets, config_file, bit_size):
    sequence_degree, radius, distance = extract_eden_parameters(config_file)
    os.system(' '.join(["EDeN", "-g DIRECTED", "-b", str(bit_size),
                        "-t", os.path.join(outdir, "input.target"),
                        "-A", os.path.join(outdir, "test.ids"),
                        "-B", os.path.join(outdir, "train.ids"),
                        "-a NEAREST_NEIGHBOR",
                        "-i", os.path.join(outdir, "input.gspan"),
                        "-f SEQUENCE", "-r", radius, "-d", distance, "-M", sequence_degree, "-x 3", "-y", outdir,
                        ">>", os.path.join(outdir, "EDeN.log")]))

    d_target_class = {v: k for k, v in d_class_targets.items()}
    d_bga = map_blockgroup_annotations(outdir)
    
    d_nn_predictions = {}
    fh_knn = open(os.path.join(outdir, "knn_target_value"))
    for c, line in enumerate(fh_knn, 1):
        # predict if it is unknown
        if d_bga[c].split(':')[0] == 'unknown':
            d_targets = defaultdict(int)
            f = line.rstrip('\n').rstrip().split()
            del f[0]
            for target in f:
                d_targets[target] += 1
            for target in d_targets.keys():
                if d_targets[target]/float(len(f)) > 0.5:
                    d_nn_predictions[c] = 'predicted_' + d_target_class[target]
        else:
            d_nn_predictions[c] = d_bga[c].split(':')[0]
        if c not in d_nn_predictions:
            d_nn_predictions[c] = 'unknown'
    fh_knn.close()
    write_predictions(outdir, os.path.join(outdir, "nearest_neighbour_predictions.txt"), d_nn_predictions)


def write_class_targets(rna_class, class_target, type, outdir):
    if not os.path.exists(os.path.join(outdir, rna_class)):
        os.makedirs(os.path.join(outdir, rna_class))
    fh_all_targets = open(os.path.join(outdir, type + '.target'))
    fh_class_targets = open(os.path.join(outdir, rna_class, rna_class + '.' + type + '.target'), 'w')
    for line in fh_all_targets:
        if line.rstrip('\n') == class_target:
            fh_class_targets.write("1\n")
        else:
            fh_class_targets.write("-1\n")
    fh_all_targets.close()
    fh_class_targets.close()
    return


def compute_nn_performance(outdir, rna_class, class_target):
    class_dir = os.path.join(outdir, rna_class)
    if not os.path.exists(class_dir):
        os.makedirs(class_dir)
    write_class_targets(rna_class, class_target, "test", outdir)

    l_true_targets = open(os.path.join(outdir, rna_class, rna_class + ".test.target")).read().splitlines()
    l_true_targets = list(map(int, l_true_targets))
    l_pred_targets = []
    l_pred_values = []
    fh_knn = open(os.path.join(outdir, "knn_target_value"))
    for line in fh_knn:
        f = line.rstrip('\n').rstrip().split()
        s = 0
        for i in range(len(f)):
            if f[i] == class_target:
                s += 1
            else:
                s += -1
        l_pred_values.append(s/float(len(f)))
        if s/float(len(f)) >= 0:
            l_pred_targets.append(1)
        else:
            l_pred_targets.append(-1)
    fh_knn.close()

    # use predicted values for auc and targets for remaining measures
    aucroc = roc_auc_score(l_true_targets, l_pred_values)
    accuracy = accuracy_score(l_true_targets, l_pred_targets)
    precision = precision_score(l_true_targets, l_pred_targets)
    recall = recall_score(l_true_targets, l_pred_targets)
    fscore = f1_score(l_true_targets, l_pred_targets)

    print(rna_class + ":")
    print('accuracy:    {:06f}'.format(accuracy))
    print('precision:   {:06f}'.format(precision))
    print('recall:      {:06f}'.format(recall))
    print('f1-score:    {:06f}'.format(fscore))
    print('aucroc:      {:06f}'.format(aucroc))
    print("")
    return


def collect_model_based_predictions(rna_class, outdir, d_temp_predictions):
    fh_pred = open(os.path.join(outdir, rna_class, "prediction"))
    for c, line in enumerate(fh_pred, 1):
        target = line.rstrip('\n').split()[0]
        if target == "1":
            prediction = "predicted_" + rna_class
            d_temp_predictions[c].append(prediction)
    fh_pred.close()
    return d_temp_predictions


def map_blockgroup_annotations(outdir):
    fh_bg_annot = open(os.path.join(outdir, "blockgroup_annotations.txt"))
    d_bga = {}
    for c, line in enumerate(fh_bg_annot, 1):
        f = line.rstrip('\n').split('\t')
        d_bga[c] = f[1] + ":" + f[0]
    fh_bg_annot.close()
    return d_bga


def model_based_predictions(outdir, model_dir, d_rna_class_targets, config_file, bit_size):
    d_temp_predictions = defaultdict(list)
    for rna_class in sorted(d_rna_class_targets.keys()):
        write_class_targets(rna_class, d_rna_class_targets[rna_class], "train", outdir)
        write_class_targets(rna_class, d_rna_class_targets[rna_class], "test", outdir)
        sequence_degree, radius, distance = extract_eden_parameters(config_file)
        os.system(' '.join(["cp", os.path.join(outdir, "discretized.gspan"),
                           os.path.join(outdir, rna_class, "discretized.gspan")]))
        os.system(' '.join(["EDeN", "-g DIRECTED", "-b", str(bit_size),
                            "-i", os.path.join(outdir, rna_class, "discretized.gspan"),
                            "-f SEQUENCE", "-M", sequence_degree, "-r", radius, "-d", distance, "-a TEST",
                            "-m", os.path.join(model_dir, rna_class + ".model"),
                            "-y", os.path.join(outdir, rna_class), ">>", os.path.join(outdir, "EDeN.log")]))
        d_temp_predictions = collect_model_based_predictions(rna_class, outdir, d_temp_predictions)

    d_bga = map_blockgroup_annotations(outdir)
    d_model_based_predictions = {}
    for bg in sorted(d_bga.keys()):
        bg_class = d_bga[bg].split(':')[0]
        if bg_class == "unknown":
            bg_class = ';'.join(d_temp_predictions[bg])
        if bg_class == "":
            bg_class = "unknown"
        d_model_based_predictions[bg] = bg_class

    write_predictions(outdir, os.path.join(outdir, "model_based_predictions.txt"), d_model_based_predictions)


def compute_model_based_performance(outdir, rna_class):
    fh_pred = open(os.path.join(outdir, rna_class, "prediction"))
    l_true_targets = open(os.path.join(outdir, rna_class, rna_class + ".test.target")).read().splitlines()
    l_true_targets = list(map(int, l_true_targets))
    l_pred_targets = []
    for c, line in enumerate(fh_pred):
        f = line.rstrip('\n').split()
        if float(f[1]) < 0:
            l_pred_targets.append(-1)
        else:
            l_pred_targets.append(1)

    accuracy = accuracy_score(l_true_targets, l_pred_targets)
    precision = precision_score(l_true_targets, l_pred_targets)
    recall = recall_score(l_true_targets, l_pred_targets)
    fscore = f1_score(l_true_targets, l_pred_targets)
    aucroc = roc_auc_score(l_true_targets, l_pred_targets)

    fh_perf = open(os.path.join(outdir, 'performance_measures'), 'a')
    print(rna_class + ":")
    print('accuracy:    {:06f}'.format(accuracy))
    print('precision:   {:06f}'.format(precision))
    print('recall:      {:06f}'.format(recall))
    print('f1-score:    {:06f}'.format(fscore))
    print('aucroc:      {:06f}'.format(aucroc))
    print("")
    fh_perf.write(rna_class + ":\n")
    fh_perf.write('accuracy:    {:06f}\n'.format(accuracy))
    fh_perf.write('precision:   {:06f}\n'.format(precision))
    fh_perf.write('recall:      {:06f}\n'.format(recall))
    fh_perf.write('f1-score:    {:06f}\n'.format(fscore))
    fh_perf.write('aucroc:      {:06f}\n'.format(aucroc))
    fh_perf.close()
    return


def filter_known_blockgroups(in_bbo, out_bbo):
    fh_in_bbo = open(in_bbo)
    d_blockgroups = defaultdict(list)
    consider_cluster = False
    current_blockgroup = None
    for line in fh_in_bbo:
        if line.startswith('>'):
            consider_cluster = False
        if line.startswith('>') and "unknown" not in line:
            current_blockgroup = line
            consider_cluster = True
        else:
            if consider_cluster:
                d_blockgroups[current_blockgroup].append(line)
    fh_in_bbo.close()

    fh_out_bbo = open(out_bbo, 'w')
    for blockgroup in sorted(d_blockgroups.keys()):
        fh_out_bbo.write(blockgroup)
        for tag in d_blockgroups[blockgroup]:
            fh_out_bbo.write(tag)
    fh_out_bbo.close()


def auc_roc(outdir):
    d_blockgroup_class = defaultdict()
    fh_annot_bbo = open(os.path.join(outdir, "annotated.bbo"))
    bg_count = 0
    for line in fh_annot_bbo:
        if line.startswith('>'):
            f = line.rstrip('\n').split('\t')
            bg_class = f[9]
            d_blockgroup_class[bg_count] = bg_class
            bg_count += 1
    fh_annot_bbo.close()

    fh_sim = open(os.path.join(outdir, 'matrix'))
    d_class_specific_aucs = defaultdict(list)
    for n_row, line in enumerate(fh_sim):
        f = line.rstrip('\n').rstrip().split()
        row_class = d_blockgroup_class[n_row]
        d_similarities = defaultdict(lambda: defaultdict(str))
        for i, s in enumerate(f):
            d_similarities[s][i] = d_blockgroup_class[i]

        l_labels = []
        l_similarities = []
        for sim in sorted(d_similarities.keys()):
            for n_col in d_similarities[sim].keys():
                col_class = d_similarities[sim][n_col]
                if row_class == col_class:
                    l_labels.append(1)
                else:
                    l_labels.append(0)
                l_similarities.append(float(sim))
        d_class_specific_aucs[row_class].append(roc_auc_score(l_labels, l_similarities))
    fh_sim.close()

    roc_auc_sum_all = 0
    known_blockgroups = 0
    uknown_blockgroups = 0
    fh_perf = open(os.path.join(outdir, 'performance_measures.txt'), 'w')
    fh_perf.write('{:<20s}{:>5s}{:>12s}\n'.format('ncRNA class',
                                                   'count',
                                                   'AUC ROC'))
    print('{:<20s}{:>5s}{:>12s}'.format('ncRNA class',
                                        'count',
                                        'AUC ROC'))
    for rna_class in sorted(d_class_specific_aucs.keys()):
        roc_auc_sum_class = 0
        if rna_class == 'unknown':
            uknown_blockgroups += 1
        else:
            for auc in d_class_specific_aucs[rna_class]:
                roc_auc_sum_class += float(auc)
                roc_auc_sum_all += float(auc)
                known_blockgroups += 1
            class_avg_roc_auc = roc_auc_sum_class/float(len(d_class_specific_aucs[rna_class]))
            fh_perf.write('{:<20s}{:>5d}{:>12.6f}\n'.format(rna_class,
                                                             len(d_class_specific_aucs[rna_class]),
                                                             class_avg_roc_auc))
            print('{:<20s}{:>5d}{:>12.6f}'.format(rna_class,
                                                 len(d_class_specific_aucs[rna_class]),
                                                 class_avg_roc_auc))

    roc_auc_avg_all = roc_auc_sum_all / float(known_blockgroups)
    fh_perf.write('{:<20s}{:>5d}{:>12.6f}\n'.format('Average',
                                                     known_blockgroups,
                                                     roc_auc_avg_all))
    print('{:<20s}{:>5d}{:>12.6f}'.format('Average',
                                          known_blockgroups,
                                          roc_auc_avg_all))
    fh_perf.close()
    return


def init():
    parser = argparse.ArgumentParser(description='Efficient clustering and classification of non-coding RNAs from short'
                                                 ' read RNA-seq profiles',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-m", '--mode', type=str, choices=["PRE", "ANALYSIS", "POST"], default='ANALYSIS', dest='mode',
                        help="Mode of operation"
                             "PRE = Preprocessing mode. convert from reads BAM to tags BED."
                             "ANALYSIS = Clustering and/or Classification mode."
                             "POST = Post processing such as plotting and annotation with known Rfam families etc.")
    parser.add_argument('-a', '--accept', action='store', dest='accept_annotations',
                        help='Annotations of known ncRNAs in BED format')
    parser.add_argument('-r', '--reject', action='store', dest='reject_annotations',
                        help='Annotations of other known transcripts (eg. protein coding) in BED format')
    parser.add_argument('-t', '--test_input', action='store', dest='test_input',
                        help='Output of preprocessing mode as input.')
    parser.add_argument('-o', '--out', action='store', dest='output_dir',
                        help='Output directory path for the whole analysis')
    parser.add_argument('-f', '--config', action='store', dest='config_file',
                        default=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                             'share', 'blockclust_data', 'blockclust.config'),
                        help='blockClust configuration file.')
    parser.add_argument("-c", '--classify', action='store_true', help="Classify the input blockgroups")
    parser.add_argument("-cm", '--clmode', type=str, choices=["NEAREST", "MODEL"], default='MODEL',
                        dest='classification_mode',
                        help="Type of classification"
                             "MODEL = Model based classification"
                             "NEAREST= Nearest neighbour classification")
    parser.add_argument('-md', '--model_dir', action='store', dest='model_dir',
                        default=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                             'share', 'blockclust_data', 'models'),
                        help='Directory containing trained models for classification')
    parser.add_argument('-cs', '--cmsearch_out', action='store', dest='cmsearch_out',
                        help='Output of cmsearch tool')
    parser.add_argument('-cbed', '--clust_bed', action='store', dest='clusters_bed',
                        help='BED file containing clusters from ANALYSIS mode')
    parser.add_argument('-bam', '--bam', action='store', dest='bam',
                        help='Input bam file')
    parser.add_argument('-tbed', '--tags_bed', action='store', dest='tags_bed',
                        help='BED file of tags')
    parser.add_argument('-tab', '--sim_tab', action='store', dest='sim_tab',
                        help='Tabular file of pairwise blockgroup similarities')
    parser.add_argument('-b', '--bit_size', action='store', type=int, default=15, dest='bit_size',
                        help='Bit size')
    parser.add_argument('-rfam', '--rfam_map', action='store', dest='rfam_map',
                        default=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                             'share', 'blockclust_data', 'rfam_map.txt'),
                        help='Mapping of Rfam families')

    parser.add_argument("-chr", '--no_chr', action='store_true',
                        help="Input blockgroups do not contain 'chr' in the begining of chromosome ids "
                             "(for eg. Ensembl database do not use 'chr'). ")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.1.0')

    args = parser.parse_args()

    if args.mode == 'PRE':
        if not args.bam or not args.tags_bed:
            sys.stderr.write("options '-bam', '-tbed' are mandatory in pre-processing mode")
            exit(1)

        fh_bam = pysam.AlignmentFile(args.bam, "rb")
        d_reads = defaultdict(lambda: defaultdict(list))
        for alignment in fh_bam.fetch():  # until_eof=True also fetches the unmapped reads
            readid = alignment.query_name
            ref_name = alignment.reference_name
            # TODO: remove chrM restriction
            if ref_name.startswith('chrM'):
                continue
            readseq = alignment.get_forward_sequence()
            alignmentref_start = str(alignment.reference_start)
            alignmentref_end = str(alignment.reference_end)
            strand = '+'
            if alignment.is_reverse:
                strand = '-'
            pos = '\t'.join([ref_name, alignmentref_start, alignmentref_end, strand])
            d_reads[readseq]['readid'].append(readid)
            d_reads[readseq]['pos'].append(pos)
        fh_bam.close()

        fh_tags_bed = open(args.tags_bed, 'w')
        for c, readseq in enumerate(sorted(d_reads.keys())):
            read_count = len(set(d_reads[readseq]['readid']))
            nr_of_loci = len(set(d_reads[readseq]['pos']))
            normalized_count = str("%.6f" % round(read_count / float(nr_of_loci), 6))
            for pos in set(sorted(d_reads[readseq]['pos'])):
                [chrom, start, end, strand] = pos.split('\t')
                fh_tags_bed.write('\t'.join([chrom, start, end,
                                             '|'.join(['tag_' + str(c), str(read_count), str(nr_of_loci)]),
                                             normalized_count,
                                             strand]) + '\n')

    elif args.mode == 'ANALYSIS':
        if not args.test_input or not args.config_file \
                or not args.accept_annotations or not args.reject_annotations or not args.output_dir:
            sys.stderr.write("Options '-t', '-a', '-r', '-f' and '-o' are mandatory in analysis mode\n")
            exit(1)

        if args.classify and (not args.model_dir or not args.classification_mode):
            sys.stderr.write("Please provide model directory (-md) and classification mode (-cm)\n")
            exit(1)

        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        accept_annotations = args.accept_annotations
        reject_annotations = args.reject_annotations
        if args.no_chr:
            accept_annotations = os.path.join(args.output_dir, "accept_annotations.bed")
            reject_annotations = os.path.join(args.output_dir, "reject_annotations.bed")
            add_chr(args.accept_annotations, accept_annotations)
            add_chr(args.reject_annotations, reject_annotations)

        # blockClust call
        os.system(' '.join(["blockclust",
                            "-a", accept_annotations,
                            "-r", reject_annotations,
                            "-c", args.config_file,
                            "-o", args.output_dir,
                            "-t", args.test_input,
                            "-b", str(args.bit_size)]))

        known_bbo = os.path.join(args.output_dir, "annotated.known.bbo")
        filter_known_blockgroups(os.path.join(args.output_dir, "annotated.bbo"), known_bbo)

        contain_known_ncrnas = False
        if os.stat(known_bbo).st_size != 0:
            contain_known_ncrnas = True
            known_dir = os.path.join(args.output_dir, "known")
            if not os.path.exists(known_dir):
                os.makedirs(known_dir)
            os.system(' '.join(["blockclust",
                                "-a", accept_annotations,
                                "-r", reject_annotations,
                                "-c", args.config_file,
                                "-o", known_dir,
                                "-t", known_bbo,
                                "-b", str(args.bit_size),
                                "-q", "true"]))
            os.system("cp " + args.config_file + " " + args.output_dir)
            print("=====================================")
            auc_roc(known_dir)
            print("=====================================")

        mcl_clustering(args.output_dir)
        if args.classify:
            print("Performing classification")
            d_test_class_targets, d_train_class_targets = write_train_test_targets(args.output_dir, args.model_dir)
            if args.classification_mode == 'NEAREST':
                print("NEAREST NEIGHBOR PREDICTION")
                nearest_neighbour_predictions(args.output_dir, d_train_class_targets, args.config_file, args.bit_size)
                if contain_known_ncrnas:
                    d_test_class_targets, d_train_class_targets = write_train_test_targets(
                        os.path.join(args.output_dir, "known"), args.model_dir)
                    nearest_neighbour_predictions(os.path.join(args.output_dir, "known"), d_train_class_targets,
                                                  args.config_file, args.bit_size)
                    for test_rna_class in sorted(d_test_class_targets.keys()):
                        compute_nn_performance(os.path.join(args.output_dir, "known"),
                                               test_rna_class, d_test_class_targets[test_rna_class])
            elif args.classification_mode == 'MODEL':
                print("MODEL BASED PREDICTION")
                model_based_predictions(args.output_dir, args.model_dir,
                                        d_train_class_targets, args.config_file, args.bit_size)
                if contain_known_ncrnas:
                    if not os.path.exists(os.path.join(args.output_dir, "known")):
                        os.makedirs((os.path.join(args.output_dir, "known")))
                    d_test_class_targets, d_train_class_targets = write_train_test_targets(
                        os.path.join(args.output_dir, "known"), args.model_dir)
                    model_based_predictions(os.path.join(args.output_dir, "known"), args.model_dir,
                                            d_train_class_targets, args.config_file, args.bit_size)
                    for test_rna_class in sorted(d_test_class_targets.keys()):
                        compute_model_based_performance(os.path.join(args.output_dir, "known"), test_rna_class)
        print("================ DONE ===============\n")

    elif args.mode == "POST":
        if not args.cmsearch_out or not args.clusters_bed or not args.sim_tab:
            sys.stderr.write("options '-cs', '-cbed' are mandatory in post-processing mode")
            exit(1)

        d_cluster_sizes = defaultdict(int)
        d_cluster_entries = defaultdict(list)
        fh_clust_bed = open(args.clusters_bed)
        for line in fh_clust_bed:
            f = line.rstrip('\n').split('\t')
            cf = f[3].split(':')
            d_cluster_sizes[cf[3]] += 1
            d_cluster_entries[cf[3]].append(cf[0] + ":" + cf[1] + ":" + cf[2])
        fh_clust_bed.close()

        d_rfam_family = {}
        fh_rfam = open(args.rfam_map)
        for line in fh_rfam:
            f = line.rstrip('\n').split('\t')
            d_rfam_family[f[0]] = f[1]
        fh_rfam.close()

        d_cms_cluster_rfam_family = defaultdict(lambda: defaultdict(int))
        d_uniq_found = defaultdict(lambda: defaultdict(int))
        # process cmsearch output
        fh_cmsearch = open(args.cmsearch_out)
        for line in fh_cmsearch:
            f = line.rstrip('\n').split()
            if float(f[15]) > 0.000001:
                continue
            d = f[17].split(':')
            cluster_id = d[-1]
            # print(cluster_id, f[0], f[3], d_rfam_family[f[3]])
            if f[0] in d_uniq_found[cluster_id]:
                continue
            d_cms_cluster_rfam_family[cluster_id][d_rfam_family[f[3]]] += 1
            d_uniq_found[cluster_id][f[0]] += 1
        fh_cmsearch.close()

        d_cluster_analysis = defaultdict(lambda: defaultdict(int))
        for cluster in sorted(d_cluster_sizes.keys()):
            cluster_size = d_cluster_sizes[cluster]
            for rfam_class in d_cms_cluster_rfam_family[cluster].keys():
                d_cluster_analysis[cluster][rfam_class] += d_cms_cluster_rfam_family[cluster][rfam_class]
            known_instances = len(d_uniq_found[cluster])
            unknown_instances = cluster_size - known_instances
            d_cluster_analysis[cluster]['unknown'] = unknown_instances

        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        fh_clust_distr = open(os.path.join(args.output_dir, "cluster_distribution.txt"), 'w')
        fh_clust_distr.write("Clusters\tRNA_type\n")
        for cluster in sorted(d_cluster_analysis.keys()):
            for rna_class in sorted(d_cluster_analysis[cluster].keys()):
                if d_cluster_analysis[cluster][rna_class] <= 0:
                    continue
                for i in range(d_cluster_analysis[cluster][rna_class]):
                    fh_clust_distr.write(cluster + '\t' + rna_class + '\n')
        fh_clust_distr.close()
        os.system(' '.join(["blockclust_plot.r hist ", os.path.join(args.output_dir, "cluster_distribution.txt"),
                            os.path.join(args.output_dir, "cluster_distribution.pdf")]))

        # select cluster representatives
        fh_sim_tab = open(args.sim_tab)
        d_similarities = {}
        for line in fh_sim_tab:
            [bg1, bg2, sim] = line.rstrip('\n').split('\t')
            d_similarities[bg1 + '\t' + bg2] = sim
        fh_sim_tab.close()

        l_cluster_candidates = []
        for cluster in d_cluster_entries.keys():
            entries = d_cluster_entries[cluster]
            candidate = ""
            candidate_sim = 0
            for i in entries:
                s = 0
                for j in entries:
                    s += float(d_similarities[i + '\t' + j])
                if s > candidate_sim:
                    candidate = i
                    candidate_sim = s
            l_cluster_candidates.append(candidate)

        fh_out_tab = open(os.path.join(args.output_dir, "cluster_candidate_sim.tab"), 'w')
        fh_out_mtx = open(os.path.join(args.output_dir, "cluster_candidate_sim.mtx"), 'w')
        fh_out_bga = open(os.path.join(args.output_dir, "candidate_bg_annotations.txt"), 'w')
        for i in l_cluster_candidates:
            [id1, rna_class, bg1] = i.split(':')
            fh_out_bga.write(bg1 + "\t" + rna_class + "\n")
            fh_out_mtx.write(bg1)
            for j in l_cluster_candidates:
                fh_out_tab.write(i + "\t" + j + "\t" + d_similarities[i + "\t" + j] + "\n")
                fh_out_mtx.write(" " + d_similarities[i + "\t" + j])
            fh_out_mtx.write('\n')
        fh_out_tab.close()
        fh_out_mtx.close()
        fh_out_bga.close()

        os.system(' '.join(["blockclust_plot.r clust ", os.path.join(args.output_dir, "cluster_candidate_sim.mtx"),
                            os.path.join(args.output_dir, "candidate_bg_annotations.txt"),
                            os.path.join(args.output_dir, "hclust_tree_clusters.pdf")]))

        os.system(' '.join(["mcl", os.path.join(args.output_dir, "cluster_candidate_sim.tab"), "--abc -pi 20 -I 20",
                            "-o", os.path.join(args.output_dir, "mcl_cand.out"), "2>",
                            os.path.join(args.output_dir, "mcl_cand.log")]))
    else:
        sys.stderr.write("mode (-m) can be one of the 'TRAIN', 'ANALYSIS', 'PRE' or 'POST' \n")
        exit(1)


if __name__ == "__main__":
    init()
