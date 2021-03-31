from summit_libs.peak_annotator import PeakAnnotator
from summit_libs.annotation_exporter import AnnotationExporter
from summit_libs.conditions_merger import ConditionsMerger

import logging as logger
import argparse
import glob
import multiprocessing as mp
import pandas as pd


def main():
    #np.set_printoptions(suppress=True)
    logger.info("Reading user arguments")
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseqs_in", required=True, type=str, nargs="+",
                        help="Fasta files for the reference sequence (space separated)")
    parser.add_argument("--wigs_in", required=True, type=str, nargs="+",
                        help="Wiggle files (space separated)")
    parser.add_argument("--min_len", default=30, type=int,
                        help="Minimum allowed annotation length")
    parser.add_argument("--max_len", default=300, type=int,
                        help="Maximum allowed annotation length")
    parser.add_argument("--step_size", default=3, type=int,
                        help="Derivative range")
    parser.add_argument("--threads", default=6, type=int,
                        help="Parallel file processors")
    parser.add_argument("--ignore_coverage", default=10, type=int,
                        help="Ignore coverage up to")
    parser.add_argument("--annotation_type", required=True, type=str,
                        help="Specify a name for the annotation type")
    parser.add_argument("--separate_conditions", default=False, action='store_true',
                        help="Generates annotations per each condition separately without merging")
    parser.add_argument("--merge_length_violation", default=False, choices=["allow", "remove", "no_merge"],
                        help="This allows the max_len to be violated if the length after merging is longer")

    parser.add_argument("--gff_out", required=True, type=str, help="")
    args = parser.parse_args()
    logger.info("Getting list of files")
    refseq_paths = []
    for rs_item in args.refseqs_in:
        for sub_rs_item in glob.glob(rs_item):
            refseq_paths.append(sub_rs_item)

    wig_paths = []
    for w_item in args.wigs_in:
        for sub_w_item in glob.glob(w_item):
            wig_paths.append(sub_w_item)

    wig_pool = mp.Pool(processes=args.threads)
    processes = []

    for wig_path in wig_paths:
        processes.append(wig_pool.apply_async(process_single_wiggle, args=(wig_path, refseq_paths, args)))

    wiggles_processed = [p.get() for p in processes]
    wig_pool.close()
    all_locs = pd.DataFrame()

    for wig in wiggles_processed:
        all_locs = all_locs.append(wig, ignore_index=True)
    all_locs.reset_index(inplace=True, drop=True)

    if not args.separate_conditions:
        all_locs, unmerged_locs = ConditionsMerger(all_locs, args).merge()
        AnnotationExporter(unmerged_locs, args).export("unmerged")
    AnnotationExporter(all_locs, args).export()


def process_single_wiggle(wig_path, refseq_paths, args):
    peak_annotator_obj = PeakAnnotator(wig_path=wig_path, refseq_paths=refseq_paths, args=args)
    peaks_df = peak_annotator_obj.predict()
    return peaks_df


main()
