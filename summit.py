from summit_libs.peak_annotator import PeakAnnotator
import logging as logger
import argparse
import glob
import multiprocessing as mp
import pandas as pd


def main():
    #np.set_printoptions(suppress=True)
    logger.info("Reading user arguments")
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseqs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--condition_names", required=True, help="", type=str, nargs="+")
    parser.add_argument("--min_len", default=30, help="", type=int)
    parser.add_argument("--max_len", default=300, help="", type=int)
    parser.add_argument("--step_size", default=3, help="", type=int)
    parser.add_argument("--threads", default=2, help="", type=int)
    parser.add_argument("--ignore_coverage", default=10, help="", type=int)
    parser.add_argument("--annotation_type", required=True, help="", type=str)
    parser.add_argument("--separate_conditions", default=False, help="", action='store_true')
    parser.add_argument("--report", default="best", help="", choices=["all", "best", "longest", "shortest"])
    parser.add_argument("--gff_out", required=True, help="", type=str)
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
    count = 0
    for wig_path in wig_paths:
        count += 1
        processes.append(wig_pool.apply_async(process_single_wiggle, args=(wig_path, refseq_paths, args)))
        if count == 2:
            break

    wiggles_processed = [p.get() for p in processes]
    wig_pool.close()
    all_locs = pd.DataFrame()

    for wig in wiggles_processed:
        all_locs = all_locs.append(wig, ignore_index=True)
    all_locs.reset_index(inplace=True, drop=True)

    if not args.separate_conditions:
        all_locs = PeakAnnotator.merge_overlaps(all_locs)

    PeakAnnotator.export_to_gff(all_locs, args)


def process_single_wiggle(wig_path, refseq_paths, args):
    peak_annotator_obj = PeakAnnotator(wig_path=wig_path, refseq_paths=refseq_paths, args=args)
    peaks_df = peak_annotator_obj.predict()
    return peaks_df


main()
