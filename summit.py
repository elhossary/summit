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
    parser.add_argument("--condition_names", default="", type=str, nargs="+",
                        help="Names of conditions (space separated), Must be part of the file name")
    parser.add_argument("--min_len", default=30, type=int,
                        help="Minimum allowed annotation length")
    parser.add_argument("--max_len", default=300, type=int,
                        help="Maximum allowed annotation length")
    parser.add_argument("--step_size", default=3, type=int,
                        help="Derivative range")
    parser.add_argument("--threads", default=1, type=int,
                        help="Parallel file processors")
    parser.add_argument("--annotation_type", required=True, type=str,
                        help="Specify a name for the annotation type")
    parser.add_argument("--invert_transcript", default=False, action='store_true',
                        help="Use if the library transcription initiation from the 3' end, example: Term-Seq")
    parser.add_argument("--ignore_coverage", default=10, type=int,
                        help="Ignore coverage up to")
    parser.add_argument("--drop_percentage", default=10, type=int,
                        help="Allow as lowest as percentage of the max point as a signal interruption")
    parser.add_argument("--merge", default=None, choices=[None, "all", "replicate"],
                        help="Annotations from different replicates will be merged according to your choice")
    parser.add_argument("--merge_length_violation", default="merge", choices=["merge", "remove"],
                        help="This gives options for reporting the annotations which violates max_len when merging")
    parser.add_argument("--gff_out", required=True, type=str, help="Path to output GFF file")
    args = parser.parse_args()
    logger.info("Getting list of files")
    refseq_paths = []
    for rs_item in args.refseqs_in:
        for sub_rs_item in glob.glob(rs_item):
            refseq_paths.append(sub_rs_item)

    parsed_wig_paths_df = parse_wig_paths(args.wigs_in)
    conditions_names = parsed_wig_paths_df["condition_name"].unique().tolist()
    output = {}
    for cond_name in conditions_names:
        cond_df = parsed_wig_paths_df[parsed_wig_paths_df["condition_name"] == cond_name]
        all_locs = pd.DataFrame()
        wig_pool = mp.Pool(processes=args.threads)
        processes = []
        for wig_path in cond_df["path"]:
            processes.append(wig_pool.apply_async(process_single_wiggle, args=(wig_path, refseq_paths, args)))
        wiggles_processed = [p.get() for p in processes]
        for wig in wiggles_processed:
            all_locs = all_locs.append(wig, ignore_index=True)
        wig_pool.close()
        all_locs.reset_index(inplace=True, drop=True)
        output[cond_name] = all_locs

    export_df = pd.DataFrame()
    if args.merge is None:
        for k in output.keys():
            export_df = export_df.append(output[k])
        AnnotationExporter(export_df, args).export(prefix="", mode="w")
    elif args.merge == "all":
        for k in output.keys():
            export_df = export_df.append(output[k])
        export_df = ConditionsMerger(export_df, args).merge()
        AnnotationExporter(export_df, args).export(prefix="", mode="w")
    elif args.merge == "replicate":
        for k in output.keys():
            export_df = export_df.append(output[k])
            export_df = ConditionsMerger(export_df, args).merge()
            AnnotationExporter(export_df, args).export(prefix=k, mode="w")
    else:
        print("Error")
        exit(1)


def process_single_wiggle(wig_path, refseq_paths, args):
    peak_annotator_obj = PeakAnnotator(wig_path=wig_path, refseq_paths=refseq_paths, args=args)
    peaks_df = peak_annotator_obj.predict()
    return peaks_df


def parse_wig_paths(wigs_paths):
    wig_info = [x.split(":") for x in wigs_paths]
    wig_info_extended = []
    for i, w in enumerate(wig_info):
        wig_paths = []
        for sub_w_item in glob.glob(w[0]):
            wig_paths.append(sub_w_item)
        wig_info[i].append(wig_paths)

    for w in wig_info:
        if len(w) == 1:
            for wi in w[-1]:
                wig_info_extended.append([wi, ""])
        elif len(w) > 1:
            for wi in w[-1]:
                wig_info_extended.append([wi, w[1]])
        else:
            print("Error")
            exit(1)

    return pd.DataFrame(data=wig_info_extended, columns=["path", "condition_name"])


main()
