from summit_libs.wiggle import Wiggle
import numpy as np
from functools import reduce
from Bio import SeqIO
import os
import logging as logger
from scipy.signal import find_peaks
from statistics import mean
import sys
import pandas as pd
import pybedtools as pybed
from io import StringIO


class PeakAnnotator:

    def __init__(self, wig_path, refseq_paths, args):
        self.refseq_paths = refseq_paths
        self.args = args
        self.arr_dict = {}
        self.cond_name = ""
        chrom_sizes = self.get_chrom_sizes(refseq_paths)
        wig_obj = Wiggle(wig_path, chrom_sizes, is_len_extended=True)
        self.transform_wiggle(wig_obj)
        self.wig_orient = wig_obj.orientation
        del wig_obj
        logger.basicConfig(filename='example.log', encoding='utf-8', level=logger.DEBUG)
        logger.getLogger().addHandler(logger.StreamHandler(sys.stdout))

    def predict(self):
        out_df = pd.DataFrame()
        for seqid_key in self.arr_dict.keys():
            # Generate location
            tmp_df = self.generate_locs(self.arr_dict[seqid_key],
                                        True if self.wig_orient == "r" else False,
                                        self.cond_name)
            print(f"Possible {tmp_df.shape[0]} peaks for {self.cond_name}")
            # Score the generated positions
            tmp_df = self.generate_position_score(tmp_df)
            # Filter low scored and Select best from candidates
            tmp_df = self.filter_best_candidates(tmp_df)
            # Group overlaps
            tmp_df = self.generate_grouping_column(tmp_df)
            # Select best from overlaps
            tmp_df = self.select_annotations(tmp_df)
            # append
            tmp_df["seqid"] = seqid_key
            out_df = out_df.append(tmp_df, ignore_index=True)
        out_df.reset_index(inplace=True, drop=True)
        round_cols_dict = {"coverage_mean": 1, "background_mean": 1, "bg_cov_diff": 1,
                           "enrichment": 1, "rp_cov": 1, "fp_cov": 1,
                           "max_cov": 1, "min_cov": 1, "min_max_ratio": 1,
                           "len_enrich_ratio": 1, "mean_max_ratio": 1}
        out_df = out_df.round(round_cols_dict)
        return out_df

    def generate_locs(self, coverage_array, is_reversed, cond_name):
        print(f"Generating all possible locations for: {cond_name}")
        #locs_df = pd.DataFrame(columns=['start', 'end', 'strand'])
        location_col = 0
        raw_coverage_col = 1
        rising_col = 2
        falling_col = 3
        ## small helper functions
        cov_range_func = lambda pos: coverage_array[pos[0] - 1:pos[1] - 1, raw_coverage_col].tolist()
        cov_func = lambda pos: coverage_array[pos - 1, raw_coverage_col]
        backward_bg_func = lambda pos: mean(coverage_array[pos - 7:pos - 3, raw_coverage_col].tolist())
        forward_bg_func = lambda pos: mean(coverage_array[pos + 1:pos + 6, raw_coverage_col].tolist())

        ## Find peaks
        rising_peaks, rising_peaks_props = find_peaks(coverage_array[:, rising_col], distance=5)
        falling_peaks, falling_peaks_props = find_peaks(coverage_array[:, falling_col], distance=5)
        rising_peaks_list = rising_peaks.tolist()
        falling_peaks_list = falling_peaks.tolist()
        falling_peaks_set = set(falling_peaks_list)
        strand = "-" if is_reversed else "+"
        possible_locs = []
        for rp_id, rp in enumerate(rising_peaks_list):
            range_start = max(rp - self.args.max_len, 0) if is_reversed else rp + self.args.min_len
            range_end = rp - self.args.min_len if is_reversed else rp + self.args.max_len
            fp_range = set(range(range_start, range_end, 1))
            possible_fp = list(falling_peaks_set.intersection(fp_range))
            if not possible_fp:
                continue
            for fp in possible_fp:
                upper_loc = int(coverage_array[fp, location_col])
                lower_loc = int(coverage_array[rp, location_col])
                if is_reversed:
                    upper_loc, lower_loc = lower_loc, upper_loc
                rp_cov = cov_func(lower_loc)
                fp_cov = cov_func(upper_loc)
                lower_loc = lower_loc - 1 if rp_cov * 0.50 < cov_func(lower_loc - 1) else lower_loc
                upper_loc = upper_loc + 1 if fp_cov * 0.50 < cov_func(upper_loc + 1) else upper_loc
                if rp_cov == 0 or fp_cov == 0:
                    continue
                if min([rp_cov, fp_cov]) / max([rp_cov, fp_cov]) < 0.10:
                    continue
                cover_range = cov_range_func([lower_loc, upper_loc])
                max_cov = max(cover_range)
                min_cov = min(cover_range)
                mean_cov = mean(cover_range)
                if min_cov <= max_cov * 0.10 or mean_cov <= 1:
                    # if the coverage signal between two points was interrupted by a very low coverage
                    continue
                pos_len = upper_loc - lower_loc + 1
                try:
                    f_bg = forward_bg_func(upper_loc)
                except Exception as e:
                    print(f"{e}")
                    print(f"at position {lower_loc}..{upper_loc}")
                    f_bg = 0
                try:
                    r_bg = backward_bg_func(lower_loc)
                except Exception as e:
                    print(f"{e}")
                    print(f"at position {lower_loc}..{upper_loc}")
                    r_bg = 0
                bg_cov_diff = max([f_bg, r_bg]) - min([f_bg, r_bg])
                bg_mean = mean([f_bg, r_bg])
                enrichment = mean_cov - bg_mean
                rising_peak_cov = cov_func(lower_loc) - r_bg
                falling_peak_cov = cov_func(upper_loc) - f_bg
                mean_max_ratio = mean_cov / max_cov * 100
                len_enrich_ratio = enrichment / pos_len
                min_max_ratio = min_cov / max_cov * 100
                if min_cov < bg_mean:
                    continue
                possible_locs.append({'start': lower_loc,
                                      'end': upper_loc,
                                      'strand': strand,
                                      "position_length": pos_len,
                                      "coverage_mean": mean_cov,
                                      "background_mean": bg_mean,
                                      "bg_cov_diff": bg_cov_diff,
                                      "enrichment": enrichment,
                                      "condition_name": cond_name,
                                      "rp_cov": rising_peak_cov,
                                      "fp_cov": falling_peak_cov,
                                      "max_cov": max_cov,
                                      "min_cov": min_cov,
                                      "min_max_ratio": min_max_ratio,
                                      "len_enrich_ratio": len_enrich_ratio,
                                      "mean_max_ratio": mean_max_ratio})
        return pd.DataFrame(data=possible_locs)

    def select_annotations(self, locs_df):
        print(f"Selecting best candidates for: {self.cond_name}")
        groups = locs_df["group"].unique().tolist()
        for group in groups:
            group_df = locs_df[locs_df["group"] == group].copy()
            if group_df.shape[0] == 1:
                continue
            locs_df.drop(group_df.index, inplace=True)
            group_df.sort_values(["score", "position_length"], inplace=True, ascending=[False, False])
            locs_df = locs_df.append(group_df.iloc[0])
        print(f"\t {locs_df.shape[0]} selected peak annotations for: {self.cond_name}")
        return locs_df

    def generate_position_score(self, df):
        print(f"Scoring annotations for: {self.cond_name}")
        df['score'] = 0
        score_col_id = df.columns.get_loc('score')

        for sort_key in ["start", "end"]:
            pos_keys = df[sort_key].unique().tolist()
            for key_pos in pos_keys:
                tmp_df = df[df[sort_key] == key_pos].copy()
                df.drop(tmp_df.index.tolist(), inplace=True)
                tmp_df.sort_values([sort_key, "len_enrich_ratio"], inplace=True, ascending=[True, False])
                tmp_df.iat[0, score_col_id] += 1
                tmp_df.sort_values([sort_key, "enrichment"], inplace=True, ascending=[True, False])
                tmp_df.iat[0, score_col_id] += 1
                tmp_df.sort_values([sort_key, "min_max_ratio"], inplace=True, ascending=[True, True])
                tmp_df.iat[0, score_col_id] += 1
                tmp_df.sort_values([sort_key, "bg_cov_diff"], inplace=True, ascending=[True, True])
                tmp_df.iat[0, score_col_id] += 1
                tmp_df.sort_values([sort_key, "rp_cov"], inplace=True, ascending=[True, False])
                tmp_df.iloc[0, score_col_id] += 1
                tmp_df.sort_values([sort_key, "fp_cov"], inplace=True, ascending=[True, False])
                tmp_df.iloc[0, score_col_id] += 1
                df = df.append(tmp_df)
        return df

    def filter_best_candidates(self, df: pd.DataFrame):
        print(f"Filtering valid candidates for: {self.cond_name}")
        for sort_key in ["start", "end"]:
            key_sorted = df[sort_key].unique().tolist()
            # group_df.sort_values(["start", "enrichment"], inplace=True, ascending=False)
            for s in key_sorted:
                tmp = df[df[sort_key] == s].copy()
                max_score = tmp['score'].max()
                x = tmp[tmp['score'] == max_score].copy()
                df.drop(tmp.index, inplace=True)
                df = df.append(x)
        print(f"\t {df.shape[0]} valid peak annotations for: {self.cond_name}")
        return df

    def _filter_best_candidates(self, df: pd.DataFrame):
        cols = df.columns.tolist()
        out_df = pd.DataFrame(columns=cols)
        print(f"Filtering valid candidates for: {self.cond_name}")
        df_arr = df.to_numpy()
        score_col_id = df.columns.get_loc("score")
        for sort_key in ["start", "end"]:
            key_sorted = df[sort_key].unique().tolist()
            # group_df.sort_values(["start", "enrichment"], inplace=True, ascending=False)
            key_sorted_id = df.columns.get_loc(sort_key)
            for s in key_sorted:
                tmp = df_arr[df_arr[:, key_sorted_id] == s]
                max_score = np.max(tmp[:, score_col_id])
                out_df = out_df.append(pd.DataFrame(data=tmp[tmp[:, score_col_id] == max_score], columns=cols),
                                       ignore_index=True)
        print(f"\t {out_df.shape[0]} valid peak annotations for: {self.cond_name}")
        del df
        return out_df

    def generate_grouping_column(self, df_in):
        print(f"Grouping overlapping annotations for: {self.cond_name}")
        df_in_subset = df_in[["start", "end"]].copy()
        df_in_subset["seqid"] = "chr1"
        df_in_subset = df_in_subset.reindex(columns=["seqid", "start", "end"])
        df_in_str = df_in_subset.to_csv(sep="\t", header=False, index=False)
        bed_locs = str(pybed.BedTool(df_in_str, from_string=True).sort().merge(d=-1))
        merged_bed_locs_df = pd.read_csv(StringIO(bed_locs), names=["seqid", "start", "end"], sep="\t")
        df_in["group"] = None
        group_counter = 1
        for i in merged_bed_locs_df.index:
            mask = df_in["start"].between(merged_bed_locs_df.at[i, "start"], merged_bed_locs_df.at[i, "end"])
            df_in.loc[mask, ["group"]] = group_counter
            group_counter += 1
        return df_in

    def group_df_rows_by_interval(self, df_in):
        print(f"Grouping overlapping annotations for: {self.cond_name}")
        interval_func = lambda row: pd.Interval(left=row["start"], right=row["end"])
        df_in.sort_values(['strand', 'start', 'end'], inplace=True, ascending=[True, True, False])
        df_in.reset_index(inplace=True, drop=True)
        df_in["interval"] = df_in.apply(func=interval_func, axis=1)
        df_in["group"] = np.nan
        df_in.reset_index(inplace=True, drop=True)
        df_len = df_in.shape[0]
        i = 0
        loop_breaker = 0
        while i < df_len:
            loop_breaker += 1
            for j in range(i + 1, df_len, 1):
                if df_in.at[j, "interval"].overlaps(df_in.at[i, "interval"]):
                    df_in.at[i, "group"], df_in.at[j, "group"] = i, i
                else:
                    i = j
                    break
            if loop_breaker > df_len:
                break
        df_in["group"].fillna(df_in.index.to_series(), downcast='infer', inplace=True)
        #df_in.drop(["interval"], inplace=True, axis=1)
        return df_in

    def transform_wiggle(self, wig_obj):
        wig_cols = ["variableStep_chrom", "location", "score"]
        wig_df = wig_obj.get_wiggle()
        self.cond_name = wig_df.iat[0, 1]
        wig_df = wig_df.loc[:, wig_cols]
        wig_df["score"] = wig_df["score"].abs()
        wig_df.loc[wig_df['score'] <= self.args.ignore_coverage, ["score"]] = 0.0
        merged_df = reduce(lambda x, y: pd.merge(x, y, on=["variableStep_chrom", "location"], how='left'),
                           [wig_df.loc[:, wig_cols],
                            wig_obj.to_step_height(self.args.step_size, "start_end").loc[:, wig_cols],
                            wig_obj.to_step_height(self.args.step_size, "end_start").loc[:, wig_cols]])
        merged_df["location"] = merged_df["location"].astype(int)

        for seqid in merged_df["variableStep_chrom"].unique():
            tmp_merged = merged_df[merged_df["variableStep_chrom"] == seqid].drop("variableStep_chrom", axis=1).copy()
            ret_arr = np.absolute(tmp_merged.to_numpy(copy=True))
            self.arr_dict[seqid] = ret_arr

    def get_chrom_sizes(self, fasta_pathes):
        ret_list = []
        for fasta_path in fasta_pathes:
            logger.info(f"Parsing reference sequence: {os.path.basename(fasta_path)}")
            f_parsed = SeqIO.parse(fasta_path, "fasta")
            logger.info(f"Parsed  reference sequence: {os.path.basename(fasta_path)}")
            for seq_record in f_parsed:
                ret_list.append({"seqid": seq_record.id,
                                 "size": len(seq_record.seq),
                                 "fasta": os.path.basename(fasta_path)})
            logger.info(f"Chrom sizes added")
        return ret_list
