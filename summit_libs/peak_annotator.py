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
            print(f"Possible {tmp_df.shape[0]} peaks for: {self.cond_name}")
            # Score the generated positions
            tmp_df = self.generate_position_score(tmp_df, ["start", "end"], False)
            # Filter low scored and Select best from candidates
            tmp_df = self.filter_best_candidates(tmp_df)
            # Group overlaps
            tmp_df = self.generate_grouping_column(tmp_df)
            # Select best from overlaps
            tmp_df = self.generate_position_score(tmp_df, ["group"], True)
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

    def _find_connected_peaks_indexes(self, coverage_array, cond_name) -> pd.DataFrame:
        print(f"Finding possible connected peaks for: {cond_name}")
        ## Find peaks
        rising_col = 2
        falling_col = 3
        if self.args.invert_transcript:
            rising_col, falling_col = falling_col, rising_col
        rising_peaks = find_peaks(coverage_array[:, rising_col])[0]
        falling_peaks = find_peaks(coverage_array[:, falling_col])[0]
        cov_range_func = lambda start, end: coverage_array[start:end, 1].tolist()
        point_info_func = lambda point, col: coverage_array[point, col]
        possible_locs = []
        drop_percent = self.args.drop_percentage / 100
        for rp_id, rp in enumerate(rising_peaks):
            fp_range = np.arange(rp + self.args.min_len, rp + self.args.max_len, 1)
            possible_fps = np.intersect1d(falling_peaks, fp_range).tolist()
            if not possible_fps:
                continue

            possible_fps_heights = [point_info_func(fp, falling_col) for fp in possible_fps]
            #perc = np.percentile(possible_fps_heights, 80)
            #mean_heights = mean(possible_fps_heights)
            for fp_id, fp in enumerate(possible_fps):
                if possible_fps_heights[fp_id] < mean(possible_fps_heights):
                    continue
                pos_range_cov = cov_range_func(rp, fp)
                min_pos_range_cov = min(pos_range_cov)
                if min_pos_range_cov <= max(pos_range_cov) * drop_percent \
                        or min_pos_range_cov <= self.args.ignore_coverage:
                    continue
                possible_locs.append((rp, fp, point_info_func(rp, rising_col), point_info_func(fp, falling_col)))
        possible_locs_df = pd.DataFrame(data=possible_locs, columns=["rp_index", "fp_index", "rp_height", "fp_height"])
        return possible_locs_df

    def generate_locs(self, coverage_array, is_reversed, cond_name) -> pd.DataFrame:
        print(f"Generating all possible locations for: {cond_name}")
        if self.args.invert_transcript and is_reversed:
            is_reversed = False
            strand = "-"
        elif self.args.invert_transcript and not is_reversed:
            is_reversed = True
            strand = "+"
        else:
            strand = "-" if is_reversed else "+"
            pass
        if is_reversed:
            coverage_array = np.flipud(coverage_array)
        connected_peaks_df = self._find_connected_peaks_indexes(coverage_array, cond_name)
        connected_peaks_df = \
            connected_peaks_df.apply(self._apply_row_operations, args=(coverage_array, strand, cond_name), axis=1)
        connected_peaks_df.drop(["rp_index", "fp_index"], inplace=True, axis=1)
        #print(connected_peaks_df.head(100).to_string())
        connected_peaks_df["start"] = connected_peaks_df["start"].astype(int)
        connected_peaks_df["end"] = connected_peaks_df["end"].astype(int)
        return connected_peaks_df

    def _apply_row_operations(self, row, coverage_array, strand, cond_name):
        rising_col = 2
        falling_col = 3
        if self.args.invert_transcript:
            rising_col, falling_col = falling_col, rising_col
        rp_index = int(row["rp_index"])
        fp_index = int(row["fp_index"])
        ## small helper functions
        cov_range_func = lambda start, end: coverage_array[start: end, 1].tolist()
        get_loc_func = lambda point: int(coverage_array[point, 0])
        point_info_func = lambda point, col: coverage_array[point, col]
        #

        # Adjust the ends recursively
        original_rp_index = rp_index
        original_fp_index = fp_index
        rp_index = self._adjust_end_position(rp_index, coverage_array, 0.25, "up", original_rp_index)
        fp_index = self._adjust_end_position(fp_index, coverage_array, 0.25, "down", original_fp_index)
        ##

        loc_pair = (get_loc_func(rp_index), get_loc_func(fp_index))
        lower_loc = min(loc_pair)
        upper_loc = max(loc_pair)
        cov_range = cov_range_func(rp_index, fp_index)
        max_cov = max(cov_range)
        min_cov = min(cov_range)
        mean_cov = mean(cov_range)
        pos_len = upper_loc - lower_loc + 1
        f_pos_s = fp_index + 2
        f_pos_e = fp_index + 5
        r_pos_s = rp_index - 5
        r_pos_e = rp_index - 2

        if self.args.invert_transcript:
            f_pos_s, r_pos_s = r_pos_s, f_pos_s
            f_pos_e, r_pos_e = r_pos_e, f_pos_e

        try:
            f_bg = mean(cov_range_func(f_pos_s, f_pos_e))
        except Exception as e:
            print(f"{e}, at position {f_pos_s}..{f_pos_e}")
            f_bg = 0
        try:
            r_bg = mean(cov_range_func(r_pos_s, r_pos_e))
        except Exception as e:
            print(f"{e}, at position {r_pos_s}..{r_pos_e}")
            r_bg = 0
        if strand == "-":
            f_bg, r_bg = r_bg, f_bg
        bg_mean = mean([f_bg, r_bg])
        enrichment = mean_cov - bg_mean
        rp_height = point_info_func(original_rp_index, rising_col)
        fp_height = point_info_func(original_fp_index, falling_col)
        #rp_height_log = log10(rp_height)
        #fp_height_log = log10(fp_height)

        row["start"] = lower_loc
        row["end"] = upper_loc
        row["strand"] = strand
        row["position_length"] = pos_len
        row["coverage_mean"] = mean_cov
        row["background_mean"] = bg_mean
        row["bg_cov_diff"] = max([f_bg, r_bg]) - min([f_bg, r_bg])
        row["enrichment"] = enrichment
        row["condition_name"] = cond_name
        #row["rp_cov"] = point_info_func(lower_loc, 1) - r_bg
        #row["fp_cov"] = point_info_func(upper_loc, 1) - f_bg
        row["max_cov"] = max_cov
        row["min_cov"] = min_cov
        row["min_max_ratio"] = min_cov / max_cov * 100
        row["len_enrich_ratio"] = enrichment / pos_len
        row["mean_max_ratio"] = mean_cov / max_cov * 100
        #row["log_diff"] = max((rp_height_log, fp_height_log)) - min((rp_height_log, fp_height_log))
        row["peak_heights_diff"] = max(rp_height, fp_height) - min(rp_height, fp_height)
        #row["perc"] = len([x for x in cov_range if x >= np.percentile(cov_range, 80)])
        row['rp_bg'] = r_bg
        row['fp_bg'] = f_bg
        row['score'] = 0
        """
        if min_cov < bg_mean:
            continue
        """
        return row

    def _adjust_end_position(self, p_index: int, coverage_array, percentage, direction, original_peak_index) -> int:
        test_peak = p_index - 1 if direction == "up" else p_index + 1
        if abs(test_peak - original_peak_index) == 3:  # Exit if recursion depth reached 3
            return p_index
        if coverage_array[test_peak, 1] / coverage_array[original_peak_index, 1] >= percentage:
            return self._adjust_end_position(test_peak, coverage_array, percentage, direction, original_peak_index)
        else:
            return p_index

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

    def generate_position_score(self, df, sort_keys, reset_scores):
        if reset_scores:
            df["score"] = 0
        print(f"Scoring annotations by {' and '.join(sort_keys)} for: {self.cond_name}")
        scoring_criteria = {"len_enrich_ratio": False,
                            "enrichment": False,
                            "min_max_ratio": True,
                            "bg_cov_diff": True,
                            "peak_heights_diff": True,
                            "rp_height": False,
                            "fp_height": False}
        """
        ,
                            "rp_bg": True,
                            "fp_bg": True
        """
        #x = {3986432, 3986431, 3986433, 2976306, 2976305, 2976307}
        for sort_key in sort_keys:
            pos_keys = df[sort_key].unique().tolist()
            for key_pos in pos_keys:
                tmp_df = df[df[sort_key] == key_pos].copy()
                """
                if len(x.intersection(tmp_df["start"].tolist())) > 0:
                    print(tmp_df.to_string())
                """
                df.drop(tmp_df.index.tolist(), inplace=True)
                for k, v in scoring_criteria.items():
                    tmp_df = self._calc_score(df=tmp_df, col=k, ascending=v)
                df = df.append(tmp_df)
        return df

    @staticmethod
    def _calc_score(df, col, ascending) -> pd.DataFrame:
        scoring_data = df[col].unique().tolist()
        scoring_range = list(range(0, len(scoring_data) + 1, 1))
        scoring_matrix = dict(zip(scoring_data, scoring_range if ascending else reversed(scoring_range)))
        for indx in df.index:
            df.at[indx, "score"] += scoring_matrix[df.at[indx, col]]
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

    def transform_wiggle(self, wig_obj):
        wig_cols = ["variableStep_chrom", "location", "score"]
        wig_df = wig_obj.get_wiggle()
        self.cond_name = wig_df.iat[0, 1]
        wig_df = wig_df.loc[:, wig_cols]
        wig_df["score"] = wig_df["score"].abs()
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
