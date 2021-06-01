import pandas as pd
import pybedtools as pybed
from itertools import product
from io import StringIO


class ConditionsMerger:

    def __init__(self, df, args):
        self.df = df
        self.args = args

    def merge(self):
        out_df = pd.DataFrame(columns=self.df.columns.tolist())
        seqid_list = self.df["seqid"].unique().tolist()
        combinations = product(seqid_list, ["+", "-"])
        basic_columns = ["seqid", "start", "end"]
        combine_columns = \
            [x for x in self.df.columns.tolist() if x not in basic_columns and x not in
             ["score", "strand", "position_length"]]
        for comb in combinations:
            comb_df = self.df[(self.df["seqid"] == comb[0]) & (self.df["strand"] == comb[1])][basic_columns]
            comb_df_str = comb_df.to_csv(sep="\t", header=False, index=False)
            comb_df_bed = pybed.BedTool(comb_df_str, from_string=True)
            merged_comb_df_bed = str(comb_df_bed.sort().merge(d=-1))
            merged_comb_df = pd.read_csv(StringIO(merged_comb_df_bed), names=basic_columns, sep="\t")
            merged_comb_df["strand"] = comb[1]
            merged_comb_df["score"] = "."
            merged_comb_df["position_length"] = merged_comb_df["end"] - merged_comb_df["start"] + 1
            merged_comb_df["position_length"] = merged_comb_df["position_length"].astype(str)
            for i in merged_comb_df.index:
                tmp_df = self.df[(self.df["start"].between(merged_comb_df.at[i, "start"], merged_comb_df.at[i, "end"]))
                                 & (self.df["seqid"] == comb[0]) & (self.df["strand"] == comb[1])]
                for cc in combine_columns:
                    merged_comb_df.at[i, cc] = '|'.join(set(tmp_df[cc].astype(str).tolist()))
            out_df = out_df.append(merged_comb_df, ignore_index=True)
            out_df.reset_index(inplace=True, drop=True)
            out_df["len"] = out_df["end"] - out_df["start"] + 1
        if self.args.merge_length_violation == "remove":
            out_df = out_df[out_df["len"] <= self.args.max_len]
        elif self.args.merge_length_violation == "merge":
            pass
        else:
            print(f"Unrecognized choice {self.args.merge_length_violation} for --merge_length_violation argument")
        out_df.drop(["len"], inplace=True, axis=1)
        return out_df