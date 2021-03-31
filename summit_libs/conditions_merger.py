import pandas as pd


class ConditionsMerger:

    def __init__(self, df, args):
        self.df = df
        self.args = args

    def merge(self):
        print("Merging annotations of all conditions")
        df_cols = self.df.columns.tolist()
        out_df = pd.DataFrame(columns=df_cols)
        unmerged_df = pd.DataFrame(columns=df_cols)
        self.df.sort_values(["seqid", "strand", "start", "end"], inplace=True, ascending=[True, True, True, False])
        self.df.reset_index(inplace=True, drop=True)
        out_df = out_df.append(self.df.iloc[0])
        interval_col_id = out_df.columns.get_loc("interval")
        start_col_id = out_df.columns.get_loc("start")
        end_col_id = out_df.columns.get_loc("end")
        none_merge_cols = ["seqid", "strand", "start", "end", "interval"]
        merge_cols_ids = [x for x in df_cols if x not in none_merge_cols]
        for col in merge_cols_ids:
            out_df[col] = out_df[col].astype(str)
        for i in self.df.index:
            if self.df.at[i, "interval"].overlaps(out_df.iat[-1, interval_col_id]):
                min_start = min([self.df.at[i, "start"], out_df.iat[-1, start_col_id]])
                max_end = max([self.df.at[i, "end"], out_df.iat[-1, end_col_id]])
                if max_end - min_start + 1 > self.args.max_len:
                    if self.args.merge_length_violation != "no_merge":
                        unmerged_df = unmerged_df.append(self.df.loc[i])
                        self.df.drop(i, inplace=True)
                        continue
                    elif self.args.merge_length_violation != "remove":
                        out_df = out_df.iloc[:-1]
                        self.df.drop(i, inplace=True)
                        continue
                    else:
                        pass
                out_df.iat[-1, start_col_id] = min_start
                out_df.iat[-1, end_col_id] = max_end
                out_df.iat[-1, interval_col_id] = pd.Interval(left=min_start, right=max_end)
                for col in merge_cols_ids:
                    col_id = out_df.columns.get_loc(col)
                    out_df.iat[-1, col_id] = f"{str(out_df.iat[-1, col_id])},{str(self.df.at[i, col])}"
            else:
                out_df = out_df.append(self.df.loc[i])
        return out_df, unmerged_df