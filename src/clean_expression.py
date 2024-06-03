#!/usr/bin/env python3
import polars as pl
import polars.selectors as cs
import argparse, sys
from motif_mining.utils.misc import lowercase_cols

def cline(args_ls = sys.argv[1:]):
    p = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = '',
    )
    p.add_argument('infile', help = 'path to expression file')
    p.add_argument('outfile', help = 'path to write out')
    p.add_argument('--meta_out', help = 'write metadata table to path, if supplied')
    p.add_argument('--obs_id', default='obs_id',
                   help = 'Create obs_id column of this name by joining group and subject id cols')
    p.add_argument('--obs_id_sep', default='_', help = 'Separator to use in creating obs_id')
    p.add_argument('--group_id', default='visit')
    p.add_argument('--rename_group_id', action='store_true', help = 'If supplied, rename group_id to "group"')
    p.add_argument('--subject_id', default='subject_id')
    p.add_argument('--var_id', default='symbol', help = '')
    p.add_argument('--value', default='npx', help = '')
    args = p.parse_args(args_ls)
    return args

def main(args):
    id_cols = [args.subject_id, args.group_id]
    npx = pl.read_csv(args.infile, null_values=['NA', 'na', 'none', 'None'])
    npx = npx.rename({x: x.lower().replace(' ', '_') for x in npx.columns})
    if args.rename_group_id:
        npx = npx.rename({args.group_id: 'group'})
        id_cols = [x.replace(args.group_id, 'group') for x in id_cols]
    npx = npx.cast({col: pl.Utf8 for col in id_cols})
    npx = npx.melt(id_vars=id_cols, variable_name=args.var_id, value_name=args.value)
    npx = npx.with_columns(
        pl.concat_str(id_cols, separator=args.obs_id_sep,).alias(args.obs_id),
        pl.col(args.var_id).str.to_uppercase() # Symbol (also all var_ids, off the top of my head) expect UPPER
    )
    npx.write_csv(args.outfile)

    if args.meta_out:
        meta = npx.select(pl.exclude([args.value, args.var_id])).unique()
        meta.write_csv(args.meta_out)

if __name__ == "__main__":
    main(cline())
