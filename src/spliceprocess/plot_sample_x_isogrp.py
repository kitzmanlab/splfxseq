import pandas as pd
from splshared.ssshared import *
import altair as alt

def make_sample_isogrp_heatmap(
    sample_tbl: pd.DataFrame,
):
    
    tblj = []
    for _, r in sample_tbl.iterrows():
        tblj.append( pd.read_table( r['per_samp_isogrp_rpt'] ) )
    tblj = pd.concat( tblj ).reset_index(drop=True)

    lpl = []

    for seqname, subtbl in tblj.groupby('seq_name'):
        base = alt.Chart( 
            subtbl
        ).properties(
            width=100 + 25 * len(subtbl['libname'].unique()),
        )

        loc = subtbl.columns.tolist()

        lisos_sorted = subtbl.groupby('isogrp_name')['isogrp_psi_bc'].mean().sort_values(ascending=False).index.tolist()

        hm = base.mark_rect().encode(
            alt.Y('isogrp_name:O', sort=lisos_sorted),
            alt.X('libname:O' ),
            color=alt.Color('isogrp_psi_bc:Q'),
        )

        hmt = base.mark_text().encode(
            alt.Y('isogrp_name:O', sort=lisos_sorted),
            alt.X('libname:O' ),
            text=alt.Text('isogrp_psi_bc:Q',format='.2f'),
            tooltip=alt.Tooltip(loc)
        )

        lpl.append( hm + hmt )
    
    plj = alt.vconcat(*lpl)

    return plj, tblj

def main():
    import argparse

    parser = argparse.ArgumentParser(description='load in per-sample isogrp tables, and make a sample x isogrp heatmap plot')

    parser.add_argument('--sample_tbl',  dest='sample_tbl' )
    parser.add_argument('--out_plot',  dest='out_plot' )
    parser.add_argument('--out_joined_tbl',  dest='out_joined_tbl' )

    args = parser.parse_args()

    p, tj = make_sample_isogrp_heatmap( pd.read_table( args.sample_tbl ) )
    p.save( args.out_plot )
    tj.to_csv( args.out_joined_tbl, sep='\t', index=False )
    
if __name__ == '__main__':
    main()