import altair as alt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .coords import pos_to_hgvspos
import splanl.merge_bcs as mbcs


def waterfall_plot(
    bctbl,
    col_read_counts,
    percentile,
    title='',
    x_ax_title='Read Count Rank Log10',
    y1_ax_title='Read Count Log10',
    y2_ax_title='Cumulative Read Count Percentile'
):

    bc_ranks = bctbl.copy()

    bc_ranks.sort_values(by=[col_read_counts],ascending=False,inplace=True)

    bc_ranks[col_read_counts+'_rank'] =  np.arange(bc_ranks.shape[0])
    bc_ranks[col_read_counts+'_log10']=np.log10(bc_ranks[col_read_counts]+.1)
    bc_ranks[col_read_counts+'_rank_log10']=np.log10(bc_ranks[col_read_counts+'_rank']+1)
    bc_ranks['cumulative_read_percentage'] = 100*bc_ranks[col_read_counts].cumsum()/bc_ranks[col_read_counts].sum()
    percentile_cutoff = bc_ranks.loc[bc_ranks['cumulative_read_percentage']>=percentile][col_read_counts+'_rank_log10'].reset_index(drop=True)[0]
    read_count_cutoff = bc_ranks.loc[bc_ranks['cumulative_read_percentage']>=percentile][col_read_counts].reset_index(drop=True)[0]

    fig, ax1 = plt.subplots()

    color = 'blue'
    ax1.set_xlabel(x_ax_title)
    ax1.set_ylabel(y1_ax_title, color=color)
    ax1.plot(bc_ranks[col_read_counts+'_rank_log10'], bc_ranks[col_read_counts+'_log10'], color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'red'
    ax2.set_ylabel(y2_ax_title, color=color)  # we already handled the x-label with ax1
    ax2.plot(bc_ranks[col_read_counts+'_rank_log10'], bc_ranks['cumulative_read_percentage'], color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    plt.axvline(percentile_cutoff, color='black')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    print('The read count cut off at the',percentile,'th percentile is',read_count_cutoff)


##########################################################################################
##########################################################################################
##########################################################################################
##
##  psi by coordinate plots
##


##########################################################################################
## altair-style, one sample

def altplot_psi_by_pos_onesamp(
    tblbyvar,
    col_iso_y,
    addl_cols_tooltip,
    fill_col='alt',
    lposrng_exon=None,
    yscale=None,
    height=100,
    width=800,
    title='',
    x_ax_title='pos',
    y_ax_title=None
):
    assert (tblbyvar.ref.str.len() == 1).all()
    assert (tblbyvar.alt.str.len() == 1).all()

    tbv = tblbyvar.copy()
    tbv['pos_display'] = tbv['pos']
    for ofs,mut in zip( (0.1,0.3,0.5,0.7), 'ACGT' ):
        tbv.loc[ tbv.alt==mut, 'pos_display' ]+=ofs

    if yscale is None:
        y = alt.Y('{}:Q'.format(col_iso_y), axis=alt.Axis(title=y_ax_title) )
    else:
        y = alt.Y('{}:Q'.format(col_iso_y), scale=yscale, axis=alt.Axis(title=y_ax_title) )

    if y_ax_title is None:
        y_ax_title=col_iso_y

    points = alt.Chart( tbv , title=title ).mark_circle(
        size=15
    ).encode(
        x=alt.X('pos_display:Q', scale=alt.Scale(zero=False) ,
        axis=alt.Axis(title=x_ax_title)),#, sort=list( pos_display_categ.categories )),
        y=y,
        tooltip=['pos','ref','alt']+addl_cols_tooltip,
        fill=alt.Color(fill_col+':N',scale=alt.Scale(scheme='category10'))
    ).properties( height=height, width=width )

    gr = points

    if lposrng_exon is not None:
        exonbox=pd.DataFrame( {'x':[lposrng_exon[i][0] for i in range(len(lposrng_exon))],
                              'x2':[lposrng_exon[i][1] for i in range(len(lposrng_exon))],
                             'y':tblbyvar[col_iso_y].min(), 'y2':tblbyvar[col_iso_y].max()} )
        ex = alt.Chart( exonbox ).mark_rect( color='yellow', opacity=0.1 ).encode( x='x', x2='x2', y='y', y2='y2' )
        gr = ex + gr

    return gr

##########################################################################################
## altair-style, one sample w/ zoom

def altplot_psi_by_pos_onesamp_wzoom(
    tblbyvar,
    col_iso_y,
    addl_cols_tooltip,
    fill_col='alt',
    lposrng_exon=None,
    yscale=None,
    height=100,
    width=800,
    zoom_frac = 0.5,
    title=''
):

    zoom_frac = max(0.1, min(0.9, zoom_frac))

    assert (tblbyvar.ref.str.len() == 1).all()
    assert (tblbyvar.alt.str.len() == 1).all()

    tbv = tblbyvar.copy()
    tbv['pos_display'] = tbv['pos']
    for ofs,mut in zip( (0.1,0.3,0.5,0.7), 'ACGT' ):
        tbv.loc[ tbv.alt==mut, 'pos_display' ]+=ofs

    if yscale is None:
        y = alt.Y('{}:Q'.format(col_iso_y))
    else:
        y = alt.Y('{}:Q'.format(col_iso_y), scale=yscale )

    selbrush = alt.selection(type='interval', encodings=['x'])

    points = alt.Chart( tbv, title=title ).mark_point(
        size=15
    ).encode(
        x=alt.X('pos_display:Q', scale=alt.Scale(zero=False, domain=selbrush)),#, sort=list( pos_display_categ.categories )),
        y=y,
        tooltip=['pos','ref','alt']+addl_cols_tooltip,
        # fill='alt:N'
        fill=alt.Color(fill_col+':N',scale=alt.Scale(scheme='category10'))
    ).properties( height=height, width=int(1-zoom_frac)*width )

    gr = points

    if lposrng_exon is not None:
        exonbox=pd.DataFrame( {'x':[lposrng_exon[i][0] for i in range(len(lposrng_exon))],
                              'x2':[lposrng_exon[i][1] for i in range(len(lposrng_exon))],
                             'y':np.min(tblbyvar.col_iso_y), 'y2':1} )
        ex = alt.Chart( exonbox ).mark_rect( color='yellow', opacity=0.1 ).encode( x='x', x2='x2', y='y', y2='y2' )
        gr = ex + gr


    pointszoom = points.properties( width=int(zoom_frac)*width ).add_selection( selbrush )

    gr_ret = pointszoom | gr

    return gr_ret

def altplot_psi_by_pos_onesamp_multiiso_wzoom(
    tblbyvar,
    l_col_iso_y,
    addl_cols_tooltip,
    lposrng_exon=None,
    yscale=None,
    height=100,
    width=800,
    zoom_frac = 0.5,
    title='' ):

    zoom_frac = max(0.1, min(0.9, zoom_frac))

    assert (tblbyvar.ref.str.len() == 1).all()
    assert (tblbyvar.alt.str.len() == 1).all()

    tbv = tblbyvar.copy()
    tbv['pos_display'] = tbv['pos']
    for ofs,mut in zip( (0.1,0.3,0.5,0.7), 'ACGT' ):
        tbv.loc[ tbv.alt==mut, 'pos_display' ]+=ofs

    selbrush = alt.selection(type='interval', encodings=['x'])

    # hack to get a title over everything
    tpl = alt.Chart( {'values':[{'text':title}]} ).mark_text( size=14 ).encode( text='text:N' )

    lpoints = [tpl]

    for col_iso_y in l_col_iso_y:
        if yscale is None:
            y = alt.Y('{}:Q'.format(col_iso_y))
        else:
            y = alt.Y('{}:Q'.format(col_iso_y), scale=yscale )

        points = alt.Chart( tbv ).mark_point(
            size=15
        ).encode(
            x=alt.X('pos_display:Q', scale=alt.Scale(zero=False, domain=selbrush)),#, sort=list( pos_display_categ.categories )),
            y=y,
            tooltip=['pos','ref','alt']+addl_cols_tooltip,
            # fill='alt:N'
            fill=alt.Color('alt:N',scale=alt.Scale(scheme='category10'))
        ).properties( height=height, width=int(1-zoom_frac)*width )

        pointszoom = points.properties( width=int(zoom_frac)*width ).add_selection( selbrush )

        if lposrng_exon is not None:
            exonbox=pd.DataFrame( {'x':[lposrng_exon[i][0] for i in range(len(lposrng_exon))],
                                  'x2':[lposrng_exon[i][1] for i in range(len(lposrng_exon))],
                                 'y':0, 'y2':1} )
            ex = alt.Chart( exonbox ).mark_rect( color='yellow', opacity=0.1 ).encode( x='x', x2='x2', y='y', y2='y2' )
            points = ex + points

        lpoints.append( pointszoom | points )

    gr = alt.vconcat( *lpoints ).configure_view( stroke=None ).configure_concat( spacing=1 )

    return gr

def altplot_scatter_onesamp(
    tblbyvar,
    col_iso_y,
    col_x,
    addl_cols_tooltip=[],
    fill_col='var_type',
    yscale=None,
    height=300,
    width=300,
    title=''
):

    tbv = tblbyvar.copy()

    if yscale is None:
        y = alt.Y('{}:Q'.format(col_iso_y))
    else:
        y = alt.Y('{}:Q'.format(col_iso_y), scale=yscale )

    points = alt.Chart( tbv , title=title ).mark_circle(
        size=20
    ).encode(
        x=alt.X(col_x+':Q', scale=alt.Scale(zero=False)),#, sort=list( pos_display_categ.categories )),
        y=y,
        tooltip=addl_cols_tooltip,
        fill=alt.Color(fill_col+':N',scale=alt.Scale(scheme='category10'))
    ).properties( height=height, width=width )

    gr = points

    return gr

def altplot_scatter_nofill(
    tblbyvar,
    col_iso_y,
    col_x,
    addl_cols_tooltip=[],
    yscale=None,
    height=300,
    width=300,
    title=''
):

    tbv = tblbyvar.copy()

    if yscale is None:
        y = alt.Y('{}:Q'.format(col_iso_y))
    else:
        y = alt.Y('{}:Q'.format(col_iso_y), scale=yscale )

    points = alt.Chart( tbv , title=title ).mark_circle(
        size=20
    ).encode(
        x=alt.X(col_x+':Q', scale=alt.Scale(zero=False)),#, sort=list( pos_display_categ.categories )),
        y=y,
        tooltip=addl_cols_tooltip
    ).properties( height=height, width=width )

    gr = points

    return gr

def altplot_violin_onesamp(
    tblbyvar,
    col_y,
    col_cat_x,
    yscale=None,
    height=400,
    width=400,
    title=''
):
    assert (tblbyvar.ref.str.len() == 1).all()
    assert (tblbyvar.alt.str.len() == 1).all()

    tbv = tblbyvar.copy()

    points = alt.Chart( tbv , title=title ).transform_density(
        col_y,
        as_=[col_y,'density'],
        groupby=[col_cat_x]
    ).mark_area(
        orient='horizontal'
    ).encode(
        y=col_y+':Q',
        color=col_cat_x+':N',
        x=alt.X(
            'density:Q',
            stack='center',
            axis=alt.Axis(labels=False, values=[0],grid=False, ticks=True),
        ),
        column=alt.Column(
            col_cat_x+':N',
            header=alt.Header(
                titleOrient='bottom',
                labelOrient='bottom',
                labelPadding=0,
            )
        ),
    ).properties( height=height, width=width )

    gr = points

    return gr

def PlotPSIByPos(vardf,
                 col_y_isoform,
                 vec_corange_cloned,
                 vec_corange_exons,
                 cdna_corange_exons,
                 shade_exons,
                 gene_name,
                 zoom=None,
                 tick_spacing=10,
                 legend_loc='best',
                 y_ax_title='PSI',
                 legend_title='Nucleotide Substitution from WT',
                 y_ax_lim=None
                ):

    tbv=vardf.copy()

    tbv.sort_values(by=['pos'],inplace=True)
    bcs_tbl=tbv.pivot(index='pos',columns='alt',values=col_y_isoform)
    bcs_tbl['hgvs_pos'] = pos_to_hgvspos( bcs_tbl.index,
                           vec_corange_cloned,
                           vec_corange_exons,
                           cdna_corange_exons
                         )

    if zoom:
        bcs_tbl=bcs_tbl.loc[ zoom[0]:zoom[1] ]

    bcs_tbl.loc[:,['A','C', 'G', 'T']].plot.bar(color=['#C00001', '#00AD4F', '#FFCF07', '#002966'],
                                                align='center',
                                                width=1,
                                                figsize=(20,7))

    plt.title(
        col_y_isoform+' '+y_ax_title+' by Position for Single Nucleotide Variants in $\it{%s}$'%gene_name,
        fontsize=24)

    if y_ax_lim:
        plt.ylim(0,y_ax_lim)

    plt.ylabel(col_y_isoform+' '+y_ax_title,fontsize=22)
    plt.yticks(fontsize=18)

    plt.xlabel('cDNA Position',fontsize=22)
    plt.xticks( [idx for idx,p in enumerate(bcs_tbl.index) if idx%tick_spacing==0],
                   [c for idx,c in enumerate(bcs_tbl.hgvs_pos) if idx%tick_spacing==0],
                   fontsize=18,
                   rotation='vertical' )

    legend = plt.legend(title=legend_title,
                        ncol=2,
                        loc=legend_loc,
                        fontsize=14)
    plt.setp(legend.get_title(),fontsize=14)

    for ex in shade_exons:
        plt.axvspan(bcs_tbl.index.get_loc(ex[0])-.5,
                    bcs_tbl.index.get_loc(ex[1])+.5,
                    facecolor='gray',
                    alpha=0.15)

    plt.tight_layout()
    plt.show()

def PlotBCsByPos(vardf,
                 vec_corange_cloned,
                 vec_corange_exons,
                 cdna_corange_exons,
                 shade_exons,
                 gene_name,
                 y_ax_lim=None):

    tbv=vardf.copy()

    tbv.sort_values(by=['pos'],inplace=True)
    bcs_tbl=tbv.pivot(index='pos',columns='alt',values='n_bc')
    bcs_tbl['hgvs_pos'] = mbcs.pos_to_hgvspos( bcs_tbl.index,
                           vec_corange_cloned,
                           vec_corange_exons,
                           cdna_corange_exons
                         )

    bcs_tbl.loc[:,['A','C', 'G', 'T']].plot.bar(stacked=True,
                                                color=['#C00001', '#00AD4F', '#FFCF07', '#002966'],
                                                figsize=(20,7))

    plt.title(
        'Number of Distinct Barcodes Present by Position for Single Nucleotide Variants in $\it{%s}$'%gene_name,
        fontsize=24)

    if y_ax_lim:
        plt.ylim(0,y_ax_lim)

    plt.ylabel('Number of Distinct Barcodes',fontsize=22)
    plt.yticks(fontsize=18)

    plt.xlabel('cDNA Position',fontsize=22)
    plt.xticks( [idx for idx,p in enumerate(bcs_tbl.index) if idx%10==0],
               [c for idx,c in enumerate(bcs_tbl.hgvs_pos) if idx%10==0],
               fontsize=18,
               rotation='vertical' )

    legend = plt.legend(title='Nucleotide Substitution from WT',
                        ncol=2,
                        loc='upper left',
                        fontsize=14)
    plt.setp(legend.get_title(),fontsize=14)

    for ex in shade_exons:
        plt.axvspan(bcs_tbl.index.get_loc(ex[0])-.5,
                    bcs_tbl.index.get_loc(ex[1])+.5,
                    facecolor='gray',
                    alpha=0.15)

    plt.tight_layout()
    plt.show()
    bcs_tbl['num_missing']=bcs_tbl.isnull().sum(axis=1)
    present = sum( bcs_tbl.reset_index().isna().sum() ) - bcs_tbl.shape[0]
    #subtraction accounts for missing values due to reference allele at each position
    print( str( 100 * ( 1- ( present / (3*bcs_tbl.shape[0] ) ) ) ) +'% of all possible mutations present' )

def plot_corr_waterfalls(benchmark_df,
                        compare_df_long,
                        corr_col,
                        benchmark_samp_name,
                        merge_idx = ['chrom','pos','ref','alt','varlist'],
                        sample_col = 'sample'):

    rank_df = benchmark_df.copy()[ merge_idx + [corr_col] ]
    rank_df[ corr_col+'_rank' ] =  np.argsort( np.argsort( np.array( -rank_df[ corr_col ] ) ) )
    rank_df = rank_df.set_index( merge_idx )

    compare_df = compare_df_long.copy().set_index( merge_idx )

    merge_df = compare_df.merge( rank_df[ corr_col+'_rank' ], left_index=True, right_index=True )

    for samp in list( set( merge_df[ sample_col ] ) ):
        print(samp)
        plt.scatter(merge_df.query('sample == "%s"' %samp)[ corr_col+'_rank' ],
                    merge_df.query('sample == "%s"' %samp)[ corr_col ],
                    s=5)

        plt.ylabel(samp + ' ' + corr_col)
        plt.xlabel( benchmark_samp_name + ' ' + corr_col + ' rank')
        plt.show()
