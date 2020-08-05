import altair as alt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as tfrms
from collections import Counter
from .coords import pos_to_hgvspos
import splanl.merge_bcs as mbcs
import splanl.post_processing as pp


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

    bc_ranks[col_read_counts+'_rank'] = np.arange(bc_ranks.shape[0])
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
                 shade_exons,
                 gene_name,
                 fig_size = (20,7),
                 coords = 'cdna',
                vec_corange_cloned = None,
                vec_corange_exons = None,
                cdna_corange_exons = None,
                 zoom=None,
                 tick_spacing=10,
                 legend_loc='best',
                 y_ax_title='',
                 legend_title='Nucleotide Substitution from WT',
                 y_ax_lim = (0,1),
                 invert_x = False,
                 tight = True,
                 print_ex_count = False,
                 scale_bar = False,
                 rev_trans = False,
                 hlines = None
                ):

    tbv=vardf.copy()

    tbv.sort_values( by = ['pos'], inplace = True )
    bcs_tbl = tbv.pivot( index = 'pos', columns = 'alt', values = col_y_isoform )

    if zoom:
        assert zoom[1] > zoom[0], 'Your final zoom coordinate must be larger than the first zoom coordinate'
        bcs_tbl = bcs_tbl.loc[ zoom[0]:zoom[1] ]

    #check if plot will lose bars due to not enough space for all the pixels
    dpi = 100
    if ( bcs_tbl.shape[0]*6 ) > ( fig_size[0]*dpi ):
        fig_height = fig_size[1]
        fig_size = ( ( bcs_tbl.shape[0]*6 / dpi ), fig_height )
        print('Adjusting figure width to accomodate all pixels...')

    if print_ex_count:
        print('This figure shows %i exonic bases.' % sum( c[1]-c[0] for c in shade_exons ) )

    #usually want to represent alternate as seen on the forward strand so reverse complement columns
    if rev_trans:
        bcs_tbl = bcs_tbl.rename( columns = { "A": "T", "C": "G", "G": "C", "T": "A" } )

    bcs_tbl.loc[:,['A', 'C', 'G', 'T']].plot.bar( color = [ '#C00001', '#00AD4F', '#FFCF07', '#002966' ],
                                                 align='center',
                                                 width=1,
                                                figsize= fig_size )

    plt.title(
        col_y_isoform+' '+y_ax_title+' by Position for Single Nucleotide Variants in $\it{%s}$'%gene_name,
        fontsize=24)

    if y_ax_lim:
        plt.ylim( y_ax_lim )

    plt.ylabel(col_y_isoform+' '+y_ax_title,fontsize=22)
    plt.yticks(fontsize=18)

    if coords.lower() == 'cdna':
        plt.xlabel('cDNA Position',fontsize=22)
        bcs_tbl['hgvs_pos'] = pos_to_hgvspos( bcs_tbl.index,
                               vec_corange_cloned,
                               vec_corange_exons,
                               cdna_corange_exons
                             )
        plt.xticks( [idx for idx,p in enumerate(bcs_tbl.index) if idx%tick_spacing==0],
                       [c for idx,c in enumerate(bcs_tbl.hgvs_pos) if idx%tick_spacing==0],
                       fontsize=18,
                       rotation='vertical' )

    elif coords.lower() == 'gdna':
        plt.xlabel('gDNA Position',fontsize=22)
        plt.xticks( [idx for idx,p in enumerate(bcs_tbl.index) if idx%tick_spacing==0],
                   [c for idx,c in enumerate(bcs_tbl.index) if idx%tick_spacing==0],
                   fontsize=18,
                   rotation='vertical' )

    elif coords.lower() == 'vector':
        plt.xlabel('Vector Position',fontsize=22)
        plt.xticks( [idx for idx,p in enumerate(bcs_tbl.index) if idx%tick_spacing==0],
                   [c for idx,c in enumerate(bcs_tbl.index) if idx%tick_spacing==0],
                   fontsize=18,
                   rotation='vertical' )

    if hlines:
        for line in hlines:
            plt.axhline( line, c = 'black', ls = '--', alpha = .6 )

    legend = plt.legend( title = legend_title,
                         ncol = 2,
                         loc = legend_loc,
                         fontsize = 14)
    plt.setp(legend.get_title(),fontsize=14)

    for ex in shade_exons:
        plt.axvspan( bcs_tbl.index.get_loc( ex[0] ) - .5,
                    bcs_tbl.index.get_loc( ex[1] ) + .5,
                    facecolor = 'gray',
                    alpha = 0.15)

    if invert_x:
        plt.gca().invert_xaxis()

    if scale_bar:
        ax = plt.gca()
        trans = tfrms.blended_transform_factory( ax.transData, ax.transAxes )
        plt.errorbar( tick_spacing, 0.96, xerr=tick_spacing/2, color='black', capsize=3, transform=trans)
        plt.text( tick_spacing, 0.94, str( tick_spacing )+' bases',  horizontalalignment='center',
        verticalalignment='top', transform=trans, fontsize = 14 )

    if tight:
        plt.tight_layout()

    plt.show()

def PlotBCsByPos(vardf,
                 vec_corange_cloned,
                 vec_corange_exons,
                 cdna_corange_exons,
                 shade_exons,
                 gene_name,
                 y_ax_lim=None,
                 fig_size = (20,7) ):

    tbv=vardf.copy()

    tbv.sort_values(by=['pos'],inplace=True)
    bcs_tbl=tbv.pivot(index='pos',columns='alt',values='n_bc_passfilt')
    bcs_tbl['hgvs_pos'] = mbcs.pos_to_hgvspos( bcs_tbl.index,
                           vec_corange_cloned,
                           vec_corange_exons,
                           cdna_corange_exons
                         )

    #check if plot will lose bars due to not enough space for all the pixels
    dpi = 100
    if ( bcs_tbl.shape[0]*1.5 ) > ( fig_size[0]*dpi ):
        fig_height = fig_size[1]
        fig_size = ( ( bcs_tbl.shape[0]*1.5 / dpi ), fig_height )
        print('Adjusting figure width to accomodate all pixels...')

    bcs_tbl.loc[:,['A','C', 'G', 'T']].plot.bar( stacked = True,
                                                color = [ '#C00001', '#00AD4F', '#FFCF07', '#002966' ],
                                                figsize = fig_size )

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

    #counts entries that are missing - one is missing per position as its reference
    missing = sum( bcs_tbl.isnull().sum() ) - bcs_tbl.shape[0]
    #counts number of entries with 0 barcodes passing the filter
    zero = 4*bcs_tbl.shape[0] - sum( bcs_tbl[[ 'A', 'C', 'G', 'T' ]].astype( bool ).sum() )
    print( 100 * ( 1- ( ( missing + zero ) / (3*bcs_tbl.shape[0] ) ) ), '% of all possible mutations present' )

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

def barplot_allisos( allisos_df,
                    isoform_df,
                    psi_cols,
                    stat,
                    title = ''):

    plot_df = allisos_df.copy()
    iso_df = isoform_df.copy()

    if stat == 'max':
        plot_df[ psi_cols ].max().plot( kind = 'bar',
                                        figsize = ( 40, 5 ),
                                        title = title,
                                        ylim = ( 0, 1 ) )

    if stat == 'min':
        plot_df[ psi_cols ].min().plot( kind = 'bar',
                                        figsize = ( 40, 5 ),
                                        title = title,
                                        ylim = ( 0, 1 ) )

    if stat == 'mean':
        plot_df[ psi_cols ].mean().plot( kind = 'bar',
                                        figsize = ( 40, 5 ),
                                        title = title,
                                        ylim = ( 0, 1 ) )

    loc, labels = plt.xticks()

    plt.xticks( loc,
                iso_df.loc[ [ iso.split('_')[1] for iso in psi_cols ] ].isoform.tolist(),\
               fontsize=14,
               rotation='vertical' )

plt.show()

def barplot_across_samples( byvartbl_long,
                            bar_col,
                            y_scale = None,
                            y_label = None,
                            title = '',
                            y_lim = None,
                            color_col = None,
                            color_dict = None):

    tbv = byvartbl_long.copy()

    assert ( color_col and color_dict ) or not( color_col or color_dict ), \
    "color_col and color_dict must either both be none or both entered"

    if color_col:
        tbv.set_index('sample')[ bar_col ].plot.bar( color=[ color_dict[i] for i in tbv[ color_col ] ] )
    else:
        tbv.set_index('sample')[ bar_col ].plot.bar()

    if y_label:
        plt.ylabel( y_label )

    if y_lim:
        plt.ylim( y_lim )

    if y_scale:
        plt.yscale( y_scale )

    plt.title( title )

    plt.show()

def per_sdv_by_thresh( byvartbl,
                     sdv_col,
                     thresh_range = None,
                     abs_vals = True,
                     num_pts = 20,
                     title = '',
                     y_lim = None,
                     vlines = ( 1.96, 3 ),
                     fig_size = ( 12, 9.5 ) ):

    tbv = byvartbl.loc[ byvartbl.n_bc_passfilt > 0  ].copy()

    samp = 'sample' in tbv.columns

    if thresh_range:
        thresh = np.arange( thresh_range[0],
                            thresh_range[1],
                            ( thresh_range[1] - thresh_range[0] ) / num_pts )
    else:
        min_thresh = 0 if tbv[ sdv_col ].min() <= 0 else tbv[ sdv_col ].min()
        thresh = np.arange( min_thresh,
                            tbv[ sdv_col ].max(),
                            ( tbv[ sdv_col ].max() - min_thresh ) / num_pts )

    #add vline points into array while maintaining sort order
    for pt in vlines:
        thresh = np.insert( thresh, thresh.searchsorted( float( pt ) ), float( pt ) )

    plt.figure( figsize = fig_size )

    if samp:
        for smpl in set( tbv['sample'] ):
            sdv_per = []
            for t in thresh:
                tbv_samp = tbv.query( 'sample == "%s"' % smpl ).copy()
                tbv_samp = pp.sdvs( tbv_samp, sdv_col, t, abs_vals = abs_vals )
                sdv_per.append( 100*( tbv_samp.sdv.sum() / tbv_samp.shape[0] ) )

            plt.plot( thresh, sdv_per, label = smpl )

            print( smpl )
            for pt in vlines:
                per = sdv_per[ np.where( thresh == pt )[0][0] ]
                print( 'At a threshold of %.2f, %.2f%% of variants are splice disrupting.' % ( pt, per ) )

        plt.legend()

    else:
        sdv_per = []
        for t in thresh:
            tbv = pp.sdvs( tbv.query( 'sample == "%s"' % smpl ), sdv_col, t, abs_vals = abs_vals )
            sdv_per.append( 100*( tbv.sdv.sum() / tbv.shape[0] ) )

        plt.plot( thresh, sdv_per, label = smpl )

        for pt in vlines:
            per = sdv_per[ np.where( thresh == pt ) ]
            print( 'At a threshold of %.2f, %.2f%% of variants are splice disrupting.' % ( pt, per ) )

    for pt in vlines:
        plt.axvline( x = pt, color = 'gray', linestyle = '--' )

    if y_lim:
        plt.ylim( y_lim )
    else:
        plt.ylim( ( 0, 100 ) )

    plt.ylabel( 'Percentage splice disrupting variants' )
    plt.xlabel( sdv_col + ' threshold' )

    plt.title( title )

    plt.show()

def per_repeat_sdv( byvartbl,
                     sdv_col = None,
                     thresh = None,
                     abs_vals = True,
                     title = '',
                     y_lim = None, ):

    tbv = byvartbl.loc[ byvartbl.n_bc_passfilt > 0  ].copy()

    n_var = len( set( tbv.varlist ) )

    assert ( sdv_col and thresh ) or ( 'sdv' in tbv.columns ), \
    'Please specify a column and threshold to determine splice disrupting variants'

    if ( sdv_col and thresh ):
        tbv = pp.sdvs( tbv, sdv_col, thresh, abs_vals = abs_vals )

    n_samp = len( set( tbv[ 'sample' ] ) )

    assert n_samp > 1, 'Please provide a dataset with more than one sample'

    sdvs = tbv.loc[ tbv.sdv ]

    #this gets you a counter { n_samples: n_sdvs }
    repeat_counts = Counter( Counter( sdvs.varlist ).values() )

    n_sdvs = sum( repeat_counts.values() )

    print( '%i of %i variants (%.2f%%) are splice disrupting in at least one sample.' \
            % ( n_sdvs, n_var, 100*( n_sdvs / n_var )  ) )

    #adds zero counts to create plot with full range of possible values
    for i in range( 1, n_samp + 1 ):
        if i not in repeat_counts:
            repeat_counts[ i ] = 0

    print( '%i of %i variants (%.2f%%) are splice disrupting in all samples.' \
            % ( repeat_counts[ n_samp ], n_var, 100*( repeat_counts[ n_samp ] / n_var )  ) )

    #gives us a list sorted by n_samples in ascending order
    repeat_counts = sorted( repeat_counts.items() )

    labels, counts = zip( *repeat_counts )

    per = [ 100*( count / n_sdvs ) for count in counts ]

    indexes = np.arange( n_samp )

    plt.bar( indexes, per )

    plt.xticks(indexes, labels)

    if y_lim:
        plt.ylim( y_lim )
    else:
        plt.ylim( ( 0, 100 ) )

    plt.ylabel( 'Percent of splice disrupting variants' )
    plt.xlabel( 'Number of samples' )

    plt.title( title )

    plt.show()
