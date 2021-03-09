import altair as alt
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.transforms as tfrms
from matplotlib.artist import Artist as art
from matplotlib import gridspec
from collections import Counter
import splanl.coords as cds
import splanl.merge_bcs as mbcs
import splanl.post_processing as pp
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve
import upsetplot as up

github_colors = '3182bd6baed69ecae1c6dbefe6550dfd8d3cfdae6bfdd0a231a35474c476a1d99bc7e9c0756bb19e9ac8bcbddcdadaeb636363969696bdbdbdd9d9d9'
light_colors = [ '#' + github_colors[i:i+6] for i in range( 0, len( github_colors ), 6 ) ]

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
        bcs_tbl['hgvs_pos'] = cds.pos_to_hgvspos( bcs_tbl.index,
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
    bcs_tbl['hgvs_pos'] = cds.pos_to_hgvspos( bcs_tbl.index,
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

def plot_stacked_psi(    var_df,
                          zcols,
                          pos_col,
                          colors,
                          color_col = 'alt',
                          fig_size = ( 20, 5 ),
                          shade_exons = False,
                          zoom=None,
                          tick_spacing = 10,
                          alt_labels = False,
                          bar_labels = False,
                          title = '',
                          y_ax_title='',
                          y_ax_lim = None,
                          x_ax_title = '',
                          legend = True,
                          legend_title = '',
                          legend_loc = 'best',
                          legend_labels = None,
                          tight = True,
                          print_ex_count = False,
                          scale_bar = False,
                          rev_trans = False,
                          hlines = None,
                          savefile = None ):

    tbv = var_df.sort_values( by = [ 'pos', 'alt' ] ).copy()

    if zoom:
        assert zoom[1] > zoom[0], 'Your final zoom coordinate must be larger than the first zoom coordinate'
        tbv_filt = tbv.set_index( 'pos' ).loc[ zoom[0]:zoom[1] ].reset_index()
    else:
        tbv_filt = tbv

    #check if plot will lose bars due to not enough space for all the pixels
    dpi = 100
    if ( tbv_filt.shape[0]*1.5 ) > ( fig_size[0]*dpi ):
        fig_height = fig_size[1]
        fig_size = ( ( tbv_filt.shape[0]*1.5 / dpi ), fig_height )
        print('Adjusting figure width to accomodate all pixels...')

    if print_ex_count:
        print('This figure shows %i exonic bases.' % sum( c[1]-c[0] for c in shade_exons ) )

    col_pivot = tbv_filt.set_index( ['pos', 'alt'] )[ zcols ]

    col_pivot.plot.bar( color = colors,
                        stacked = True,
                        align = 'center',
                        width = 1,
                        figsize = fig_size )

    ax = plt.gca()

    if alt_labels or bar_labels:

        rects = ax.patches

        heights = [ rect.get_height() for rect in rects ]

        n_bars = int( len( heights ) / len( zcols ) )

        pos_ht_sum = [ sum( heights[ i ] for i in range( j, len( heights ), n_bars ) if heights[ i ] > 0 )
                       for j in range( n_bars ) ]
        neg_ht_sum = [ sum( heights[ i ] for i in range( j, len( heights ), n_bars ) if heights[ i ] < 0 )
                       for j in range( n_bars ) ]

    if alt_labels:

        labels = [ idx[1] for idx in col_pivot.index ]

        min_ht = min( neg_ht_sum )

        for rect, label in zip( rects, labels ):
            height = rect.get_height()
            ax.text( rect.get_x() + rect.get_width() / 2,
                     min_ht - 1.5,
                     label,
                     fontsize = 14,
                     fontweight = 'bold',
                     ha='center',
                     va='bottom' )

    if bar_labels:

        for iidx, colmark in enumerate( bar_labels ):

            col, marker = colmark

            #true false vector for where to put the marker
            locs = tbv_filt[ col ]

            for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc:

                    ax.text( rect.get_x() + rect.get_width() /2,
                             pos_ht_sum[ jidx ] + .1 + .9*iidx,
                             marker,
                             fontsize = 10,
                             ha = 'center',
                             va = 'bottom' )

    plt.title( title, fontsize = 24 )

    if not( y_ax_lim ) and ( bar_labels or alt_labels ):

        col_pivot[ 'y_max' ] = col_pivot[ col_pivot > 0 ].sum( axis = 1 )
        col_pivot[ 'y_min' ] = col_pivot[ col_pivot < 0 ].sum( axis = 1 )

        if bar_labels and alt_labels:
            y_ax_lim = ( col_pivot.y_min.min() - 1.6,
                         col_pivot.y_max.max() + .1 + .9*len( bar_labels ) )
        elif bar_labels:
            y_ax_lim = ( col_pivot.y_min.min()*1.01,
                         col_pivot.y_max.max() + .1 + .9*len( bar_labels ) )
        else:
            y_ax_lim = ( col_pivot.y_min.min() - 1.6,
                         col_pivot.y_max.max()*1.01 )

    if y_ax_lim:
        plt.ylim( y_ax_lim )

    plt.ylabel( y_ax_title, fontsize = 22 )
    plt.yticks( fontsize = 18 )

    plt.xlabel( x_ax_title, fontsize = 22 )
    plt.xticks( [ idx for idx,p in enumerate( tbv_filt.index ) if idx%( 3*tick_spacing ) == 1 ],
                [ c for idx,c in enumerate( tbv_filt[ pos_col ] ) if idx%( 3*tick_spacing ) == 1 ],
                fontsize=18,
                rotation='vertical' )

    if hlines:
        for line in hlines:
            plt.axhline( line, c = 'black', ls = '--', alpha = .6 )

    if legend:

        if legend_labels:
            legend = plt.legend( title = legend_title,
                                 ncol = 2,
                                 loc = legend_loc,
                                 labels = legend_labels,
                                 fontsize = 14 )
        else:
            legend = plt.legend( title = legend_title,
                                 ncol = 2,
                                 loc = legend_loc,
                                 fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )
    else:
        ax.legend_ = None
        plt.draw()

    if shade_exons:
        for ex in shade_exons:
            plt.axvspan( col_pivot.index.get_loc( ex[0] ).start - .5,
                         col_pivot.index.get_loc( ex[1] ).stop,
                         facecolor = 'gray',
                         alpha = 0.15 )

    if scale_bar:
        trans = tfrms.blended_transform_factory( ax.transData, ax.transAxes )
        plt.errorbar( 3*tick_spacing, 0.96, xerr = ( 3*tick_spacing / 2 ), color='black', capsize=3, transform=trans)
        txt = str( tick_spacing ) + ' bases' if tick_spacing > 1 else str( tick_spacing ) + ' base'
        plt.text( 3*tick_spacing, 0.94, txt,  horizontalalignment='center',
        verticalalignment='top', transform=trans, fontsize = 14 )

    if tight:
        plt.tight_layout()

    if savefile:
        plt.savefig( savefile )

    plt.show()

def subplot_psi_by_alt(   var_df,
                          y_col,
                          pos_col,
                          colors,
                          color_col = 'alt',
                          edge_col = None,
                          ax = None,
                          shade_by_base = False,
                          shade_exons = False,
                          tick_spacing = 10,
                          bar_labels = False,
                          bar_label_loc = 'middle',
                          darken_bars = None,
                          darken_edges = None,
                          labels_legend = False,
                          labels_legend_title = '',
                          labels_legend_loc = 'best',
                          wt_labels = False,
                          title = '',
                          y_ax_lim = None,
                          y_ax_title='',
                          x_ax_title = '',
                          legend = True,
                          legend_title = '',
                          legend_loc = 'best',
                          legend_labels = None,
                          tight = True,
                          print_ex_count = False,
                          scale_bar = False,
                          scale_bar_loc = 'left',
                          rev_trans = False,
                          hlines = None ):

    tbv = var_df.copy()

    if print_ex_count:
        print('This figure shows %i exonic bases.' % sum( c[1]-c[0] for c in shade_exons ) )

    if ax:
        ax = ax
    else:
        ax = plt.gca()

    col_pivot = tbv.pivot( index = [ 'pos', 'ref', pos_col ], columns = color_col, values = y_col )
    #col_pivot = tbv.pivot( index = [ 'pos', 'ref', pos_col ], columns = 'alt', values = y_col )

    col_pivot.plot.bar( ax = ax,
                        color = colors,
                        edgecolor = edge_col,
                        align = 'center',
                        width = 1,
                      )

    if bar_labels or wt_labels or shade_by_base or darken_bars or darken_edges:

        rects = ax.patches

        heights = [ rect.get_height() for rect in rects ]
        neg_ht = [ h for h in heights if h < 0 ]

        if bar_label_loc == 'above':
            trans = tfrms.blended_transform_factory( ax.transData,
                                                     ax.transAxes )
        elif bar_label_loc == 'middle':
            trans = tfrms.blended_transform_factory( ax.transData,
                                                     ax.transData )

    if darken_bars:

        column, colors = darken_bars

        #true false vector for where to put the marker
        _locs = tbv.pivot( index = [ 'pos' ], columns = 'alt', values = column )

        #bars are ordered with all A's first then all C's... this creates a colors vector to match that pattern
        colors_vec = [ c for c in colors for i in range( _locs.shape[0] ) ]

        locs = []
        for col in _locs.columns:
            locs.extend( _locs[ col ].tolist() )

        for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc and not np.isnan( loc ):

                        color = colors_vec[ jidx ]

                        rect.set_facecolor( color )

    if darken_edges:

        column, color = darken_edges

        #true false vector for where to put the marker
        _locs = tbv.pivot( index = [ 'pos' ], columns = 'alt', values = column )

        locs = []
        for col in _locs.columns:
            locs.extend( _locs[ col ].tolist() )

        for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc and not np.isnan( loc ):

                        rect.set_edgecolor( color )

    if wt_labels:

        labels = col_pivot.index.get_level_values( 'ref' )

        min_ht = min( neg_ht ) if len( neg_ht ) > 0 else 0
        neg_ht_mean = np.mean( neg_ht ) if len( neg_ht ) > 0 else -4

        #the way I'm messing with rects should get the label in the middle of the four bars for the pos
        for rect, label in zip( rects[ 2*len( labels ): ], labels ):
            ax.text( rect.get_x(),
                     min_ht + 4*neg_ht_mean,
                     label,
                     fontsize = 18,
                     fontweight = 'bold',
                     ha='center',
                     va='bottom' )

    if bar_labels:

        if bar_label_loc == 'middle':

            max_rect_ht = y_ax_lim[ 1 ] if y_ax_lim else ax.get_ylim()[ 1 ]

            mid_rect_ht = max_rect_ht / 2

            mid_label = int( len( bar_labels ) / 2 )

            #sets up constant y label positions
            #first label is the lowest and they are centered in the middle of the positive y axis
            label_pos = [ mid_rect_ht + ( i - ( mid_label ) )*.2*max_rect_ht for i in range( len( bar_labels ) ) ]

        elif bar_label_loc == 'above':

            label_pos = [ 1 + .1*i for i in range( len( bar_labels ) ) ]

        for iidx, colmark in enumerate( bar_labels ):

            col, marker, _label, fsize = colmark

            #true false vector for where to put the marker
            _locs = tbv.pivot( index = [ 'pos' ], columns = color_col, values = col )

            locs = []
            for col in _locs.columns:
                locs.extend( _locs[ col ].tolist() )

            for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc and not np.isnan( loc ):

                    ax.text( rect.get_x() + rect.get_width() / 2,
                             label_pos[ iidx ],
                             marker,
                             fontsize = fsize,
                             transform = trans,
                             ha = 'center',
                             va = 'bottom' )

    plt.title( title, fontsize = 24 )

    ax.set_ylabel( y_ax_title, fontsize = 18 )
    ax.tick_params( axis = 'y', which='major', labelsize = 24 )

    plt.xlabel( x_ax_title, fontsize = 22 )
    plt.xticks( [ idx for idx,p in enumerate( col_pivot.index.get_level_values( 'pos' ) ) if idx%( tick_spacing ) == 0 ],
                [ c for idx,c in enumerate( col_pivot.index.get_level_values( pos_col ) ) if idx%( tick_spacing ) == 0 ],
                fontsize=18,
                rotation='vertical' )

    if hlines:
        for line in hlines:
            plt.axhline( line, c = 'black', ls = '--', alpha = .6 )

    if shade_by_base:

        for i in range( len( set( tbv.pos ) ) ):

            #shade only every other position
            if i % 2 == 0:

                rect = rects[ i ]

                #this plot has four bases at every position
                ax.axvspan( rect.get_x(),
                             rect.get_x() + 4*rect.get_width(),
                             facecolor = 'gray',
                             alpha = 0.15 )

    if shade_exons:

        for ex in shade_exons:

            plt.axvspan( col_pivot.index.get_loc( ex[0] ).start - .5,
                         col_pivot.index.get_loc( ex[1] ).stop - .5,
                         facecolor = 'gray',
                         alpha = 0.15 )

    if scale_bar:

        assert scale_bar_loc.lower() in [ 'left', 'middle', 'right' ], \
        'Scale bar location can be "left", "middle", or "right".'

        if scale_bar_loc == 'left':
            x_loc = 3*tick_spacing
        elif scale_bar_loc == 'middle':
            x_loc = tbv.shape[0] / 2
        else:
            x_loc = tbv.shape[0] - 3*tick_spacing

        trans = tfrms.blended_transform_factory( ax.transData, ax.transAxes )
        ax.errorbar( x_loc,
                      0.92,
                      xerr = ( tick_spacing / 2 ),
                      color = 'black',
                      capsize=3,
                      transform=trans )

        txt = str( tick_spacing ) + ' bases' if tick_spacing > 1 else str( tick_spacing ) + ' base'
        plt.text( x_loc,
                  0.90,
                  txt,
                  horizontalalignment = 'center',
                  verticalalignment = 'top',
                  transform = trans,
                  fontsize = 14 )

    if legend:

        if legend_labels:
            legend = ax.legend( title = legend_title,
                                 ncol = 2,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 labels = legend_labels,
                                 fontsize = 14 )
        else:
            legend = ax.legend( title = legend_title,
                                 ncol = 2,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )
    elif labels_legend:

        fake_handles = [ plt.Rectangle( (0, 0), 0, 0, fill=False, edgecolor='none', visible = 'False')
                        for i in range( len( bar_labels ) ) ]

        legend = ax.legend(  handles = fake_handles,
                             title = labels_legend_title,
                             bbox_to_anchor = ( 0, 1 ), #if you turn that on you can put legend outside of plot
                             loc = labels_legend_loc,
                             labels = ( symbol + ' '*6 + label for _col, symbol, label, _fsize in bar_labels ),
                             fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )
    else:
        ax.legend_ = None
        plt.draw()

    if tight:
        plt.tight_layout()

    return ax

def subplots_wrapper( var_df,
                      y_cols,
                      pos_col,
                     colors,
                     fig_size = ( 30, 10 ),
                     share_x = True,
                     share_y = True,
                     color_col = 'alt',
                     edge_col = None,
                     shade_by_base = False,
                     shade_exons = False,
                     zoom=None,
                     tick_spacing = 10,
                     wt_labels = False,
                     bar_labels = False,
                     bar_label_loc = 'middle',
                     darken_bars = None,
                     darken_edges = None,
                     labels_legend = False,
                     labels_legend_title = '',
                     labels_legend_loc = 'best',
                     title = '',
                     y_ax_title='',
                     y_ax_lim = None,
                     x_ax_title = '',
                     legend = True,
                     legend_title = '',
                     legend_loc = 'best',
                     legend_labels = None,
                     tight = True,
                     print_ex_count = False,
                     scale_bar = False,
                     scale_bar_loc = 'left',
                     rev_trans = False,
                     hlines = None,
                     savefile = None,
                     ):

    tbv = var_df.sort_values( by = [ 'pos', 'alt' ] ).copy()

    if zoom:
        assert zoom[1] > zoom[0], 'Your final zoom coordinate must be larger than the first zoom coordinate'
        tbv_filt = tbv.set_index( 'pos' ).loc[ zoom[0]:zoom[1] ].reset_index()
    else:
        tbv_filt = tbv

    fig, axes = plt.subplots( len( y_cols),
                              sharex = share_x,
                              sharey = share_y,
                              figsize = fig_size )

    if len( y_cols ) == 1:
        axes = [ axes ]

    if bar_label_loc == 'above':
        all_bar_labels = False
    elif bar_label_loc == 'middle':
        all_bar_labels = bar_labels

    #this whole loop is basically to give each its own y axis labels
    #and to only have the legend & scale bar on the top plot
    #and to only have wt labels on the bottom plot
    for idx, col in enumerate( y_cols ):

        if idx == 0 and len( y_cols ) > 1:
            subplot_psi_by_alt(   tbv_filt,
                                  col,
                                  pos_col,
                                  colors,
                                  ax = axes[ idx ],
                                  color_col = color_col,
                                 edge_col = edge_col,
                                 shade_by_base = shade_by_base,
                                 shade_exons = shade_exons,
                                 tick_spacing = tick_spacing,
                                 bar_labels = bar_labels,
                                 bar_label_loc = bar_label_loc,
                                 darken_bars = darken_bars,
                                 darken_edges = darken_edges,
                                 labels_legend = labels_legend,
                                 labels_legend_title = labels_legend_title,
                                 labels_legend_loc = labels_legend_loc,
                                 title = title,
                                 y_ax_lim = y_ax_lim,
                                 y_ax_title = y_ax_title[ idx ],
                                 legend = legend,
                                 legend_title = legend_title,
                                 legend_loc = legend_loc,
                                 legend_labels = legend_labels,
                                 tight = tight,
                                 print_ex_count = print_ex_count,
                                 scale_bar = scale_bar,
                                 scale_bar_loc = scale_bar_loc,
                                 rev_trans = rev_trans,
                                 hlines = hlines, )
        elif idx == ( len( y_cols ) - 1 ):
            subplot_psi_by_alt(   tbv_filt,
                                  col,
                                  pos_col,
                                  colors,
                                  ax = axes[ idx ],
                                  color_col = color_col,
                                 edge_col = edge_col,
                                 shade_by_base = shade_by_base,
                                 shade_exons = shade_exons,
                                 tick_spacing = tick_spacing,
                                 bar_labels = all_bar_labels,
                                 darken_bars = darken_bars,
                                 darken_edges = darken_edges,
                                 wt_labels = wt_labels,
                                 legend = False,
                                 y_ax_lim = y_ax_lim,
                                 y_ax_title = y_ax_title[ idx ],
                                 tight = tight,
                                 rev_trans = rev_trans,
                                 hlines = hlines,
                                 x_ax_title = x_ax_title, )
        else:
            subplot_psi_by_alt(   tbv_filt,
                                  col,
                                  pos_col,
                                  colors,
                                  ax = axes[ idx ],
                                  color_col = color_col,
                                 edge_col = edge_col,
                                 shade_by_base = shade_by_base,
                                 shade_exons = shade_exons,
                                 tick_spacing = tick_spacing,
                                 bar_labels = all_bar_labels,
                                 darken_bars = darken_bars,
                                 darken_edges = darken_edges,
                                 y_ax_lim = y_ax_lim,
                                 legend = False,
                                 y_ax_title = y_ax_title[ idx ],
                                 tight = tight,
                                 rev_trans = rev_trans,
                                 hlines = hlines, )

    if wt_labels:

        ht_max = tbv_filt[ y_cols ].max().max()
        ht_min = tbv_filt[ y_cols ].min().min()

        y_max = y_ax_lim[ 1 ] if y_ax_lim else 1.01*ht_max
        y_min = ht_min + ht_min if ht_min < 0 else -20
        y_ax_lim = ( y_min,
                     y_max )

    if y_ax_lim:
        plt.ylim( y_ax_lim )

    if savefile:
        plt.savefig( savefile,
                     dpi = 300,
                     bbox_inches = 'tight' )

    plt.show()

def pr_curves( var_df,
               truth_col,
               pred_cols,
               colors,
               fig_size = ( 5, 5 ),
               grid = False,
               x_ax_label = 'Recall\n(%)',
               y_ax_label = 'Precision\n(%)',
               add_point = False,
               savefile = None,
               **kwargs
             ):

    #pr curve function hates missing values
    tbv = var_df.dropna(subset = [ truth_col ] ).copy()

    if tbv.shape[ 0 ] != var_df.shape[ 0 ]:
        print( 'Missing values in truth column.', str( var_df.shape[ 0 ] - tbv.shape[ 0 ] ), 'rows removed.' )

    plt.figure( figsize = fig_size )

    for column, color in zip( pred_cols, colors ):

        col_df = tbv.dropna( subset = [ column ] ).copy()

        if tbv.shape[ 0 ] != col_df.shape[ 0 ]:
            print( 'Missing values in', column, 'column.', \
                   str( tbv.shape[ 0 ] - col_df.shape[ 0 ] ), 'rows removed.' )

        precision, recall, _ = precision_recall_curve( col_df[ truth_col ],
                                                       col_df[ column ] )

        plt.plot( 100*recall,
                  100*precision,
                  color = color,
                  **kwargs )

        print( column, auc( recall, precision ) )

    if add_point:
        plt.plot( add_point[ 0 ],
                  add_point[ 1 ],
                  color = 'black',
                  marker = add_point[ 2 ],
                  markersize = add_point[ 3 ]
                 )

    plt.xlabel( x_ax_label, fontsize = 24 )
    plt.xticks( fontsize = 20 )

    plt.ylabel( y_ax_label, fontsize = 24 )
    plt.yticks( fontsize = 20 )

    plt.grid( grid )

    if savefile:
        plt.savefig( savefile,
                     dpi = 300,
                     bbox_inches = 'tight' )

    plt.show()

def plot_RBP_chg( var_df,
              y_col,
              x_col,
              arrow_col,
              colors,
              color_col = 'alt',
              arrow_filt = None,
              ax = None,
              title = '',
              y_ax_lim = None,
              x_ax_lim = None,
              y_ax_title='',
              x_ax_title = '',
              marker_col = None,
              marker = 'o',
              **kwargs ):

    tbv = var_df.copy()

    if ax:
        ax = ax
    else:
        ax = plt.gca()

    if not marker_col:

        ax.scatter( tbv[ x_col ],
                     tbv[ y_col ],
                     color = tbv[ color_col ],
                     **kwargs
                   )

    else:

        for m, v in marker.items():

            tbv_filt = tbv.loc[ tbv[ marker_col ] == v[ 0 ] ]

            ax.scatter( tbv_filt[ x_col ],
                     tbv_filt[ y_col ],
                     color = tbv_filt[ color_col ],
                    marker = m,
                    s = v[ 1 ],
                     **kwargs
                   )

    arrow_x_start = tbv[ x_col ].values
    arrow_y_start = tbv[ y_col ].values
    arrow_len = tbv[ arrow_col ].values

    #add arrows to plot
    for i, alt in enumerate( tbv.alt ):

        if arrow_filt and arrow_x_start[ i ] < arrow_filt:
            continue

        ax.arrow( arrow_x_start[ i ],
                  arrow_y_start[ i ],
                  0,       #change in x
                  arrow_len[ i ],                      #change in y
                  head_width = 0.3*( np.abs( arrow_len[ i ] ) / 2 ),         #arrow head width
                  head_length = 0.1*( np.abs( arrow_len[ i ] ) / 2 ),        #arrow head length
                  width = 0.03,              #arrow stem width
                  fc = colors[ alt ],             #arrow fill color
                  ec = colors[ alt ]
                )             #arrow edge color

    plt.title( title, fontsize = 24 )

    ax.set_xlabel( x_ax_title, fontsize = 18 )
    ax.tick_params( axis = 'x', which='major', labelsize = 12 )

    ax.set_ylabel( y_ax_title, fontsize = 18 )
    ax.tick_params( axis = 'y', which='major', labelsize = 12 )

    if x_ax_lim:
        ax.set_xlim( x_ax_lim )

    if y_ax_lim:
        ax.set_ylim( y_ax_lim )

    plt.show()

def subplot_psi_wt(   var_df,
                          y_col,
                          pos_col,
                          colors,
                          color_col = 'ref',
                          edge_col = None,
                          ax = None,
                           zoom = None,
                          shade_exons = False,
                          tick_spacing = 10,
                          title = '',
                          y_ax_lim = None,
                          y_ax_title='',
                          x_ax_title = '',
                          legend = True,
                          legend_title = '',
                          legend_loc = 'best',
                          legend_labels = None,
                          tight = True,
                          print_ex_count = False,
                          scale_bar = False,
                          scale_bar_loc = 'left',
                          rev_trans = False,
                          hlines = None,
                          snap = False):

    tbv = var_df.copy()

    if print_ex_count:
        print('This figure shows %i exonic bases.' % sum( c[1]-c[0] for c in shade_exons ) )

    if ax:
        ax = ax
    else:
        ax = plt.gca()

    d2c = dict( zip( [ 'A', 'C', 'G', 'T' ], colors) )

    if zoom:
        assert zoom[1] > zoom[0], 'Your final zoom coordinate must be larger than the first zoom coordinate'
        tbv = tbv.set_index( 'pos' ).loc[ zoom[0]:zoom[1] ].reset_index()

    tbv.plot.bar( x = 'pos',
                  y = y_col,
                        ax = ax,
                        color = map( d2c.get, tbv[ color_col ] ),
                        edgecolor = edge_col,
                        align = 'center',
                        width = 1,
                      )

    plt.title( title, fontsize = 24 )

    ax.set_ylabel( y_ax_title, fontsize = 18 )
    ax.tick_params( axis = 'y', which='major', labelsize = 36 )

    plt.xlabel( x_ax_title, fontsize = 36 )
    plt.xticks( [ idx for idx,p in enumerate( tbv.pos ) if idx%( tick_spacing ) == 0 ],
                [ c for idx,c in enumerate( tbv[ pos_col ] ) if idx%( tick_spacing ) == 0 ],
                fontsize=36,
                rotation='vertical' )

    if hlines:
        for line in hlines:
            x_cds, col, style, lw = line
            ax.axhline( x_cds, c = col, ls = style, linewidth = lw, alpha = .6 )

    if shade_exons:

        for ex in shade_exons:

            plt.axvspan( col_pivot.index.get_loc( ex[0] ).start - .5,
                         col_pivot.index.get_loc( ex[1] ).stop,
                         facecolor = 'gray',
                         alpha = 0.15 )

    if scale_bar:

        assert scale_bar_loc.lower() in [ 'left', 'middle', 'right' ], \
        'Scale bar location can be "left", "middle", or "right".'

        if scale_bar_loc == 'left':
            x_loc = 3*tick_spacing
        elif scale_bar_loc == 'middle':
            x_loc = tbv.shape[0] / 2
        else:
            x_loc = tbv.shape[0] - 3*tick_spacing

        trans = tfrms.blended_transform_factory( ax.transData, ax.transAxes )
        ax.errorbar( x_loc,
                      0.92,
                      xerr = ( tick_spacing / 2 ),
                      color = 'black',
                      capsize=3,
                      transform=trans )

        txt = str( tick_spacing ) + ' bases' if tick_spacing > 1 else str( tick_spacing ) + ' base'
        plt.text( x_loc,
                  0.90,
                  txt,
                  horizontalalignment = 'center',
                  verticalalignment = 'top',
                  transform = trans,
                  fontsize = 14 )

    if legend:

        legend = ax.legend( title = legend_title,
                                 ncol = 2,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )

        plt.setp( legend.get_title(), fontsize=14 )
    else:
        ax.legend_ = None
        plt.draw()

    if y_ax_lim:
        ax.set_ylim( y_ax_lim )

    art.set_snap( ax, snap )

    if tight:
        plt.tight_layout()

    return ax

def subplot_psi_by_alt(   var_df,
                          y_col,
                          pos_col,
                          colors,
                          color_col = 'alt',
                          edge_col = None,
                          ax = None,
                          shade_by_base = False,
                          shade_exons = False,
                          tick_spacing = 10,
                          bar_labels = False,
                          bar_label_loc = 'middle',
                          bar_labels_offset = True,
                          darken_bars = None,
                          darken_bars2 = None,
                          darken_edges = None,
                          labels_legend = False,
                          labels_legend_title = '',
                          labels_legend_loc = 'best',
                          wt_labels = False,
                          title = '',
                          y_ax_lim = None,
                          y_ax_title='',
                          x_ax_title = '',
                          legend = True,
                          legend_title = '',
                          legend_loc = 'best',
                          legend_labels = None,
                          tight = True,
                          print_ex_count = False,
                          scale_bar = False,
                          scale_bar_loc = 'left',
                          rev_trans = False,
                          hlines = None,
                          snap = False ):

    tbv = var_df.copy()

    if print_ex_count:
        print('This figure shows %i exonic bases.' % sum( c[1]-c[0] for c in shade_exons ) )

    if ax:
        ax = ax
    else:
        ax = plt.gca()


    col_pivot = tbv.pivot( index = [ 'pos', 'ref', pos_col ], columns = 'alt', values = y_col )

    col_pivot.plot.bar( ax = ax,
                        color = colors,
                        edgecolor = edge_col,
                        align = 'center',
                        width = 1,
                      )

    if bar_labels or wt_labels or shade_by_base or darken_bars or darken_edges:

        rects = ax.patches

        heights = [ rect.get_height() for rect in rects ]
        neg_ht = [ h for h in heights if h < 0 ]

        if bar_label_loc == 'above':
            trans = tfrms.blended_transform_factory( ax.transData,
                                                     ax.transAxes )
        elif bar_label_loc == 'middle':
            trans = tfrms.blended_transform_factory( ax.transData,
                                                     ax.transData )

    if darken_bars2:

        column, colors = darken_bars2

        #true false vector for where to put the marker
        _locs = tbv.pivot( index = [ 'pos' ], columns = 'alt', values = column )

        #bars are ordered with all A's first then all C's... this creates a colors vector to match that pattern
        colors_vec = [ c for c in colors for i in range( _locs.shape[0] ) ]

        locs = []
        for col in _locs.columns:
            locs.extend( _locs[ col ].tolist() )

        for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc and not np.isnan( loc ):

                        color = colors_vec[ jidx ]

                        rect.set_facecolor( color )

    if darken_bars:

        column, colors = darken_bars

        #true false vector for where to put the marker
        _locs = tbv.pivot( index = [ 'pos' ], columns = 'alt', values = column )

        #bars are ordered with all A's first then all C's... this creates a colors vector to match that pattern
        colors_vec = [ c for c in colors for i in range( _locs.shape[0] ) ]

        locs = []
        for col in _locs.columns:
            locs.extend( _locs[ col ].tolist() )

        for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc and not np.isnan( loc ):

                    color = colors_vec[ jidx ]

                    rect.set_facecolor( color )

    if darken_edges:

        column, color = darken_edges

        #true false vector for where to put the marker
        _locs = tbv.pivot( index = [ 'pos' ], columns = 'alt', values = column )

        locs = []
        for col in _locs.columns:
            locs.extend( _locs[ col ].tolist() )

        for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc and not np.isnan( loc ):

                        rect.set_edgecolor( color )

    if wt_labels:

        labels = col_pivot.index.get_level_values( 'ref' )

        min_ht = min( neg_ht ) if len( neg_ht ) > 0 else 0
        neg_ht_mean = np.mean( neg_ht ) if len( neg_ht ) > 0 else -4

        #the way I'm messing with rects should get the label in the middle of the four bars for the pos
        for rect, label in zip( rects[ 2*len( labels ): ], labels ):
            ax.text( rect.get_x(),
                     min_ht + 4*neg_ht_mean,
                     label,
                     fontsize = 18,
                     fontweight = 'bold',
                     ha='center',
                     va='bottom' )

    if bar_labels:

        if bar_labels_offset:

            if bar_label_loc == 'middle':

                max_rect_ht = y_ax_lim[ 1 ] if y_ax_lim else ax.get_ylim()[ 1 ]

                mid_rect_ht = max_rect_ht / 2

                mid_label = int( len( bar_labels ) / 2 )

                #sets up constant y label positions
                #first label is the lowest and they are centered in the middle of the positive y axis
                label_pos = [ mid_rect_ht + ( i - ( mid_label ) )*.2*max_rect_ht for i in range( len( bar_labels ) ) ]

            elif bar_label_loc == 'above':

                label_pos = [ 1 + .1*i for i in range( len( bar_labels ) ) ]

        else:

            if bar_label_loc == 'middle':

                pos_ht_mean = np.mean( pos_ht_sum )

                max_rect_ht = y_ax_lim[ 1 ] if y_ax_lim else ax.get_ylim()[ 1 ]

                mid_rect_ht = max_rect_ht / 2

                #sets up constant y label positions
                #they are centered in the middle of the positive y axis
                label_pos = [ mid_rect_ht for i in range( len( bar_labels ) ) ]

            elif bar_label_loc == 'above':

                label_pos = [ .885 if marker == r"$\bullet$" else 1.01 for _col, marker, _label, _fsize in bar_labels ]

        for iidx, colmark in enumerate( bar_labels ):

            col, marker, _label, fsize = colmark

            #true false vector for where to put the marker
            _locs = tbv.pivot( index = [ 'pos' ], columns = 'alt', values = col )

            locs = []
            for col in _locs.columns:
                locs.extend( _locs[ col ].tolist() )

            for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                xpos = rect.get_x() + rect.get_width() / 2

                if marker == r"$\bullet$":
                    xpos = .99*( xpos )

                if loc and not np.isnan( loc ):

                    ax.text( xpos,
                             label_pos[ iidx ],
                             marker,
                             fontsize = fsize,
                             transform = trans,
                             ha = 'center',
                             va = 'bottom' )

    plt.title( title, fontsize = 24 )

    ax.set_ylabel( y_ax_title, fontsize = 18 )
    ax.tick_params( axis = 'y', which='major', labelsize = 36 )

    plt.xlabel( x_ax_title, fontsize = 36 )
    plt.xticks( [ idx for idx,p in enumerate( col_pivot.index.get_level_values( 'pos' ) ) if idx%( tick_spacing ) == 0 ],
                [ c for idx,c in enumerate( col_pivot.index.get_level_values( pos_col ) ) if idx%( tick_spacing ) == 0 ],
                fontsize=36,
                rotation='vertical' )

    if hlines:
        for line in hlines:
            x_cds, col, style, lw = line
            ax.axhline( x_cds, c = col, ls = style, linewidth = lw, alpha = .6 )

    if shade_by_base:

        for i in range( len( set( tbv.pos ) ) ):

            #shade only every other position
            if i % 2 == 0:

                rect = rects[ i ]

                #this plot has four bases at every position
                ax.axvspan( rect.get_x(),
                             rect.get_x() + 4*rect.get_width(),
                             facecolor = 'gray',
                             alpha = 0.15 )

    if shade_exons:

        for ex in shade_exons:

            plt.axvspan( col_pivot.index.get_loc( ex[0] ).start - .5,
                         col_pivot.index.get_loc( ex[1] ).stop,
                         facecolor = 'gray',
                         alpha = 0.15 )

    if scale_bar:

        assert scale_bar_loc.lower() in [ 'left', 'middle', 'right' ], \
        'Scale bar location can be "left", "middle", or "right".'

        if scale_bar_loc == 'left':
            x_loc = 3*tick_spacing
        elif scale_bar_loc == 'middle':
            x_loc = tbv.shape[0] / 2
        else:
            x_loc = tbv.shape[0] - 3*tick_spacing

        trans = tfrms.blended_transform_factory( ax.transData, ax.transAxes )
        ax.errorbar( x_loc,
                      0.92,
                      xerr = ( tick_spacing / 2 ),
                      color = 'black',
                      capsize=3,
                      transform=trans )

        txt = str( tick_spacing ) + ' bases' if tick_spacing > 1 else str( tick_spacing ) + ' base'
        plt.text( x_loc,
                  0.90,
                  txt,
                  horizontalalignment = 'center',
                  verticalalignment = 'top',
                  transform = trans,
                  fontsize = 14 )

    if legend:

        if legend_labels:
            legend = ax.legend( title = legend_title,
                                 ncol = 2,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 labels = legend_labels,
                                 fontsize = 14 )
        else:
            legend = ax.legend( title = legend_title,
                                 ncol = 2,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )
    elif labels_legend:

        fake_handles = [ plt.Rectangle( (0, 0), 0, 0, fill=False, edgecolor='none', visible = 'False')
                        for i in range( len( bar_labels ) ) ]

        legend = ax.legend(  handles = fake_handles,
                             title = labels_legend_title,
                             bbox_to_anchor = ( 0, 1 ), #if you turn that on you can put legend outside of plot
                             loc = labels_legend_loc,
                             labels = ( symbol + ' '*6 + label for _col, symbol, label, _fsize in bar_labels ),
                             fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )
    else:
        ax.legend_ = None
        plt.draw()

    if y_ax_lim:
        ax.set_ylim( y_ax_lim )

    art.set_snap( ax, snap )

    if tight:
        plt.tight_layout()

    return ax

def subplots_wrapper( var_df,
                      y_cols,
                      pos_col,
                     colors,
                     fig_size = ( 30, 10 ),
                     share_x = True,
                     share_y = True,
                     color_col = 'alt',
                     edge_col = None,
                     shade_by_base = False,
                     shade_exons = False,
                     zoom=None,
                     tick_spacing = 10,
                     wt_labels = False,
                     bar_labels = False,
                     bar_label_loc = 'middle',
                     bar_labels_offset = True,
                     darken_bars = None,
                     darken_bars2 = None,
                     darken_edges = None,
                     labels_legend = False,
                     labels_legend_title = '',
                     labels_legend_loc = 'best',
                     title = '',
                     y_ax_title='',
                     y_ax_lim = None,
                     x_ax_title = '',
                     legend = True,
                     legend_title = '',
                     legend_loc = 'best',
                     legend_labels = None,
                     tight = True,
                     print_ex_count = False,
                     scale_bar = False,
                     scale_bar_loc = 'left',
                     rev_trans = False,
                     hlines = None,
                     savefile = None,
                     snap = False
                     ):

    tbv = var_df.sort_values( by = [ 'pos', 'alt' ] ).copy()

    if zoom:
        assert zoom[1] > zoom[0], 'Your final zoom coordinate must be larger than the first zoom coordinate'
        tbv_filt = tbv.set_index( 'pos' ).loc[ zoom[0]:zoom[1] ].reset_index()
    else:
        tbv_filt = tbv

    fig, axes = plt.subplots( len( y_cols),
                              sharex = share_x,
                              sharey = share_y,
                              figsize = fig_size )

    if len( y_cols ) == 1:
        axes = [ axes ]

    if bar_label_loc == 'above':
        all_bar_labels = False
    elif bar_label_loc == 'middle':
        all_bar_labels = bar_labels

    if not darken_bars:
        dbar_dict = { i: None for i in range( len( y_cols) ) }
    elif len( darken_bars ) == 1:
        dbar_dict = { i: darken_bars[ 0 ] for i in range( len( y_cols) ) }
    elif len( darken_bars ) == len( y_cols ):
        dbar_dict = { i: darken_bars[ i ] for i in range( len( y_cols) ) }
    else:
        print( 'Something is wrong with how you entered the darken bars argument!' )

    if not y_ax_lim:
        ylim_dict = { i: None for i in range( len( y_cols) ) }
    elif len( y_ax_lim ) == 1:
        ylim_dict = { i: y_ax_lim[ 0 ] for i in range( len( y_cols) ) }
    elif len( y_ax_lim ) > 1 and share_y:
        print( 'You specified too many y axis limits to have them all sharey' )
    elif len( y_ax_lim ) == len( y_cols ):
        ylim_dict = { i: y_ax_lim[ i ] for i in range( len( y_cols) ) }
    else:
        print( 'Something is wrong with how you entered the y axis limits argument!' )

    if wt_labels:

        ht_max = tbv_filt[ y_cols ].max().max()
        ht_min = tbv_filt[ y_cols ].min().min()

        y_max = y_ax_lim[ -1 ][ 1 ] if y_ax_lim else 1.01*ht_max
        y_min = ht_min + ht_min if ht_min < 0 else -20
        ylim_dict[ len( y_cols ) - 1 ] = ( y_min, y_max )

    #this whole loop is basically to give each its own y axis labels
    #and to only have the legend & scale bar on the top plot
    #and to only have wt labels on the bottom plot
    for idx, col in enumerate( y_cols ):

        if idx == 0:
            subplot_psi_by_alt(   tbv_filt,
                                  col,
                                  pos_col,
                                  colors,
                                  ax = axes[ idx ],
                                  color_col = color_col,
                                 edge_col = edge_col,
                                 shade_by_base = shade_by_base,
                                 shade_exons = shade_exons,
                                 tick_spacing = tick_spacing,
                                 bar_labels = bar_labels,
                                 bar_label_loc = bar_label_loc,
                                 bar_labels_offset = bar_labels_offset,
                                 darken_bars = dbar_dict[ idx ],
                                 darken_bars2 = darken_bars2,
                                 darken_edges = darken_edges,
                                 labels_legend = labels_legend,
                                 labels_legend_title = labels_legend_title,
                                 labels_legend_loc = labels_legend_loc,
                                 title = title,
                                 y_ax_lim = ylim_dict[ idx ],
                                 y_ax_title = y_ax_title[ idx ],
                                 legend = legend,
                                 legend_title = legend_title,
                                 legend_loc = legend_loc,
                                 legend_labels = legend_labels,
                                 tight = tight,
                                 print_ex_count = print_ex_count,
                                 scale_bar = scale_bar,
                                 scale_bar_loc = scale_bar_loc,
                                 rev_trans = rev_trans,
                                 hlines = hlines,
                                 snap = snap )
        elif idx == ( len( y_cols ) - 1 ):
            subplot_psi_by_alt(   tbv_filt,
                                  col,
                                  pos_col,
                                  colors,
                                  ax = axes[ idx ],
                                  color_col = color_col,
                                 edge_col = edge_col,
                                 shade_by_base = shade_by_base,
                                 shade_exons = shade_exons,
                                 tick_spacing = tick_spacing,
                                 bar_labels = all_bar_labels,
                                 darken_bars = dbar_dict[ idx ],
                                 darken_bars2 = darken_bars2,
                                 darken_edges = darken_edges,
                                 wt_labels = wt_labels,
                                 legend = False,
                                 y_ax_lim = ylim_dict[ idx ],
                                 y_ax_title = y_ax_title[ idx ],
                                 tight = tight,
                                 rev_trans = rev_trans,
                                 hlines = hlines,
                                 x_ax_title = x_ax_title,
                                 snap = snap )
        else:
            subplot_psi_by_alt(   tbv_filt,
                                  col,
                                  pos_col,
                                  colors,
                                  ax = axes[ idx ],
                                  color_col = color_col,
                                 edge_col = edge_col,
                                 shade_by_base = shade_by_base,
                                 shade_exons = shade_exons,
                                 tick_spacing = tick_spacing,
                                 bar_labels = all_bar_labels,
                                 darken_bars = dbar_dict[ idx ],
                                 darken_bars2 = darken_bars2,
                                 darken_edges = darken_edges,
                                 #y_ax_lim = ylim_dict[ idx ],
                                 legend = False,
                                 y_ax_title = y_ax_title[ idx ],
                                 tight = tight,
                                 rev_trans = rev_trans,
                                 hlines = hlines,
                                 snap = snap )

    #if y_ax_lim:
        #plt.ylim( y_ax_lim )

    if savefile:
        plt.savefig( savefile,
                     #dpi = 300,
                     bbox_inches = 'tight'
                   )

    plt.show()

def split_ax_bcs( var_df,
                  bc_col,
                  pos_col,
                  colors,
                  y_lims,
                  index_cols = [ 'pos', 'alt' ],
                  fig_size = ( 30, 5 ),
                  hratios = [ 1, 1 ],
                  zoom=None,
                  tick_spacing = 10,
                  title = '',
                  y_ax_title='',
                  y_ax_lim = None,
                  x_ax_title = '',
                  legend = True,
                  legend_title = '',
                  legend_loc = 'best',
                  legend_labels = None,
                  tight = True,
                  savefile = None
                  ):

    tbv = var_df.sort_values( by = index_cols ).copy()

    if zoom:
        assert zoom[1] > zoom[0], 'Your final zoom coordinate must be larger than the first zoom coordinate'
        tbv_filt = tbv.set_index( 'pos' ).loc[ zoom[0]:zoom[1] ].reset_index()
    else:
        tbv_filt = tbv

    #check if plot will lose bars due to not enough space for all the pixels
    dpi = 100
    if ( tbv_filt.shape[0]*1.5 ) > ( fig_size[0]*dpi ):
        fig_height = fig_size[1]
        fig_size = ( ( tbv_filt.shape[0]*1.5 / dpi ), fig_height )
        print('Adjusting figure width to accomodate all pixels...')

    fig = plt.figure( figsize = fig_size )
    gs = gridspec.GridSpec( nrows = len( y_lims ),
                            ncols = 1,
                            height_ratios = hratios )

    col_pivot = tbv_filt.pivot( index = 'pos', columns = 'alt', values = bc_col )

    axes = []
    for i in range( len( y_lims ) ):
        axes.append( plt.subplot( gs[ i ],
                                  #sharex = True
                                )
                    )
        col_pivot.plot.bar( ax = axes[ -1 ],
                            color = colors,
                            stacked = True,
                            align = 'center',
                            width = 1, )
        axes[ - 1 ].set_ylim( y_lims[ i ] )
        axes[ -1 ].tick_params( axis = 'y', which='major', labelsize = 36 )

    for i, ax in enumerate( axes ):
        if i != ( len( y_lims ) - 1 ):
            ax.spines[ 'bottom' ].set_visible( False )
            ax.tick_params( axis='x',          # changes apply to the x-axis
                            which='both',      # both major and minor ticks are affected
                            bottom=False,      # ticks along the bottom edge are off
                            top=False,         # ticks along the top edge are off
                            labelbottom=False) # labels along the bottom edge are off
            ax.set_xlabel( '' )
        else:
            ax.set_xlabel( x_ax_title, fontsize = 36 )
            ax.tick_params( axis = 'x', which='major', labelsize = 36, labelrotation = 90 )

        if i != 0:
            ax.spines[ 'top' ].set_visible( False )

    pos_pivot = tbv_filt.pivot( index = 'pos', columns = 'alt', values = pos_col )
    pos_pivot = pos_pivot.fillna( '' )
    pos_pivot[ pos_col ] = [ a if a != '' else c for a,c in zip( pos_pivot.A, pos_pivot.C ) ]
    plt.xticks( [ idx for idx,p in enumerate( col_pivot.index ) if idx % tick_spacing == 1 ],
                [ c for idx,c in enumerate( pos_pivot[ pos_col ] ) if idx % tick_spacing == 1 ],)

    if legend:

        if legend_labels:
            legend = plt.legend( title = legend_title,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 labels = legend_labels,
                                 fontsize = 14 )
        else:
            legend = plt.legend( title = legend_title,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )
    else:
        for ax in axes:
            ax.legend_ = None
            plt.draw()

    if tight:
        plt.tight_layout()

    if savefile:
        plt.savefig( savefile,
                     #dpi = 300,
                     bbox_inches = 'tight'
                   )

    plt.show()

def plot_stacked_psi(     var_df,
                          zcols,
                          pos_col,
                          colors,
                          index_cols = [ 'pos', 'alt' ],
                          edge_col = None,
                          fig_size = ( 30, 5 ),
                          shade_exons = False,
                          shade_by_base = False,
                          zoom=None,
                          tick_spacing = 10,
                          alt_labels = False,
                          bar_labels = False,
                          bar_label_loc = 'middle',
                          bar_labels_offset = True,
                          darken_bars = False,
                          darken_edges = False,
                          labels_legend = False,
                          labels_legend_title = '',
                          labels_legend_loc = 'best',
                          title = '',
                          y_ax_title='',
                          y_ax_lim = None,
                          x_ax_title = '',
                          legend = True,
                          legend_title = '',
                          legend_loc = 'best',
                          legend_labels = None,
                          tight = True,
                          print_ex_count = False,
                          scale_bar = False,
                          scale_bar_loc = 'left',
                          rev_trans = False,
                          hlines = None,
                          savefile = None ):

    tbv = var_df.sort_values( by = index_cols ).copy()

    if zoom:
        assert zoom[1] > zoom[0], 'Your final zoom coordinate must be larger than the first zoom coordinate'
        tbv_filt = tbv.set_index( 'pos' ).loc[ zoom[0]:zoom[1] ].reset_index()
    else:
        tbv_filt = tbv

    #check if plot will lose bars due to not enough space for all the pixels
    dpi = 100
    if ( tbv_filt.shape[0]*1.5 ) > ( fig_size[0]*dpi ):
        fig_height = fig_size[1]
        fig_size = ( ( tbv_filt.shape[0]*1.5 / dpi ), fig_height )
        print('Adjusting figure width to accomodate all pixels...')

    if print_ex_count:
        print('This figure shows %i exonic bases.' % sum( c[1]-c[0] for c in shade_exons ) )

    #col_pivot = tbv_filt.set_index( index_cols )[ zcols ]

    col_pivot = tbv_filt.pivot( index = 'pos', columns = 'alt', values = zcols )

    col_pivot.plot.bar( color = colors,
                        edgecolor = edge_col,
                        stacked = True,
                        align = 'center',
                        width = 1,
                        figsize = fig_size )

    ax = plt.gca()

    if alt_labels or bar_labels or shade_by_base or scale_bar or darken_bars or darken_edges:

        rects = ax.patches

        heights = [ rect.get_height() for rect in rects ]

        n_bars = int( len( heights ) / len( zcols ) )

        pos_ht_sum = [ sum( heights[ i ] for i in range( j, len( heights ), n_bars ) if heights[ i ] > 0 )
                       for j in range( n_bars ) ]
        neg_ht_sum = [ sum( heights[ i ] for i in range( j, len( heights ), n_bars ) if heights[ i ] < 0 )
                       for j in range( n_bars ) ]

        if bar_label_loc == 'above':
            trans = tfrms.blended_transform_factory( ax.transData,
                                                     ax.transAxes )
        else:
            trans = tfrms.blended_transform_factory( ax.transData,
                                                     ax.transData )

    if darken_bars:

        column, colors = darken_bars

        locs = tbv_filt[ column ]

        for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc:

                    for iidx in range( len( zcols ) ):

                        color = colors[ iidx ]

                        rects[ jidx + iidx*n_bars ].set_facecolor( color )

    if darken_edges:

        column, color = darken_edges

        locs = tbv_filt[ column ]

        for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc:

                    for iidx in range( len( zcols ) ):

                        rects[ jidx + iidx*n_bars ].set_edgecolor( color )

    if alt_labels:

        labels = col_pivot.index.get_level_values( 'alt' )

        min_ht = min( neg_ht_sum )
        neg_ht_mean = np.mean( neg_ht_sum )

        for rect, label in zip( rects, labels ):
            ax.text( rect.get_x() + rect.get_width() / 2,
                     min_ht + .75*neg_ht_mean,
                     label,
                     fontsize = 16,
                     fontweight = 'bold',
                     ha='center',
                     va='bottom' )

    if bar_labels:

        if bar_labels_offset:

            if bar_label_loc == 'middle':

                pos_ht_mean = np.mean( pos_ht_sum )

                max_rect_ht = y_ax_lim[ 1 ] if y_ax_lim else ax.get_ylim()[ 1 ]

                mid_rect_ht = max_rect_ht / 2

                mid_label = int( len( bar_labels ) / 2 )

                #sets up constant y label positions
                #first label is the lowest and they are centered in the middle of the positive y axis
                label_pos = [ mid_rect_ht + ( i - ( mid_label ) )*.2*max_rect_ht for i in range( len( bar_labels ) ) ]

            elif bar_label_loc == 'above':

                label_pos = [ 1 + .1*i for i in range( len( bar_labels ) ) ]

        else:

            if bar_label_loc == 'middle':

                pos_ht_mean = np.mean( pos_ht_sum )

                max_rect_ht = y_ax_lim[ 1 ] if y_ax_lim else ax.get_ylim()[ 1 ]

                mid_rect_ht = max_rect_ht / 2

                #sets up constant y label positions
                #they are centered in the middle of the positive y axis
                label_pos = [ mid_rect_ht for i in range( len( bar_labels ) ) ]

            elif bar_label_loc == 'above':

                label_pos = [ 1.01 for i in range( len( bar_labels ) ) ]

        for iidx, colmark in enumerate( bar_labels ):

            col, marker, _label, fsize = colmark

            #true false vector for where to put the marker
            locs = tbv_filt[ col ]

            for jidx, rloc in enumerate( zip( rects, locs ) ):

                rect, loc = rloc

                if loc:

                    ax.text( rect.get_x() + rect.get_width() /2,
                             label_pos[ iidx ],
                             marker,
                             fontsize = fsize,
                             transform = trans,
                             ha = 'center',
                             va = 'bottom' )

    plt.title( title, fontsize = 24 )

    if alt_labels:

        col_pivot[ 'y_min' ] = col_pivot[ col_pivot < 0 ].sum( axis = 1 )
        col_pivot[ 'y_max' ] = col_pivot[ col_pivot > 0 ].sum( axis = 1 )

        y_max = y_ax_lim[ 1 ] if y_ax_lim else 1.01*col_pivot.y_max.max()
        y_ax_lim = ( col_pivot.y_min.min() + neg_ht_mean,
                     y_max )

    if y_ax_lim:
        plt.ylim( y_ax_lim )

    plt.ylabel( y_ax_title, fontsize = 36 )
    plt.yticks( fontsize = 36 )

    plt.xlabel( x_ax_title, fontsize = 36 )

    pos_pivot = tbv_filt.pivot( index = 'pos', columns = 'alt', values = pos_col )

    pos_pivot = pos_pivot.fillna( '' )

    pos_pivot[ pos_col ] = [ a if a != '' else c for a,c in zip( pos_pivot.A, pos_pivot.C ) ]

    plt.xticks( [ idx for idx,p in enumerate( col_pivot.index ) if idx % tick_spacing == 1 ],
                [ c for idx,c in enumerate( pos_pivot[ pos_col ] ) if idx % tick_spacing == 1 ],
                fontsize=36,
                rotation='vertical' )

    if hlines:
        for line in hlines:
            plt.axhline( line, c = 'black', ls = '--', alpha = .6 )

    if shade_by_base:

        for i in range( len( set( tbv_filt.pos ) ) ):

            #shade only every other position
            if i % 2 == 0:

                #three bases are non-missing at each position
                rect = rects[ 3*i ]

                plt.axvspan( rect.get_x(),
                             rect.get_x() + 3*rect.get_width(),
                             facecolor = 'gray',
                             alpha = 0.2 )

    if shade_exons:
        for ex in shade_exons:
            plt.axvspan( col_pivot.index.get_loc( ex[0] ).start - .5,
                         col_pivot.index.get_loc( ex[1] ).stop,
                         facecolor = 'gray',
                         alpha = 0.15 )

    if scale_bar:

        assert scale_bar_loc.lower() in [ 'left', 'middle', 'right' ], \
        'Scale bar location can be "left", "middle", or "right".'

        if scale_bar_loc == 'left':
            x_loc = 3*tick_spacing
        elif scale_bar_loc == 'middle':
            x_loc = tbv_filt.shape[0] / 2
        else:
            x_loc = tbv_filt.shape[0] - 3*tick_spacing

        trans = tfrms.blended_transform_factory( ax.transData,
                                                 ax.transAxes )
        plt.errorbar( x_loc,
                      0.96,
                      #xerr = ( 3*tick_spacing / 2 ),
                      xerr = ( 3*tick_spacing*rects[0].get_width() / 2 ),
                      color = 'black',
                      capsize=3,
                      transform=trans )

        txt = str( tick_spacing ) + ' bases' if tick_spacing > 1 else str( tick_spacing ) + ' base'
        plt.text( x_loc,
                  0.94,
                  txt,
                  horizontalalignment = 'center',
                  verticalalignment = 'top',
                  transform = trans,
                  fontsize = 14 )

    if legend:

        if legend_labels:
            legend = plt.legend( title = legend_title,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 labels = legend_labels,
                                 fontsize = 14 )
        else:
            legend = plt.legend( title = legend_title,
                                 bbox_to_anchor = ( 1, 1 ), #if you turn that on you can put legend outside of plot
                                 loc = legend_loc,
                                 fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )
    elif labels_legend:

        fake_handles = [ plt.Rectangle( (0, 0), 0, 0, fill=False, edgecolor='none', visible = 'False')
                        for i in range( len( bar_labels ) ) ]

        legend = plt.legend( handles = fake_handles,
                             title = labels_legend_title,
                             bbox_to_anchor = ( 0, 1 ), #if you turn that on you can put legend outside of plot
                             loc = labels_legend_loc,
                             labels = ( symbol + ' '*6 + label for _col, symbol, label, _fsize in bar_labels ),
                             fontsize = 14 )
        plt.setp( legend.get_title(), fontsize=14 )
    else:
        ax.legend_ = None
        plt.draw()

    if tight:
        plt.tight_layout()

    if savefile:
        plt.savefig( savefile,
                     #dpi = 300,
                     #bbox_inches = 'tight'
                   )

    plt.show()

def swarm_plot( tblbyvar,
                x_col,
                y_col,
                x_ax_label = '',
                x_tick_labels = None,
                y_ax_label = '',
                hlines = None,
                savefile = None,
                **kwargs ):

    sns.set_style( 'ticks' )

    tbv = tblbyvar.copy()

    ax = sns.swarmplot( x = x_col,
                        y = y_col,
                        data = tbv,
                        **kwargs )

    ax.set_xlabel( x_ax_label,
                   fontsize = 18 )

    ax.set_ylabel( y_ax_label,
                   fontsize = 18 )

    if x_tick_labels:
        ax.set_xticklabels( x_tick_labels,
                            fontsize = 14 )

    plt.yticks( fontsize = 14 )

    if hlines:
        for hl in hlines:
            ax.axhline( hl, color = 'black', linestyle = '--' )

    if savefile:
        plt.savefig( savefile )

    plt.show()

def violin( tbl_byvar,
            cat_col,
            val_col,
            ylim = None,
            x_ax_label = '',
            y_ax_label = '',
            x_tick_labels = None,
            savefig = None,
            **kwargs ):

    tbv = tbl_byvar.copy()

    ax = sns.violinplot( x = cat_col,
                        y = val_col,
                        data = tbv,
                        **kwargs
                         )

    if ylim:
        ax.set_ylim( ylim )

    ax.set_xlabel( x_ax_label,
                   fontsize = 18 )

    ax.set_ylabel( y_ax_label,
                   fontsize = 18 )

    if x_tick_labels:
        ax.set_xticklabels( x_tick_labels,
                            fontsize = 14 )

    plt.yticks( fontsize = 14 )

    if savefig:
        plt.savefig( savefig,
                     bbox_to_inches = 'tight' )

    plt.show()

def upset_plot( byvartbl,
                cat_cols,
                fig_size = ( 10, 20 ),
                dist_cols = None,
                dist_cols_kind = None,
                dist_cols_col = None,
                violin_inner = None,
                savefig = None,
                **kwargs ):

    tbv = byvartbl.set_index( cat_cols ).copy()

    up_df = up.UpSet( tbv, **kwargs )

    inner = {  }


    if dist_cols:

        for i, col in enumerate( dist_cols ):

            if dist_cols_kind[ i ] == 'violin':

                up_df.add_catplot( value = col,
                                   kind = dist_cols_kind[ i ],
                                   color = dist_cols_col[ i ],
                                   inner = 'points',
                                   #jitter = True,
                                   #s = 1.5,
                                   elements = 3,
                                   cut = 0,
                                   scale = 'width',
                                   #linewidth = 1,
                                   width = .5,
                                   #bw = 2
                                 )
            else:

                up_df.add_catplot( value = col,
                                   kind = dist_cols_kind[ i ],
                                   color = dist_cols_col[ i ],
                                   #s = 3,
                                   elements = 2,
                                   #linewidth = 0,
                                 )

    #fig1 = plt.figure( figsize = fig_size )

    fig1 = plt.gcf()

    fig1.set_size_inches(20, 40)

    up_df.plot( fig = fig1 )

    #g.fig.set_figwidth(20)
    #g.fig.set_figheight(20)

    if savefig:

        plt.savefig( savefig,
                     bbox_inches = 'tight' )

    plt.show()

def waterfall_plot(
                    vartbl,
                    waterfall_col,
                    percentile = None,
                    title = '',
                    x_ax_title = '',
                    y_ax_title = '',
                    ):

    ranks = vartbl.copy()

    ranks = ranks.sort_values( by = [ waterfall_col ],
                               ascending = False )

    ranks[ waterfall_col + '_rank' ] = np.arange( ranks.shape[0] )

    fig, ax1 = plt.subplots()

    color = 'blue'
    ax1.set_xlabel( x_ax_title )
    ax1.set_ylabel( y_ax_title, color=color )
    ax1.plot( ranks[ waterfall_col + '_rank' ],
              ranks[ waterfall_col ], color=color)

    if percentile:
        percentile_cutoff = ranks.iloc[ int( ranks.shape[0]*percentile*0.01 ) ][ waterfall_col ]
        plt.axvline( int( ranks.shape[0]*percentile*0.01 ),
                     color = 'black' )
        print( ' To call %d%% of variants as SDV, use the threshold %.2f ' % ( percentile, percentile_cutoff ) )

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

def barplot_per_repeat( table_by_var,
                        sdv_col,
                        nsamp,
                        ylim = None,
                        ylabel = '',
                        xlabel = '',
                        savefig = False,
                        **kwargs ):

    tbv = table_by_var.copy()

    vals = tbv[ sdv_col ].value_counts()

    total = vals.sum()

    xvals = [ i for i in range( 1, nsamp + 1 ) ]

    yvals = []
    for x in xvals:

        if x in vals:
            yvals.append( 100*( int( vals[ x ] ) / total ) )
        else:
            yvals.append( 0 )

    plt.bar( xvals, yvals, **kwargs )

    ax = plt.gca()

    rects = ax.patches

    trans = tfrms.blended_transform_factory( ax.transData, ax.transData )

    for rect, y in zip( rects, yvals ):

        if y > 0:

            ax.text( rect.get_x() + rect.get_width() /2,
                     rect.get_height(),
                     str( int( ( y / 100 )*total ) ),
                     fontsize = 16,
                     transform = trans,
                     ha = 'center',
                     va = 'bottom' )

    if ylim:
        plt.ylim( ylim )

    plt.ylabel( ylabel, fontsize = 18 )
    plt.yticks( fontsize = 18 )

    plt.xlabel( xlabel, fontsize = 18 )
    plt.xticks( fontsize = 18 )

    plt.tight_layout()

    if savefig:
        plt.savefig( savefig,
                     bbox_inches = 'tight' )

    plt.show()

def pr_curves_bytruth( var_df,
                     truth_cols,
                     pred_col,
                     colors,
                   fig_size = ( 5, 5 ),
                   grid = False,
                   x_ax_label = 'Recall\n(%)',
                   y_ax_label = 'Precision\n(%)',
                   add_point = False,
                   savefile = None,
                   **kwargs
             ):

    #pr curve function hates missing values
    tbv = var_df.dropna(subset = [ pred_col ] ).copy()

    if tbv.shape[ 0 ] != var_df.shape[ 0 ]:
        print( 'Missing values in predictor column.', str( var_df.shape[ 0 ] - tbv.shape[ 0 ] ), 'rows removed.' )

    plt.figure( figsize = fig_size )

    for truth, color in zip( truth_cols, colors ):

        #pr curve function hates missing values
        vdf = tbv.dropna(subset = [ truth ] ).copy()

        if vdf.shape[ 0 ] != tbv.shape[ 0 ]:
            print( 'Missing values in truth column.', str( tbv.shape[ 0 ] - vdf.shape[ 0 ] ), 'rows removed.' )

        precision, recall, _ = precision_recall_curve( vdf[ truth ],
                                                       vdf[ pred_col ] )

        plt.plot( 100*recall,
                  100*precision,
                  color = color,
                  **kwargs )

        print( truth, auc( recall, precision ) )

        if add_point:
            plt.plot( add_point[ truth ][ 0 ],
                      add_point[ truth ][ 1 ],
                      color = 'black',
                      marker = add_point[ truth ][ 2 ],
                      markersize = add_point[ truth ][ 3 ]
                    )

    plt.xlabel( x_ax_label, fontsize = 24 )
    plt.xticks( fontsize = 20 )

    plt.ylabel( y_ax_label, fontsize = 24 )
    plt.yticks( fontsize = 20 )

    plt.grid( grid )

    if savefile:
        plt.savefig( savefile,
                     dpi = 300,
                     bbox_inches = 'tight' )

    plt.show()
