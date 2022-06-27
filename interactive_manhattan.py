#!/usr/bin/env python3
#
# This file is part of the XXX distribution (https://github.com/LudvigEk/interactive_manhattan).
# Copyright (c) 2022 Ludvig Ekdahl.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import pandas as pd
import numpy as np
import sys
import argparse


def make_manhattan(input_file, output_file, refgenome, color1='#e69d00', color2='#0071b2', UCSCcolor="#000000", title="manhattan"):
    """Function to create interactive manhattan plots.

	Creates manhattan plots using bokeh that shows chromosome, position, pvalue and rsid on hovering over snps.
	When clicking points in the plot, the user is brought to UCSC genome browser at that location with a highlight on the snp.
	To be easier on the browser script engine the manhattan i downsampled prior to exporting the html.
	SNPs are starting to be filtered above p-values of 0.00005 (0.05 * 10e-3)

	Example:

	        $ python interactive_manhattan.py silly_height_gwas.txt height_gwaS_manhattan.html hg38

		To see full usage options use

			$ python interactive_manhattan.py -h


	Attributes:
	    input_file (str): file path to gwas track, required columns are: 'pval', 'chromosome', 'position', 'rsid'
		
		output_file (str): file path to output html file

		refgenome (str): coordinate system of gwas track, should be either of 'hg19' or 'hg38'
			NOTE : NO CONVERSION IS HAPPENING IN THE SCRIPT, THIS IF FOR VISUALISATION PURPOSES ONLY

		color1 (str): color of the snps that are located on uneven chromsomes

		color2 (str): color of the snps that are located on even chromsomes

		UCSCcolor (str): color of the UCSC highlight bar, default is black

	"""

    if refgenome not in ["hg19", "hg38"]:
        reportStr = "ERROR: refgenome must be specified, either of 'hg19' or 'hg38'"

    hexcolor = UCSCcolor[1:]  # remove prepended hashtag

    # input gwas track, has to have columns 'pval', 'chromosome' and 'position
    gwas_file = str(input_file)
    # Output interactive manhattan
    output_file_path = str(output_file)

    # chunk it
    # chunksize = 10 ** 6
    # for chunk in pd.read_csv(gwas_file, chunksize=chunksize, sep="\t"):
    #	gwas_track = chunk
    #	break
    gwas_track = pd.read_csv(gwas_file, sep="\t")

    # Assert correct labels for neccessary positions´
    assert all([x in gwas_track.columns.tolist() for x in ['pval', 'chromosome', 'position', 'rsid']])

    # Create a -log10(p) col
    gwas_track['log_p'] = -np.log10(gwas_track['pval'])

    # Down sampling, select every 100th point over 0.0005, every 50th point over 0.00005
    first_filter = gwas_track[(gwas_track['pval'] > 0.005)].index
    i_to_keep = []
    for i in np.arange(0, len(first_filter), 1):
        if i % 100 == 0:
            i_to_keep.append(i)
    first_filter = first_filter[i_to_keep]
    # 2nd
    second_filter = gwas_track[(gwas_track['pval'] < 0.005) & (gwas_track['pval'] > 0.0005)].index
    i_to_keep = []
    for i in np.arange(0, len(second_filter), 1):
        if i % 50 == 0:
            i_to_keep.append(i)
    second_filter = second_filter[i_to_keep]
    # 3rd
    third_filter = gwas_track[(gwas_track['pval'] < 0.0005) & (gwas_track['pval'] > 0.00005)].index
    i_to_keep = []
    for i in np.arange(0, len(third_filter), 1):
        if i % 20 == 0:
            i_to_keep.append(i)
    third_filter = third_filter[i_to_keep]
    # Keep all below
    keep_all = gwas_track[(gwas_track['pval'] < 0.00005)].index

    final_filter = first_filter.union(second_filter)
    final_filter = final_filter.union(keep_all)
    final_filter = final_filter.union(third_filter)
    assert (len(first_filter) + len(second_filter) + len(keep_all) + len(third_filter)) == len(final_filter)

    gwas_track_downsampled = gwas_track.loc[final_filter]

    reportStr = "Manhattan down sampled to " + str(len(gwas_track_downsampled)) + " pvalues\n"
    sys.stderr.write(reportStr)

    # Create a global position index
    gwas_track_downsampled['global_pos'] = 0
    to_add = 0
    chrom_midpoints = []
    for chrom in np.arange(1, 23, 1):
        # print(chrom)
        chr_idx = gwas_track_downsampled[gwas_track_downsampled['chromosome'] == chrom].index
        gwas_track_downsampled.loc[chr_idx, 'global_pos'] = gwas_track_downsampled.loc[chr_idx, 'position'] + to_add
        midPoint = int( (max(gwas_track_downsampled.loc[chr_idx, 'global_pos']) +
                         min(gwas_track_downsampled.loc[chr_idx, 'global_pos']) ) / 2)
        chrom_midpoints.append(midPoint)
        to_add = max(gwas_track_downsampled.loc[chr_idx, 'global_pos']) + 30000000

    from bokeh.plotting import figure, output_file, show
    from bokeh.models import ColumnDataSource, CustomJS, HoverTool, TapTool
    from bokeh.models.tickers import FixedTicker
    from bokeh.transform import factor_cmap
    from bokeh.io import save

    # Create nice alternating palette for chromosomes
    palette = []

    for i in np.arange(1, 25, 1):
        if i % 2 == 0:
            color = color1
        else:
            color = color2
        palette.append(color)

    print(palette)

    # Init plotting object
    source = ColumnDataSource(data=dict(x=gwas_track_downsampled['global_pos'],
                                        y=gwas_track_downsampled['log_p'],
                                        rsid=gwas_track_downsampled['rsid'],
                                        pos=gwas_track_downsampled['position'],
                                        pval=gwas_track_downsampled['pval'],
                                        chrom=gwas_track_downsampled['chromosome'].astype(str)))

    # Create javascript callback
    callback_obj = CustomJS(args=dict(source=source, refgenome=refgenome, color=hexcolor),
                            code="""
	
	        // Select the data
	        var inds = source.selected.indices;
	        var d1 = source.data;
	        var idx = []
	        var chrom = []
	        var pos = []
	        for (var i = 0; i < inds.length; i++) {
	            idx.push(d1['rsid'][inds[i]])
	            chrom.push(d1['chrom'][inds[i]])
	            pos.push(d1['pos'][inds[i]])
	        }
	
    
	        if (pos[0] !== undefined || chrom[0] !== undefined) {

	            var url = sprintf("https://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=chr%s%%3A%s-%s&highlight=%s.chr%s%%3A%s-%s%%23%s", refgenome, chrom[0], pos[0]-500, pos[0]+500, refgenome, chrom[0], pos[0], pos[0], color); //AA0000
	            window.open(url);
	        }
	        //Advanced UCSC link for reference with two highlights
	        //https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr9%3A136136597-136139844&highlight=hg19.chr9%3A136138630-136139650%23AA0000%7Chg19.chr9%3A136136630-136137650%23AA0000
	        
	""")

    output_file(output_file_path)

    p = figure(title=title, plot_width=1200, plot_height=600, tools="tap,box_select,box_zoom,reset,help")

    colors = factor_cmap('chrom', palette=palette, factors=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"])

    p.xaxis.axis_label = "Chromosome"
    p.yaxis.axis_label = "-log10(P-value)"

    p.xaxis.ticker = FixedTicker(ticks=chrom_midpoints)
    label_dict = {}
    i = 1
    for tick in chrom_midpoints:
        label_dict[tick] = str(i)
        i += 1
    p.xaxis.major_label_overrides = label_dict

    p.scatter(source=source,
              x='x',
              y='y',
              fill_color=colors,
              line_color=colors)

    taptool = p.select(type=TapTool)
    taptool.callback = callback_obj

    p.add_tools(HoverTool(tooltips=[("rsid", "@rsid"), ('chrom', '@chrom'), ('pos', '@pos'), ('pval', '@pval')]))

    reportStr = "Saving results to " + output_file_path + "\n"
    sys.stderr.write(reportStr)
    save(p, output_file_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--evencol", help="Manhattan color for even chromosomes, web format; #e69d00",
                        default="#e69d00", type=str)
    parser.add_argument("-u", "--unevencol", help="Manhattan color for uneven chromosomes, web format; #0071b2",
                        default="#0071b2", type=str)
    parser.add_argument("-hl", "--highlight", help="Color for highlight in UCSC browser", default="#000000", type=str)
    parser.add_argument("-t", "--title", help="Title for the manhattan plot", default="manhattan", type=str)
    parser.add_argument("<input>", help="path to input table", type=str )
    parser.add_argument("<output>", help="path to output html file", type=str)
    parser.add_argument("<ref>",
                        help="Coordinates of input data ('hg19'/'hg38') WARNING: there's no checks or conversion "
                             "happening through this command (!), it affects UCSC links", type=str)
    args = vars(parser.parse_args())

    make_manhattan(args["<input>"], args["<output>"], args["<ref>"], color1=args["evencol"], color2=args["unevencol"],
                   UCSCcolor=args["highlight"], title=args["title"])
