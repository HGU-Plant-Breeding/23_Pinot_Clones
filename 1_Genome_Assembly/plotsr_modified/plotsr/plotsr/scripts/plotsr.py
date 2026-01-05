#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Script Name: plotsr.py
Description: Main logic for plotting multi-genome structural annotations.
             Orchestrates data loading, filtering, and calling plotting functions.
Author: Manish Goel (Original script modified by Paolo Callipo)
Date: 30.12.2021 (Updated 2025)
"""

import argparse
import sys
import os
import logging
from collections import deque, OrderedDict
from math import ceil

# Third-party imports
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from pandas import concat as pdconcat
from pandas import unique

# Local imports
from plotsr import __version__
from plotsr.scripts.func import (
    setlogconfig, readbasecfg, readsyriout, readbedout,
    filterinput, validalign2fasta, selectchrom,
    selectregion, createribbon, drawax, pltchrom,
    pltsv, drawmarkers, readtrack, drawtracks,
    getfilehandler, definelogger,
    calculate_snp_density, get_density_statistics
)

def plotsr(args):
    ## Define loggers
    setlogconfig(args.log)
    filehandler = getfilehandler(args.logfin.name, args.log)
    global getlogger
    getlogger = definelogger(filehandler)
    logger = getlogger("Plotsr")

    ###################################################################
    # Check python version
    ###################################################################
    logger.debug('checking arguments')
    try:
        assert sys.version_info.major == 3
        assert sys.version_info.minor >= 8
    except AssertionError:
        logger.warning(f'\nPlotsr is tested for Python >=3.8. Currently using Python {sys.version_info.major}.{sys.version_info.minor}. This may result in errors.')
    except KeyboardInterrupt:
        raise()
    except Exception as E:
        sys.exit(E)

    ## Validate input
    if args.sr is None and args.bp is None:
        logger.error("No structural annotations provided. Use --sr or -bp to provide path to input files")
        sys.exit()

    if args.sr is not None and args.bp is not None:
        logger.error("Both --sr and --bp cannot be used. Use single file type for all input structural annotations files.")
        sys.exit()

    if args.chr is not None and args.reg is not None:
        logger.error("Both --chr and --reg are provided. Only one parameter can be provided at a time. Exiting.")
        sys.exit()

    if args.chr is not None and args.chrord is not None:
        logger.error("Both --chr and --chrord are provided. Only one parameter can be provided at a time. Exiting.")
        sys.exit()

    if args.rtr and args.reg is None:
        logger.error("Cannot use --rtr without --reg. Exiting.")
        sys.exit()

    # Check if SNP density arguments are valid
    if args.snp_density and (args.snp_norm_max_perc < 0 or args.snp_norm_max_perc > 100):
        logger.error("--snp-norm-max-perc must be between 0 and 100. Exiting.")
        sys.exit()

    ###################################################################
    # Declare variable using argument values
    ###################################################################
    logger.info('Starting')
    FS = args.f             # Font size
    H = args.H              # Height
    W = args.W              # Width
    O = args.o              # Output file name
    D = args.d              # Output file DPI
    R = args.R              # Create ribbons
    V = args.v              # Vertical chromosomes
    S = args.S              # Space between homologous chromosomes
    B = None if args.markers is None else args.markers.name
    TRACKS = None if args.tracks is None else args.tracks.name
    REG = None if args.reg is None else args.reg.strip().split(":")
    RTR = args.rtr
    CHRS = args.chr
    ITX = args.itx
    CHRNAME = args.chrname.name if args.chrname is not None else None

    ## Get config
    cfg = readbasecfg('', V) if args.cfg is None else readbasecfg(args.cfg.name, V)
    if S < 0.1 or S > 0.75:
        logger.warning('Value for S outside of normal range 0.1-0.75.')

    ## Check output file extension
    if len(O.split('.')) == 1:
        logger.warning("Output filename has no extension. Plot would be saved as a pdf")
        O = O + ".pdf"
    elif O.split('.')[-1] not in ['pdf', 'png', 'svg']:
        logger.warning("Output file extension is not in {'pdf','png', 'svg'}. Plot would be saved as a pdf")
        O = O.rsplit(".", 1)[0] + ".pdf"

    ## Set matplotlib backend
    try:
        matplotlib.use(args.b)
    except:
        sys.exit('Matplotlib backend cannot be selected')

    # Read alignment coords
    alignments = deque()
    chrids = deque()
    if args.sr is not None:
        for f in args.sr:
            fin = f.name
            al, cid = readsyriout(fin)
            alignments.append([os.path.basename(fin), al])
            chrids.append((os.path.basename(fin), cid))
    elif args.bp is not None:
        for f in args.bp:
            fin = f.name
            al, cid = readbedout(fin)
            alignments.append([os.path.basename(fin), al])
            chrids.append((os.path.basename(fin), cid))

    # Get groups of homologous chromosomes
    cs = set(unique(alignments[0][1]['achr']))
    if args.chrord is None:
        chrs = [k for k in chrids[0][1].keys() if k in alignments[0][1]['achr'].unique()]
    else:
        chrs = deque()
        with open(args.chrord.name, 'r') as fin:
            for line in fin:
                c = line.strip()
                if c not in cs:
                    logger.error("Chromosome {} in {} is not a chromosome in alignment file {}. Exiting.".format(c, args.chrord.name, alignments[0][0]))
                    sys.exit()
                chrs.append(c)
        chrs = list(chrs)
        if len(chrs) != len(cs):
            logger.error("Number of chromosomes in {} is less than the number of chromosomes in the alignment file {}. Exiting.".format(args.chrord.name, alignments[0][0]))
            sys.exit()

    chrgrps = OrderedDict()
    for c in chrs:
        cg = deque([c])
        cur = c
        for i in range(len(chrids)):
            n = chrids[i][1][cur]
            cg.append(n)
            cur = n
        chrgrps[c] = cg

    # Check chromosome IDs and sizes
    chrlengths, genomes = validalign2fasta(alignments, args.genomes.name)

    # Calculate SNP density if file is provided
    snp_density_data = None
    snp_norm_max = None
    if args.snp_density:
        logger.info("Calculating SNP density...")
        snp_density_data = calculate_snp_density(args.snp_density.name,
                                                 args.snp_window_size,
                                                 chrlengths)
        snp_norm_max = get_density_statistics(snp_density_data, args.snp_norm_max_perc)
        logger.info(f"SNP density color scale maximum set to: {snp_norm_max:.2f} (the {args.snp_norm_max_perc}th percentile).")

    # Filter alignments
    for i in range(len(alignments)):
        alignments[i][1] = filterinput(args, alignments[i][1], chrids[i][1], ITX)

    # Select only chromosomes selected by --chr
    if CHRS is not None:
        alignments, chrs, chrgrps, chrlengths = selectchrom(CHRS, cs, chrgrps, alignments, chrlengths, chrids)

    if REG is not None:
        alignments, chrs, chrgrps = selectregion(REG, RTR, chrlengths, alignments, chrids)

    # Combine Ribbon if selected
    if R:
        for i in range(len(alignments)):
            alignments[i][1] = createribbon(alignments[i][1])

    # Invert coord for inverted query genome
    for i in range(len(alignments)):
        df = alignments[i][1].copy()
        invindex = ['INV' in i for i in df['type']]
        g = set(df.loc[invindex, 'bstart'] < df.loc[invindex, 'bend'])
        if len(g) == 2:
            logger.error("Inconsistent coordinates in input file {}. Mixing inverted coordinates is not permitted.".format(alignments[i][0]))
            sys.exit()
        elif False in g:
            continue
        df.loc[invindex, 'bstart'] = df.loc[invindex, 'bstart'] + df.loc[invindex, 'bend']
        df.loc[invindex, 'bend'] = df.loc[invindex, 'bstart'] - df.loc[invindex, 'bend']
        df.loc[invindex, 'bstart'] = df.loc[invindex, 'bstart'] - df.loc[invindex, 'bend']
        alignments[i][1] = df.copy()

    plt.rcParams['font.size'] = FS
    try:
        if H is None and W is None:
            H = len(chrs)
            W = 3
            fig = plt.figure(figsize=[W, H])
        elif H is not None and W is None:
            fig = plt.figure(figsize=[H, H])
        elif H is None and W is not None:
            fig = plt.figure(figsize=[W, W])
        else:
            fig = plt.figure(figsize=[W, H])
    except Exception as e:
        logger.error("Error in initializing figure.\n{}".format(e.with_traceback()))
        sys.exit()
    ax = fig.add_subplot(111, frameon=False)

    allal = pdconcat([alignments[i][1] for i in range(len(alignments))])
    if ITX:
        minl = 0
        MCHR = cfg['marginchr']
        maxchr = max([sum(chrlengths[i][1].values()) for i in range(len(chrlengths))])
        maxl = int(maxchr/(MCHR + 1 - (MCHR*len(chrgrps))))
    elif cfg['maxl'] != -1:
        minl, maxl = 0, cfg['maxl']
    elif REG is None:
        minl, maxl = 0, -1
    else:
        minl = min(allal[['astart', 'bstart']].apply(min))
        maxl = max(allal[['aend', 'bend']].apply(max))

    labelcnt = 0
    if 'SYN' in allal['type'].array: labelcnt += 1
    if 'INV' in allal['type'].array: labelcnt += 1
    if 'TRA' in allal['type'].array or 'INVTR' in allal['type'].array: labelcnt += 1
    if 'DUP' in allal['type'].array or 'INVDP' in allal['type'].array: labelcnt += 1

    ## Draw Axes
    ax = drawax(ax, chrgrps, chrlengths, V, S, cfg, ITX, minl=minl, maxl=maxl, chrname=CHRNAME)

    ## Draw Chromosomes
    ax, indents, chrlabels = pltchrom(ax, chrs, chrgrps, chrlengths, V, S, genomes, cfg, ITX,
                                      minl=minl, maxl=maxl,
                                      snp_density_data=snp_density_data,
                                      snp_window_size=args.snp_window_size,
                                      snp_colormap=args.snp_colormap,
                                      snp_norm_max=snp_norm_max,
                                      snp_bar_thickness=args.snp_bar_thickness)

    if cfg['genlegcol'] < 1:
        ncol = ceil(len(chrlengths)/labelcnt)
    else:
        ncol = int(cfg['genlegcol'])

    # --- Plot structural annotations ---
    ax, svlabels = pltsv(ax, alignments, chrs, V, chrgrps, chrlengths, indents, S, cfg, ITX, maxl)

    # --- Legends ---
    if cfg['legend']:
        bbox_to_anchor = cfg['bbox']
        
        # 1. Genome Legend (Skip if SNP Density or ITX is active)
        genomes_legend_drawn = False
        if not snp_density_data and not ITX:
            l1 = plt.legend(handles=chrlabels, loc='lower left', bbox_to_anchor=bbox_to_anchor, 
                            ncol=ncol, mode=None, borderaxespad=0., frameon=False, title='Genomes')
            l1._legend_box.align = "left"
            plt.gca().add_artist(l1)
            genomes_legend_drawn = True

        # 2. Annotation Legend
        # Only shift right if the Genome legend was actually drawn
        if genomes_legend_drawn:
            bbox_to_anchor[0] += cfg['bboxmar']
            
        plt.legend(handles=svlabels, loc='lower left', bbox_to_anchor=bbox_to_anchor, 
                   ncol=1, mode='expand', borderaxespad=0., frameon=False, title='Annotations')._legend_box.align = "left"

    # Plot markers
    if B is not None:
        ax = drawmarkers(ax, B, V, chrlengths, indents, chrs, chrgrps, S, cfg, ITX, minl=minl, maxl=maxl)

    # Draw tracks
    if TRACKS is not None:
        tracks = readtrack(TRACKS, chrlengths)
        ax = drawtracks(ax, tracks, S, chrgrps, chrlengths, V, ITX, cfg, minl=minl, maxl=maxl)

    # Add color bar for SNP density
    if snp_density_data:
        cax_position = [0.15, 0.9, 0.3, 0.025] # [left, bottom, width, height]
        cax = fig.add_axes(cax_position)

        cmap = cm.get_cmap(args.snp_colormap)
        norm = mcolors.Normalize(vmin=0, vmax=snp_norm_max)
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

        cb.set_label(f'SNP Density (per {int(args.snp_window_size/1000)}kb)', size='small')
        cax.xaxis.set_label_position('top')
        cb.ax.tick_params(labelsize='x-small')
        cax.xaxis.set_ticks_position('top')

        ticks = cb.get_ticks()
        labels = [f'{int(t)}' for t in ticks]
        if labels: labels[-1] = f'{int(ticks[-1])}+'
        cb.set_ticks(ticks)
        cb.set_ticklabels(labels)

    # Save the plot
    try:
        pad_inches = 0.1 if snp_density_data else 0.01
        fig.savefig(O, dpi=D, bbox_inches='tight', pad_inches=pad_inches)
        logger.info("Plot {O} generated.".format(O=O))
    except Exception as e:
        sys.exit('Error in saving the figure. Try using a different backend.' + '\n' + str(e))
    logger.info('Finished')

def main():
    from matplotlib.rcsetup import non_interactive_bk as bklist
    parser = argparse.ArgumentParser("Plotting structural rearrangements between genomes", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    other = parser._action_groups.pop()
    inputfiles = parser.add_argument_group("Input/Output files")
    inputfiles.add_argument('--sr', help='Structural annotation mappings (syri.out) identified by SyRI', action='append', type=argparse.FileType('r'))
    inputfiles.add_argument('--bp', help='Structural annotation mappings in BEDPE format', action='append', type=argparse.FileType('r'))
    inputfiles.add_argument('--genomes', help='File containing path to genomes', type=argparse.FileType('r'), required=True)
    inputfiles.add_argument('--markers', help='File containing path to markers (bed format)', type=argparse.FileType('r'))
    inputfiles.add_argument('--tracks', help='File listing paths and details for all tracks to be plotted', type=argparse.FileType('r'))
    inputfiles.add_argument('--chrord', help='File containing reference (first genome) chromosome IDs in the order in which they are to be plotted.', type=argparse.FileType('r'))
    inputfiles.add_argument('--chrname', help='File containing reference (first genome) chromosome names to be used in the plot.', type=argparse.FileType('r'))
    inputfiles.add_argument('-o', help='Output file name. Acceptable format: pdf, png, svg', default="plotsr.pdf")

    filtering = parser.add_argument_group("Data filtering")
    filtering.add_argument('--itx', help='Use inter-chromosomal plotting mode (experimental)', default=False, action='store_true')
    filtering.add_argument('--chr', help='Select specific chromosome on reference (first genome).', type=str, action='append')
    filtering.add_argument('--reg', help='Plots a specific region. GenomeID:ChromosomeID:Start-End.', type=str)
    filtering.add_argument('--rtr', help='Plot all SRs within boundaries of homologous regions.', default=False, action='store_true')
    filtering.add_argument('--nosyn', help='Do not plot syntenic regions', default=False, action='store_true')
    filtering.add_argument('--noinv', help='Do not plot inversions', default=False, action='store_true')
    filtering.add_argument('--notr', help='Do not plot translocations regions', default=False, action='store_true')
    filtering.add_argument('--nodup', help='Do not plot duplications regions', default=False, action='store_true')
    filtering.add_argument('-s', help='minimum size of a SR to be plotted', type=int, default=10000)

    plotting = parser.add_argument_group("Plot adjustment")
    plotting.add_argument('--cfg', help='Path to config file containing parameters to adjust plot.', type=argparse.FileType('r'))
    plotting.add_argument('-R', help='Join adjacent syntenic blocks.', default=False, action="store_true")
    plotting.add_argument('-f', help='font size', type=int, default=6)
    plotting.add_argument('-H', help='height of the plot', type=float)
    plotting.add_argument('-W', help='width of the plot', type=float)
    plotting.add_argument('-S', help='Space for homologous chromosome (0.1-0.75).', default=0.7, type=float)
    plotting.add_argument('-d', help='DPI for the final image', default="300", type=int)
    plotting.add_argument('-b', help='Matplotlib backend to use', default="agg", type=str, choices=bklist)
    plotting.add_argument('-v', help='Plot vertical chromosome', default=False, action='store_true')

    snp_heatmap = parser.add_argument_group("SNP Density Heatmap")
    snp_heatmap.add_argument('--snp-density', help='Path to SNP data file to draw density heatmaps instead of solid chromosome bars. Format: GenomeID<TAB>ChrID<TAB>Position.', type=argparse.FileType('r'))
    snp_heatmap.add_argument('--snp-window-size', help='Window size in base pairs for SNP density calculation.', type=int, default=100000)
    snp_heatmap.add_argument('--snp-bar-thickness', help='Thickness of the SNP density heatmap bars.', type=float, default=0.08)
    snp_heatmap.add_argument('--snp-norm-max-perc', help='Set the color scale maximum to this percentile of the data.', type=float, default=98.0)
    snp_heatmap.add_argument('--snp-colormap', help='Matplotlib colormap to use for the SNP density heatmap.', type=str, default='viridis')

    other.add_argument("--lf", dest="logfin", help="Name of log file", type=argparse.FileType("w"), default="plotsr.log")
    other.add_argument('--log', help='Log-level', choices=['DEBUG', 'INFO', 'WARN'], default='WARN', type=str)
    other.add_argument('--version', action='version', version='{version}'.format(version=__version__))
    parser._action_groups.append(other)

    args = parser.parse_args()
    plotsr(args)

if __name__ == '__main__':
    main()