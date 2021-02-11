"""Plots simulation and optimization results.

Methods:
    run(conf)
    getData(conf)
    addExtraInfo(row, map_cod_family, map_cod_color, map_cod_marker)
    plot_data(conf, data)
    plot_exp_vs_sim(df, plotSettings, prop_code, name, plotDir)
    add_legend(name, loc, bbox_pos)
"""

# from matplotlib.ticker import MaxNLocator, MultipleLocator, AutoMinorLocator
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import copy
import math
import os

import myExceptions


def run(conf):
    """Plots the following plots:
        - Experimental vs simulation results.
        - Evolution of target function.
        - Evolution of force-field parameters.

    :param conf: (configuration.Conf) Configuration object.
    """

    data = getData(conf)

    plot_data(conf, data)


def getData(conf):
    """Returns experimental data and simulation results read from file.

    :param conf: (configuration.Conf) Configuration object.
    :return:
        data: (pandas DataFrame) Table with data.
    """

    it = conf.it

    fileName = 'ana_{0}/all_sum_{0}.out'.format(it)
    if not os.path.exists(fileName):
        raise myExceptions.NoSuchFile(fileName)

    data = pd.read_csv(fileName, sep='\s+', comment='#',
                       names=['cod', 'frm', 'run', 'pre', 'tem', 'syn1', 'w_dns', 'dns_exp', 'dns_sim', 'dns_dev',
                              'dns_err', 'syn2', 'dns_sim_err', 'dns_unit', 'syn4', 'syn5', 'w_hvp', 'hvp_exp',
                              'hvp_sim', 'hvp_dev', 'hvp_err', 'syn6', 'hvp_sim_err', 'hvp_unit', 'syn8', 'syn9'],
                       usecols=['cod', 'frm', 'run', 'pre', 'tem', 'w_dns', 'dns_exp', 'dns_sim', 'dns_dev', 'dns_err',
                                'dns_sim_err', 'w_hvp', 'hvp_exp', 'hvp_sim', 'hvp_dev', 'hvp_err', 'hvp_sim_err'])

    data = data.replace('*', np.nan)
    data = data.where(data != 0.0, np.nan)

    plotCnf = conf.plotConf

    map_cod_family = plotCnf.map_cod_family
    map_cod_color = plotCnf.map_cod_color
    map_cod_marker = plotCnf.map_cod_marker
    data = data.apply(addExtraInfo, args=(map_cod_family, map_cod_color, map_cod_marker,), axis=1)

    return data


def addExtraInfo(row, map_cod_family, map_cod_color, map_cod_marker):
    """Adds the following information to each row:
         - Number of carbon atoms.
         - First letter code.
         - Short code.
         - Family code.
         - Color used in plot.
         - Marker used in plot..

    :param row: (pandas Series) Row on which operation will be performed.
    :param map_cod_family: (OrderedDict) Dictionary which maps molecule cod to family code.
    :param map_cod_color: (OrderedDict) Dictionary which maps molecule code to colors.
    :param map_cod_marker: (OrderedDict) Dictionary which maps molecule cod to marker.
    """

    row['letter'] = row['cod'][0]
    row['N'] = row['cod'][1]
    row['short_cod'] = '{}_{}'.format(row['cod'][0], row['cod'][2])
    row['family'] = map_cod_family[row['short_cod']]
    row['color'] = map_cod_color[row['short_cod']]
    row['marker'] = map_cod_marker[row['short_cod']]
    return row


def plot_data(conf, data):
    """Plots experimental data versus simulation results.

    :param conf: (configuration.Conf) Configuration object.
    :param data: (pandas DataFrame) Table with data.
    """

    plotCnf = conf.plotConf
    plotDir = plotCnf.plotDir

    plotSettings = plotCnf.settings

    property_codes = [col.replace('_exp', '') for col in data.columns.to_list() if 'exp' in col]

    for prop_code in property_codes:
        df = data[~pd.isnull(data['w_' + prop_code])]
        df = df.replace('-', np.nan)
        df[prop_code + '_exp'] = df[prop_code + '_exp'].astype(float)
        df[prop_code + '_sim'] = df[prop_code + '_sim'].astype(float)

        plot_exp_vs_sim(df, plotSettings, prop_code, 'all', plotDir)

        # Y = 100 * (df[prop_code + '_sim'] - df[prop_code + '_exp']) / df[prop_code + '_exp']
        # Y = df[prop_code + '_sim'] - df[prop_code + '_exp']


def plot_exp_vs_sim(df, plotSettings, prop_code, name, plotDir):
    """Helper function to plot experimental data versus simulation results.

    :param df: (pandas DataFrame) Table with data.
    :param plotSettings: (PlotConf) PlotConf object.
    :param prop_code: (str) Property code (dns, hvp).
    :param name: (str) Name of plot (all).
    :param plotDir: (str) Directory of plots.
    """

    figSizeX = 5.5
    figSizeY = 5.0
    if 'figSizeX' in plotSettings:
        figSizeX = plotSettings['figSizeX']
        figSizeX = float(figSizeX)
    if 'figSizeY' in plotSettings:
        figSizeY = plotSettings['figSizeY']
        figSizeY = float(figSizeY)

    if df.shape[0] == 0:
        return

    fig = plt.figure(figsize=(figSizeX, figSizeY))
    ax = fig.add_subplot(111)

    # Plot data with different markers depending of the family group.
    for idx, row in df.iterrows():
        lab = row['family']
        c = row['color']
        m = row['marker']
        x = row[prop_code + '_exp']
        y = row[prop_code + '_sim']
        plt.plot(x, y, label=lab, color=c, marker=m, markersize=4, linewidth=0.0)

    minX = (min(min(ax.get_xlim()), min(ax.get_ylim()))) - 100
    maxX = (max(max(ax.get_xlim()), max(ax.get_ylim()))) + 100
    if minX < 0:
        minX = 0
    maxX = math.ceil(maxX)
    ax.plot((minX, maxX), (minX, maxX), ls="-", c=".3")

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    delta = 80.0
    if 'delta_' + prop_code in plotSettings:
        delta = plotSettings['delta_' + prop_code]
        delta = float(delta)
    ax.plot((minX, maxX), (minX + delta, maxX + delta), ls="--", c=".3")
    ax.plot((minX, maxX), (minX - delta, maxX - delta), ls="--", c=".3")

    if 'plot_' + prop_code + '_min' in plotSettings:
        minX = float(plotSettings['plot_' + prop_code + '_min'])
        maxX = float(plotSettings['plot_' + prop_code + '_max'])
        step = float(plotSettings['plot_' + prop_code + '_step'])

        ax.set_xlim([minX, maxX])
        ax.set_ylim([minX, maxX])
        ax.xaxis.set_ticks(np.arange(minX, maxX + 1, step))
        ax.yaxis.set_ticks(np.arange(minX, maxX + 1, step))

    if 'add_legend_' + prop_code in plotSettings:
        loc = plotSettings['loc_legend'].replace('__', ' ')
        bbox_to_anchor = plotSettings['bbox_to_anchor']
        boxX = bbox_to_anchor.split('__')[0]
        boxY = bbox_to_anchor.split('__')[1]
        bbox_to_anchor = (float(boxX), float(boxY))

        add_legend(name, loc, bbox_to_anchor)

    # Add letter (a, b, c, d) to left top of plot depending,
    if prop_code == 'dns':
        ax.text(-0.200, 0.975, '(a)', transform=ax.transAxes, size=14)
        # ax.text(-0.200, 0.975, '(c)', transform=ax.transAxes, size=14)
    elif prop_code == 'hvp':
        ax.text(-0.175, 0.975, '(b)', transform=ax.transAxes, size=14)
        # ax.text(-0.175, 0.975, '(d)', transform=ax.transAxes, size=14)

    plt.tight_layout()

    figName = '{}/{}_{}_exp_sim.png'.format(plotDir, name, prop_code)
    fig.savefig(figName)
    plt.close(fig)


def add_legend(name, loc, bbox_pos):
    """Adds legend to plot.

    :param name: (str) Name of plot (all).
    :param loc: (str) Legend location.
    :param bbox_pos: (tuple) Bbox position.
    """

    handles, labels = plt.gca().get_legend_handles_labels()
    newHandles = []
    if name == 'all':
        for h in handles:
            newH = copy.copy(h)
            newH.set_linewidth(5)
            newH.set_markersize(4)
            newHandles.append(newH)

        by_label = OrderedDict(zip(labels, newHandles))

        plt.legend(by_label.values(), by_label.keys(), loc=loc, fontsize=10, numpoints=1, bbox_to_anchor=bbox_pos,
                   borderaxespad=0, frameon=False)

        return

    by_label = OrderedDict(zip(labels, newHandles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=12, numpoints=1, frameon=False)
