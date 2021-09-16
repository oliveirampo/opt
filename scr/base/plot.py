"""Plots simulation and optimization results.

Methods:
    run(conf)
    getData(conf)
    addExtraInfo(row, map_cod_family, map_cod_color, map_cod_marker)
    plot_data(conf, data)
    plot_exp_vs_sim(df, plotSettings, prop_code, name, plotDir)
    add_legend(name, loc, bbox_pos)
    plot_target_function(it, plotConf)
    plot_target_function_helper(actual, predicted, plotDir)
    plot_prm(it, plotConf, plotDir)
    plot_prm_NB(data, plotConf, plotDir)
    getColors()
"""

from matplotlib.ticker import MaxNLocator
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import copy
import math
import sys
import os

from scr.base import myExceptions


def run(conf):
    """Plots the following plots:
        - Experimental vs simulation results.
        - Evolution of target function.
        - Evolution of force-field parameters.

    :param conf: (configuration.Conf) Configuration object.
    """

    it = conf.it
    plotConf = conf.plotConf

    data = getData(it, plotConf)

    # plot_data(plotConf, data)

    # plot_target_function(it, plotConf)
    plot_prm(it, plotConf)


def getData(it, plotConf):
    """Returns experimental data and simulation results read from file.

    :param it: (int) Iteration number.
    :param plotConf: (PlotConf) Plot configuration object.
    :return:
        data: (pandas DataFrame) Table with data.
    """

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

    map_cod_family = plotConf.map_cod_family
    map_cod_color = plotConf.map_cod_color
    map_cod_marker = plotConf.map_cod_marker

    if len(map_cod_family) == 0:
        sys.exit('\n\tERROR: No map of code to family.\n')
    if len(map_cod_color) == 0:
        sys.exit('\n\tERROR: No map of code to color.\n')
    if len(map_cod_marker) == 0:
        sys.exit('\n\tERROR: No map of code to marker.\n')

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


def plot_data(plotConf, data):
    """Plots experimental data versus simulation results.

    :param plotConf: (PlotConf) Plot configuration object.
    :param data: (pandas DataFrame) Table with data.
    """

    plotDir = plotConf.plotDir

    plotSettings = plotConf.settings

    # property_codes = [col.replace('_exp', '') for col in data.columns.to_list() if 'exp' in col]
    property_codes = [col.replace('_exp', '') for col in data.columns.tolist() if 'exp' in col]

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
    :param plotSettings: (Dictionary) Settings for plot.
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

    propLabel = plotSettings[prop_code + '_label'].replace('__', ' ')
    propUnit = plotSettings[prop_code + '_unit'].replace('__', ' ')
    ax.set_ylabel('{} {}'.format(propLabel, propUnit), fontsize=18)
    propLabel = propLabel.replace('sim', 'exp')
    ax.set_xlabel('{} {}'.format(propLabel, propUnit), fontsize=18)

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


def plot_target_function(it, plotConf):
    """Plots evolution of target function.

    :param it: (int) Iteration Number.
    :param plotConf: plotConf: (PlotConf) Plot configuration object.
    :return:
    """

    plt.rcParams['mathtext.fontset'] = 'cm'

    plotDir = plotConf.plotDir

    data = {}
    for i in range(1, it + 1):
        optFile = 'opt_{0}/opt_{0}.out'.format(i)
        if not os.path.exists(optFile):
            raise myExceptions.NoSuchFile(optFile)

        d = pd.read_csv(optFile, sep='\s+', comment='#', names=['step', 'value'])
        data[i] = d.copy()

    actual = []
    predicted = []
    for i in data.keys():
        ini = data[i]['value'].iloc[0]
        actual.append([int(i) - 1, ini])
    for i in list(data.keys()):
        fin = data[i]['value'].iloc[-1]
        predicted.append([i, fin])

    actual = np.asarray(actual)
    predicted = np.asarray(predicted)
    plot_target_function_helper(actual, predicted, plotDir)


def plot_target_function_helper(actual, predicted, plotDir):
    """Helper function to plot evolution of target function.

    :param actual: (array) Actual values of target function.
    :param predicted: (array) Predicted values of target function.
    :param plotDir: (str) Directory where plot will be saved.
    :return:
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.plot(actual[:, 0], actual[:, 1], color='red',  marker='o', label=r'$Q_i^{real}$', markersize=8, linewidth=3.0)
    plt.plot(predicted[:, 0], predicted[:, 1], color='blue', marker='o', label=r'$Q_i^{pred}$')

    ax.set_xlabel(r'$i$', fontsize=14)
    ax.set_ylabel(r'$Q_i$', fontsize=14)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # plt.title('HAL-CAL')

    ax.legend(loc='upper right', frameon=False, fontsize=14)
    plt.tight_layout()

    file_name = '{}/target_function.png'.format(plotDir)
    fig.savefig(file_name)
    plt.close(fig)


def plot_prm(it, plotConf):
    """Plot evolution of force-field parameter along iteration number.

    :param it: (int) Iteration Number.
    :param plotConf: plotConf: (PlotConf) Plot configuration object.
    :return:
    """

    plotDir = plotConf.plotDir

    data = pd.DataFrame()

    it = int(it)
    for i in range(it + 1):
        prmFile = '00_inp/prm_' + str(i) + '.dat'
        if not os.path.exists(prmFile):
            raise myExceptions.NoSuchFile(prmFile)

        d = pd.read_csv(prmFile, sep='\s+', comment='#',
                        names=['iac', 'typ', 'name', 'sig', 'rng_sig', 'eps', 'rng_eps', 'hrd', 'rng_hrd', 'eln', 'rng_eln',
                               'sig_2', 'rng_sig_2', 'eps_2', 'rng_eps_2'],
                        usecols=['iac', 'name', 'sig', 'eps', 'hrd', 'eln'])

        d['it'] = i
        data = data.append(d)

    # for iac in [13, 14, 15, 16]:
    #     data = data[data.iac != iac]

    plot_prm_NB(data, plotConf, plotDir)


def plot_prm_NB(data, plotConf, plotDir):
    """Helper function to plot evolution of force-field parameters.

    :param data: (pandas DataFrame) Table with data.
    :param plotConf: plotConf: (PlotConf) Plot configuration object.
    :param plotDir: (str) Directory where plot will be saved.
    :return:
    """

    plt.rcParams['mathtext.fontset'] = 'cm'
    color = getColors()

    prmIni = 0
    prmFin = 2
    nPrms = prmFin - prmIni
    prmTypes = ['sig', 'eps', 'hrd', 'eln', 'eln']
    prmNames = [r'$\sigma \,$ [nm]', r'$\epsilon \,$ [kJ $\cdot$ mol$^{-1}$]', r'$\eta \,$ [$e^{-1} \cdot$ V]',
                r'$\chi \,$ [V]']
    # prmNames = [r'$\tilde{\sigma} \,$ [nm]', r'$\tilde{\epsilon} \,$ [kJ $\cdot$ mol$^{-1}$]',
    #             r'$\eta \,$ [$e^{-1} \cdot$ V]', r'$\chi \,$ [V]']

    if nPrms == 4:
        fig, axes = plt.subplots(nrows=1, ncols=nPrms, figsize=(12, 8))
    elif nPrms == 2:
        fig, axes = plt.subplots(nrows=1, ncols=nPrms, figsize=(10, 8))
    else:
        print('\n\tnPrms = 2 or 4.')
        sys.exit(1)

    iac = data.iac.unique()
    iacName = iac

    if len(plotConf.map_iac_name) == 0:
        print('Define map from iac code to name.')
        sys.exit(1)

    iacName = plotConf.map_iac_name  # dict

    idx = 0
    for i in iacName:
        c = color[idx]
        nam = iacName[i]
        idx += 1

        i = int(i)
        X = data[data['iac'] == i]['it']
        for j in range(nPrms):
            Y = data[data['iac'] == i][prmTypes[j + prmIni]]
            if prmTypes[j + prmIni] in ['eln', 'hrd']:
                Y = 1.44 * Y
            axes[j].plot(X, Y, marker='o', markersize=4, color=c, label=nam)

    for j in range(nPrms):
        axes[j].grid()
        axes[j].set_ylabel(prmNames[j + prmIni], size=16)
        axes[j].xaxis.set_major_locator(MaxNLocator(integer=True))
        axes[j].tick_params(axis="x", labelsize=12)
        axes[j].tick_params(axis="y", labelsize=12)

        # axes[j].xaxis.set_ticks(np.arange(0, 7, 2))

    if nPrms == 4:
        axes[nPrms - 1].legend(loc='upper right', fontsize=16, bbox_to_anchor=(1.875, 1.0), frameon=False)
        plt.subplots_adjust(wspace=0.5)
        fig.text(0.475, 0.02, r'$i$', ha='center', size=16)
        fig.subplots_adjust(left=0.10, right=0.875)

    elif nPrms == 2:
        # LJ
        # axes[nPrms - 1].legend(loc='upper right', fontsize=12, bbox_to_anchor=(1.45, 1.01), frameon=False)
        # fig.text(0.45, 0.020, 'Iteration', ha='center', size=16)
        # fig.subplots_adjust(left=0.10, right=0.85)
        # plt.subplots_adjust(wspace=0.30)

        # EE
        # axes[nPrms - 1].legend(loc='upper right', fontsize=12, bbox_to_anchor=(1.55, 1.01), frameon=False)
        # fig.text(0.45, 0.020, 'Iteration', ha='center', size=16)
        # plt.subplots_adjust(wspace=0.30)

        # OXY + HB EEM
        plt.subplots_adjust(wspace=0.25)
        axes[nPrms - 1].legend(loc='upper right', fontsize=12, bbox_to_anchor=(1.50, 1.01), frameon=False)
        fig.subplots_adjust(left=0.0825, right=0.835)
        # axes[0].set_ylim([0, 1])
        # axes[1].set_ylim([0, 1])
        # axes[1].set_yticks(np.arange(0, 56, 5))

        # OXY + HB LJ
        # plt.subplots_adjust(wspace=0.25)
        # axes[nPrms - 1].legend(loc='upper right', fontsize=12, bbox_to_anchor=(1.400, 1.00), frameon=False)
        # fig.subplots_adjust(left=0.0825, right=0.85)
        # axes[0].set_ylim([0.22, 0.36])
        # axes[1].set_ylim([0.10, 0.90])

    file_name = '{}/prm_NB.png'.format(plotDir)
    fig.savefig(file_name)
    plt.close(fig)


def getColors():
    """Returns array with color names."""

    color = [
        'red',
        'green',
        'blue',
        'olive',
        'gray',
        'coral',
        'navy',
        'maroon',
        'magenta',
        'deeppink',
        'peru',
        'deepskyblue',
        'darkcyan',
        'darkorchid',
        'orangered',
        'lime',
        'aqua',
        'darkseagreen',
        'teal',
        'crimson',
        'turquoise',
        'gold',
        'yellow',
        'greenyellow',
        'fuchsia',
        'mediumslateblue',
        'darkorange',
        'aliceblue',
        'antiquewhite',
        'aquamarine',
        'bisque',
        'blueviolet',
        'blueviolet',
        'darkgray',
        'darkmagenta',
        'darkolivegreen',
        'darkred',
        'gainsboro',
        'honeydew',
        'indianred',
        'khaki',
        'linen',
        'rosybrown',
        'moccasin',
        ]
    return color
