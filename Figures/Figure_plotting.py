import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from neo.rawio import AxonRawIO
import math
import os

'Fonts/font size may be slightly off depending on the plot. Font size can be adjusted.'

def fig5a():
    data = pd.read_excel('Fig5a5b.xlsx')
    data = data.drop('avg_logFC', axis=1)
    data = data.drop('pct.1', axis=1)
    data = data.drop('pct.2', axis=1)
    data = data.drop('p_val_adj', axis=1)
    data = data.set_index('Gene')
    data = data.loc[:, ~data.columns.str.contains('^Unnamed')]
    print(data)
    aldoc_data = data.loc[data['Aldoc+'] > data['Plcb4+']]
    plcb_data = data.loc[data['Plcb4+'] > data['Aldoc+']]
    fig, ax = plt.subplots(1,2, figsize=(1.5,8))
    sns.heatmap(aldoc_data,cmap='Blues', ax=ax[0], vmin=0,
                cbar_kws={'orientation': 'horizontal'})
    sns.heatmap(plcb_data,cmap='Reds', ax=ax[1], vmin=0,
                cbar_kws={'orientation': 'horizontal'})
    ax[0].set_frame_on(True)
    ax[0].get_yaxis().set_ticks([])
    ax[0].set_ylabel(' ')
    ax[0].get_xaxis().set_visible(False)
    ax[1].set_frame_on(True)
    ax[1].get_yaxis().set_ticks([])
    ax[1].set_ylabel(' ')
    ax[1].get_xaxis().set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'Fig5a.pdf'), bbox_inches='tight')
    return fig, ax

def fig5b():
    data = pd.read_excel('Fig5a5b.xlsx')
    data = data.drop('avg_logFC', axis=1)
    data = data.drop('pct.1', axis=1)
    data = data.drop('pct.2', axis=1)
    data = data.drop('p_val_adj', axis=1)
    
    data = data.set_index('Gene')
    data = data.loc[:, ~data.columns.str.contains('^Unnamed')]
    print(data)
    aldoc_data = data.loc[data['Aldoc+'] > data['Plcb4+']]
    plcb_data = data.loc[data['Plcb4+'] > data['Aldoc+']]
    aldoc_data = aldoc_data[0:25]
    plcb_data = plcb_data[0:25]
    aldoc_data = aldoc_data.transpose()
    plcb_data = plcb_data.transpose()
    
    #fig, ax = plt.subplots(2,1, figsize=(1.0,3.0))
    fig, ax = plt.subplots(2,1, figsize=(10,1.0))
    sns.heatmap(plcb_data,cmap='Reds', ax=ax[0], vmin=0, annot=False, annot_kws={'color':'black'}, lw=0.5, cbar=True,
               cbar_kws = {'use_gridspec':False, 'location':'left','ticks':[0, 3.75]})
    
    sns.heatmap(aldoc_data,cmap='Blues', ax=ax[1], vmin=0, annot=False, annot_kws={'color':'white'}, lw=0.5, cbar=True,
               cbar_kws = {'use_gridspec':False, 'location':'left','ticks':[0, 3.05]})

    for a in ax.flatten():
        a.tick_params(axis='x', which='both', bottom=False)
        a.tick_params(axis='y', which='both', top=False)
        a.tick_params(axis='y', which='both', right=False)
        a.tick_params(axis='y', which='both', left=False)
        a.tick_params(axis='y', which='both', right=False)
        a.set_xlabel(' ')
        a.set_yticklabels(['Aldoc+','Plcb4+'], rotation = 45)
    ax[0].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax[1].yaxis.set_label_position('left')
    #ax[1].tick_params(axis='y', which='both', right=False)
    ax[1].yaxis.tick_right()
    #ax[0].invert_xaxis()
    #ax[1].invert_xaxis()
    plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=90)
    plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'Fig5b.pdf'), bbox_inches='tight')
    return fig, ax

def AST(data):
    if data['P-value'] < 0.001:
        return '***'
    elif data['P-value'] < 0.01:
        return '**'
    elif data['P-value'] < 0.05:
        return '*'
    return

def fig5e():
    #Data processing
    aldoc_data = pd.read_excel('Fig5e_Aldoc.xlsx')
    plcb4_data = pd.read_excel('Fig5b_Plcb4.xlsx')
    aldoc_data = aldoc_data.loc[aldoc_data['Fold enrichment'] > 10]
    aldoc_data = aldoc_data.loc[aldoc_data['P-value'] < 0.05]
    plcb4_data = plcb4_data.loc[plcb4_data['P-value'] < 0.05]
    aldoc_data['AST'] = aldoc_data.apply(AST, axis=1) 
    plcb4_data['AST'] = plcb4_data.apply(AST, axis=1) 
    
    t1 = aldoc_data.Term.str.split('~', expand=True)
    aldoc_data['Pathways'] = t1[1]
    t2 = plcb4_data.Term.str.split('~', expand=True)
    plcb4_data['Pathways'] = t2[1]
    #aldoc_data = aldoc_data.loc[aldoc_data['Pathways'].str.len() < 30]
    
    fig, ax = plt.subplots(1,2, figsize=(5,2), gridspec_kw={'width_ratios': [5, 1]})
    sns.set_color_codes()
    sns.barplot(x='Pathways',y='Fold enrichment',data=aldoc_data, ax=ax[0], hue='Category', dodge=False, edgecolor='black',
               linewidth=0.5, palette=['g','b','y'])
    sns.barplot(x='Pathways',y='Fold enrichment',data=plcb4_data, ax=ax[1], hue='Category', dodge=False, edgecolor='black',
               linewidth=0.5, palette=['b','y'])
    ax[0].set_xticklabels(aldoc_data.Pathways, rotation=45, ha='right')
    ax[1].set_xticklabels(plcb4_data.Pathways, rotation=45, ha='right')
    ax[0].set_yscale('log') #optional
    ax[1].set_yscale('log') #optional
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_xlabel('')
        a.legend().set_visible(False)
    ax[0].set_ylim(0,1000)
    ax[1].set_ylim(0,100)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'Fig5e.pdf'), bbox_inches='tight')
    return

def fig3b():
    #LTD 'data' generation
    fig, ax = plt.subplots(figsize = (2.5, 2))
    datap = [1.2, 1.4, 1.6]
    for d in datap:
        T1 = np.array([np.arange(0,300,1)])
        no_event_data1 = np.full((1,200), 1)
        event_data = np.full((1,50), d)
        no_event_data2 = np.linspace(1-(d-1), 1, 50)
        no_event_data3 = np.full((0,200),1)
        data = np.append(no_event_data1, event_data)
        data = np.append(data, no_event_data2)
        data = np.append(data, no_event_data3)
        data = pd.DataFrame(data = {'T':T1.flatten(), 'data':data.flatten()})
        sns.regplot(x='T',y='data', data=data, ax=ax, order=15, scatter=False, ci=False, line_kws={'linewidth':0.5})
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #timescale may be entirely incorrect
    ax.set_xlim(0,300)
    ax.set_ylim(0.4,1.8)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Simple spike rate (Normalized)')
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'Fig43.pdf'))
    return fig, ax

def fig3c():
    #LTD 'data' generation
    fig, ax = plt.subplots(figsize = (2.5, 2))
    datap = [[0.7,0.75, '#EC1C24'], [0.85,0.9, '#2E3191'] ,[1.25,1.3, '#EC1C24'], [1.1,1.15, '#2E3191']]
    for d in datap:
        T1 = np.array([np.arange(-20,0,1)])
        T2 = np.array([np.arange(10,45,1)])
        pre_data = np.random.uniform(0.999,1.0,20)
        post_data = np.random.uniform(d[0], d[1] ,35)
        pre = pd.DataFrame(data = {'T':T1.flatten(), 'data':pre_data.flatten()})
        post = pd.DataFrame(data = {'T':T2.flatten(), 'data':post_data.flatten()})
        data = pre.append(post)
        sns.regplot(x='T',y='data', data=data, ax=ax, order=10, scatter=False, 
                    truncate=True, fit_reg=True, ci=False, color=d[2], line_kws={'linewidth':0.5})
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-20,40)
    ax.set_ylim(0.6,1.4)
    ax.set_yticks([0.6,0.8,1.0,1.2,1.4])
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('Synapse strength')
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'Fig3c.pdf'))
    return fig, ax

def fig4b():
    fig, ax = plt.subplots(figsize = (2.5, 2))
    datap = [[0.7,0.75, '#EC1C24'], [1.25,1.3, '#2E3191']]
    for d in datap:
        T1 = np.array([np.arange(-20,0,1)])
        T2 = np.array([np.arange(10,45,1)])
        pre_data = np.random.uniform(0.999,1.0,20)
        post_data = np.random.uniform(d[0], d[1] ,35)
        pre = pd.DataFrame(data = {'T':T1.flatten(), 'data':pre_data.flatten()})
        post = pd.DataFrame(data = {'T':T2.flatten(), 'data':post_data.flatten()})
        data = pre.append(post)
        sns.regplot(x='T',y='data', data=data, ax=ax, order=10, scatter=False, 
                    truncate=True, fit_reg=True, ci=False, line_kws={'linewidth':0.5})
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-20,40)
    ax.set_ylim(0.6,1.4)
    ax.set_yticks([0.6,0.8,1.0,1.2,1.4])
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('Simple spike rate (Normalized)')
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'Fig4b.pdf'))
    return fig, ax

def fig4c():
    fig, ax = plt.subplots(figsize = (2.5, 2))
    neg_data = [[89.9,90, 66,68, '#EC1C24', '--'],[89.9,90, 112,114, '#EC1C24', None]]
    pos_data = [[59.9,60, 48,50, '#2E3191', '--'],[59.9,60, 70,72, '#2E3191', None]]
    for d in neg_data:
        T1 = np.array([np.arange(-20,0,1)])
        T2 = np.array([np.arange(10,45,1)])
        pre_data = np.random.uniform(d[0], d[1], 20)
        post_data = np.random.uniform(d[2], d[3] ,35)
        pre = pd.DataFrame(data = {'T':T1.flatten(), 'data':pre_data.flatten()})
        post = pd.DataFrame(data = {'T':T2.flatten(), 'data':post_data.flatten()})
        data = pre.append(post)
        sns.regplot(x='T',y='data', data=data, ax=ax, order=10, scatter=False, color=d[4], marker = d[5],
                    truncate=True, fit_reg=True, ci=False, line_kws={'linewidth':0.5, 'linestyle':d[5]})
    for d in pos_data:
        T1 = np.array([np.arange(-20,0,1)])
        T2 = np.array([np.arange(10,45,1)])
        pre_data = np.random.uniform(d[0], d[1], 20)
        post_data = np.random.uniform(d[2], d[3] ,35)
        pre = pd.DataFrame(data = {'T':T1.flatten(), 'data':pre_data.flatten()})
        post = pd.DataFrame(data = {'T':T2.flatten(), 'data':post_data.flatten()})
        data = pre.append(post)
        sns.regplot(x='T',y='data', data=data, ax=ax, order=10, scatter=False, color=d[4], 
                    truncate=True, fit_reg=True, ci=False, line_kws={'linewidth':0.5, 'linestyle':d[5]})
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-20,40)
    ax.set_ylim(40,140)
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('Simple spike rate')
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'Fig4c.pdf'))
    return fig, ax

def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header['signal_channels'][0][2])
    data = a.get_analogsignal_chunk() * a.header['signal_channels'][0][5]
    return data[:,0], sr

def fig4d():
    data, sr = read_abf('spikes.ABF')
    plt.plot(data[3*sr+13500:3*sr+21500], color='black', lw=1.5)
    plt.savefig(os.path.join(output_loc, 'Fig4d_l.pdf'), bbox_inches='tight')
    plt.plot(data[3*sr+20200:3*sr+20600], color='black', lw=1.5)
    plt.savefig(os.path.join(output_loc, 'Fig4d_r.pdf'), bbox_inches='tight')
    return

if __name__ == "__main__":
    output_loc = '' #location where image pdfs are stored
    plt.rcParams['font.sans-serif'] = 'Calibri'
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.5
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.5
    color = ['#2E3191','#EC1C24']
    fonts = {"font.size":7, "axes.labelsize":7, "ytick.labelsize":7, "xtick.labelsize":7} 
    plt.rcParams.update(fonts)
