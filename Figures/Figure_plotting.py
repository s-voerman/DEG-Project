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

def fig1b():
    data = pd.read_excel('genelist1.xlsx')
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
    fig.savefig(os.path.join(output_loc, 'Fig1c.pdf'), bbox_inches='tight')
    return fig, ax

def fig1c():
    data = pd.read_excel('genelist1.xlsx')
    sorted_data = data.sort_values('avg_logFC', ascending=False).reset_index()
    fig, ax = plt.subplots(1,2, figsize=(3.1,1.6), gridspec_kw={'width_ratios': [1, 1]})
    print(len(sorted_data))
    aldoc_data = sorted_data.loc[data.avg_logFC > 0]
    plcb4_data = sorted_data.loc[data.avg_logFC < 0]
    aldoc_data = aldoc_data.reset_index()
    plcb4_data = plcb4_data.reset_index()
    sns.scatterplot(ax=ax[0], x='Aldoc+', y='Plcb4+', data=aldoc_data, color='#2E3191', hue='avg_logFC', palette='Blues', 
                    s=5)
    sns.scatterplot(ax=ax[0], x='Aldoc+', y='Plcb4+', data=plcb4_data, color='#EC1C24', hue='avg_logFC', palette='Reds_r', 
                    s=5) 
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].set_xlim(-0.1,5)
    ax[0].set_ylim(-0.1,5)
    ax[0].set_yticks([0,1,2,3,4,5])
    ax[0].set_xticks([0,1,2,3,4,5])
    ax[0].legend().set_visible(False)
    
    #colorbar
    norm = TwoSlopeNorm(vmin=data['avg_logFC'].min(), vcenter=0, vmax=data['avg_logFC'].max())
    sm = plt.cm.ScalarMappable(cmap='RdBu', norm=norm)
    sm.set_array([])
    fig.colorbar(sm, label='Mean Log Fold Change', ax=ax[0])
    
    #axes
    ax[0].legend().set_visible(False)
    ax[1].axvline(x=25, lw=0.5, color='black', ymin=0.803)
    ax[1].axvline(x=50, lw=0.5, color='black', ymin=0.745)
    ax[1].axvline(x=75, lw=0.5, color='black', ymin=0.71)
    ax[1].axvline(x=450, lw=0.5, color='black', ymax=0.37)
    ax[1].axvline(x=425, lw=0.5, color='black', ymax=0.43)
    ax[1].axvline(x=400, lw=0.5, color='black', ymax=0.457)
    ax[1].axvline(x=246, lw=0.5, color='black', linestyle='--')
    ax[1].axhline(y=0.25, lw=0.5, color='black', linestyle='--')
    ax[1].axhline(y=-0.25, lw=0.5, color='black', linestyle='--')
    
    sns.lineplot(x='level_0', y='avg_logFC', data=aldoc_data, ax=ax[1], lw=0.5, color='#2E3192')
    sns.lineplot(x='level_0', y='avg_logFC', data=plcb4_data, ax=ax[1], lw=0.5, color='#ED1C24')
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].set_xlabel('Genes')
    ax[1].set_ylim(np.min(plcb4_data.avg_logFC),np.max(aldoc_data.avg_logFC))
    ax[1].set_xticks([])
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'Fig1c.pdf'))
    return fig, ax    

def fig1d():
    data = pd.read_excel('genelist1.xlsx')
    data = data.drop('avg_logFC', axis=1)
    data = data.drop('pct.1', axis=1)
    data = data.drop('pct.2', axis=1)
    data = data.drop('p_val_adj', axis=1)
    data = data.set_index('Gene')
    data = data.loc[:, ~data.columns.str.contains('^Unnamed')]
    print(data)
    aldoc_data = data.loc[data['Aldoc+'] > data['Plcb4+']]
    plcb_data = data.loc[data['Plcb4+'] > data['Aldoc+']]
    fig, ax = plt.subplots(1,2, figsize=(1.0,3.0))
    sns.heatmap(plcb_data[0:25],cmap='Reds', ax=ax[0], vmin=0, annot=False, annot_kws={'color':'black'}, cbar=False,
               cbar_kws = {'use_gridspec':False, 'location':'top','ticks':[0, 3.75]})
    sns.heatmap(aldoc_data[0:25],cmap='Blues', ax=ax[1], vmin=0, annot=False, annot_kws={'color':'white'}, cbar=False,
               cbar_kws = {'use_gridspec':False, 'location':'top','ticks':[0, 3.05]})
    for a in ax.flatten():
        a.set_ylabel('')
        a.tick_params(axis='x', which='both', bottom=False)
        a.tick_params(axis='y', which='both', left=False)
        a.tick_params(axis='y', which='both', right=False)
        a.set_xticklabels(['Aldoc+','Plcb4+'], rotation = 45, ha='right')
    ax[1].yaxis.set_label_position('right')
    ax[1].tick_params(axis='y', which='both', right=False)
    ax[1].yaxis.tick_right()
    ax[0].invert_xaxis()
    ax[1].invert_xaxis()
    plt.setp(ax[0].yaxis.get_majorticklabels(), rotation=0)
    plt.setp(ax[1].yaxis.get_majorticklabels(), rotation=0)
    fig.savefig(os.path.join(output_loc, 'Fig1d.pdf'), bbox_inches='tight')
    return fig, ax

def AST(data):
    if data['P-value'] < 0.001:
        return '***'
    elif data['P-value'] < 0.01:
        return '**'
    elif data['P-value'] < 0.05:
        return '*'
    return

def fig2c():
    #Data processing
    aldoc_data = pd.read_excel('EASE_Aldoc.xlsx')
    plcb4_data = pd.read_excel('EASE_Plcb4.xlsx')
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
    fig.savefig(os.path.join(output_loc, 'Fig2c.pdf'), bbox_inches='tight')
    return

def fig4b():
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
    fig.savefig(os.path.join(output_loc, 'Fig4b.pdf'))
    return fig, ax

def fig4c():
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
    fig.savefig(os.path.join(output_loc, 'Fig4c.pdf'))
    return fig, ax

def fig5b():
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
    fig.savefig(os.path.join(output_loc, 'Fig5b.pdf'))
    return fig, ax

def fig5c():
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
    fig.savefig(os.path.join(output_loc, 'Fig5c.pdf'))
    return fig, ax

def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header['signal_channels'][0][2])
    data = a.get_analogsignal_chunk() * a.header['signal_channels'][0][5]
    return data[:,0], sr

def fig5d():
    data, sr = read_abf('spikes.ABF')
    plt.plot(data[3*sr+13500:3*sr+21500], color='black', lw=1.5)
    plt.savefig(os.path.join(output_loc, 'Fig5d_l.pdf'), bbox_inches='tight')
    plt.plot(data[3*sr+20200:3*sr+20600], color='black', lw=1.5)
    plt.savefig(os.path.join(output_loc, 'Fig5d_r.pdf'), bbox_inches='tight')
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