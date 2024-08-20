#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

def density(df, xlabel, ylabel, title, path, **kwargs):
    """
    Creates a simple density plot.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.kdeplot(df, fill=True, ax=ax, **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    fig.suptitle(title, fontsize=65, weight='bold', x=0.36, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def density_x(df_list, xlabel, ylabel, xscale, title, colors, crit_list, path, xvline=0, yminvline=0, ymaxvline=0, **kwargs):
    """
    Creates a density plot with a number x of dfs to represent on the same ax.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    for i, df in enumerate(df_list):
        sns.kdeplot(data=df, fill=True, ax=ax, color=colors[i], **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    ax.set_xscale(xscale)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.vlines(x=xvline, ymin=yminvline, ymax=ymaxvline, 
            linestyles='dashed', colors='black')

    legend_list = []
    for i, crit in enumerate(crit_list):
        legend_element = mpatches.Patch(color=colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1),
                fontsize=20)

    fig.suptitle(title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=500)

def density_x_mod(df_list, xlabel, ylabel, xscale, title, colors, crit_list, path, xvline=0, yminvline=0, ymaxvline=0, **kwargs):
    """
    Creates a density plot with a number x of dfs to represent on the same ax.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    #fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    for i, df in enumerate(df_list):
        sns.kdeplot(df, fill=True, color=colors[i], **kwargs)
    #plt.xlabel(xlabel, fontdict={'fontsize': 35})
    #plt.ylabel(ylabel, fontdict={'fontsize': 35})
    #ax.set_xscale(xscale)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    #plt.spines['right'].set_linewidth(0)
    #plt.spines['top'].set_linewidth(0)
    #ax.vlines(x=xvline, ymin=yminvline, ymax=ymaxvline, 
    #        linestyles='dashed', colors='black')

    legend_list = []
    for i, crit in enumerate(crit_list):
        legend_element = mpatches.Patch(color=colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1),
                fontsize=20)

    plt.title(title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=500)

def density_percentile(df, xlabel, ylabel, title, path, percentile_list, percent_colors, percent_label, **kwargs):
    """
    Creates a density plot with vertical lines 
    correpsonding to percentiles in percentile_list.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.kdeplot(df, fill=True, ax=ax, **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    legend_list = []
    for i, perc in enumerate(percentile_list):
        ax.vlines(x=perc, ymin=0, ymax=0.02, linestyles='dashed', colors=percent_colors[i])
        legend_element = mpatches.Patch(color=percent_colors[i], label=percent_label[i])
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(1,1),
                fontsize=20)
    fig.suptitle(title, fontsize=65, weight='bold', x=0.36, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def lineplot(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path, **kwargs):
    """ 
    Create a vertical connected dot plot or lineplot.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    rc = {'ytick.labelsize': 25, 'xtick.labelsize': 25}
    plt.rcParams.update(**rc)    
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.lineplot(df, x=x_col, y=y_col, hue=hue_col, palette=color_dict, 
                    marker='o', markeredgecolor='grey', **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 30})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 30})
    ax.set_ylim(0, 1.05)
    plt.legend(fontsize=35)
    for handle in plt.gca().get_legend().legendHandles:
        handle.set_linewidth(16)
    
    fig.suptitle(title, fontsize=30, x=0.5, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def lineplot_errorbars(df, x_col, y_col, std_col, hue_col, xlabel, ylabel, title, color_dict, path, **kwargs):
    """ 
    Create a vertical connected dot plot or lineplot with errorbars.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    rc = {'ytick.labelsize': 20, 'xtick.labelsize': 20}
    plt.rcParams.update(**rc)    
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.lineplot(df, x=x_col, y=y_col, hue=hue_col, palette=color_dict, 
                    marker='o', markeredgecolor='grey', **kwargs)
    lower = df[y_col] - df[std_col]
    upper = df[y_col] + df[std_col]                
    ax.plot(df[x_col], lower, color='tab:blue', alpha=0.1)
    ax.plot(df[x_col], upper, color='tab:blue', alpha=0.1)
    ax.fill_between(df[x_col], lower, upper, color='lightblue', alpha=0.2)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 25})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 25})
    ax.set_ylim(0, 1.05)
    fig.suptitle(title, fontsize=30, x=0.5, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

def pie_multiple(y, x, count_list, labels, colors, ax_title, title, legend_title, path, **kwargs):
    """
    Creates 8 pie charts from a list of list (count_list) where each global
    element corresponds to a list of local elements (ex: percent per rank across
    8 model intersection).
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(y, x, figsize=(25, 14))
    plt.subplots_adjust(hspace=0.5)
    ax = axes.flatten()
    for i, element in enumerate(count_list):
        count_per_element = count_list[i][:]
        ax[i].set_title(ax_title[i], fontdict={'fontsize': 25}, x=0.5, y=0.8)
        ax[i].pie(count_per_element, colors=colors, textprops={'fontsize': 19},
                    **kwargs)
        ax[i].axis('equal')
        white_circle = plt.Circle((0, 0), 0.4, color='white') #to create a donut chart
        ax[i].add_artist(white_circle) #to create a donut chart

    fig.suptitle(title, x=0.5, y=0.9, fontsize=25)
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(1.1, 0.5),
                prop={'size': 20}, title=legend_title)
    plt.savefig(path, dpi=600)

def donut(count_list, labels, colors, ax_title, title, legend_title, path, **kwargs):
    """
    Creates a donut chart from a list.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(25, 14))
    plt.subplots_adjust(hspace=0.5)
    ax.set_title(ax_title, fontdict={'fontsize': 25}, x=0.5, y=0.8)
    ax.pie(count_list, colors=colors, textprops={'fontsize': 19},
                **kwargs)
    ax.axis('equal')
    white_circle = plt.Circle((0, 0), 0.4, color='white') #to create a donut chart
    ax.add_artist(white_circle) #to create a donut chart

    fig.suptitle(title, x=0.5, y=0.9, fontsize=25)
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(1.1, 0.5),
                prop={'size': 20}, title=legend_title)
    plt.savefig(path, dpi=600)

def stacked_bar(lists, x_tick_labels, labels, title, xlabel, ylabel, colors, min_y, max_y, optional_annot, path, **kwargs):
    """
    Create a stacked bar chart from a list of lists ('lists').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    df = pd.DataFrame(lists, index=x_tick_labels, columns=labels)
    print(df)
    ax = df.plot.bar(stacked=True, figsize=(30, 12), color=colors, **kwargs)
    ax.set_xticklabels(x_tick_labels, rotation=0)
    
    # Add optional annotation above bars (ex: number of sno in each bar)
    n = len(labels)
    
    def chunker(list_, size):
        # Get chunks of n elements from list_ where n=size
        return (list_[pos:pos + size] for pos in range(0, len(list_), size))
    
    # Get the cumulative height of each bar until the before last stacked bar
    previous_heigths = [0] * len(x_tick_labels)
    for i, chunks_index in enumerate(chunker(list(range(0, n*len(x_tick_labels))[:-len(x_tick_labels)]), len(x_tick_labels))):  # from 0 to the before last bars
        for index_ in chunks_index:
            bar = list(ax.patches)[index_]
            previous_heigths[index_ - len(x_tick_labels) * i] += bar.get_height()
    
    # Add the cumulative height of previous bars to the last stacked bar of each stack
    last_bars = [bar_ for j, bar_ in enumerate(ax.patches) if j in list(range(0, n*len(x_tick_labels))[-len(x_tick_labels):])]  # last (i.e. topmost) bars
    for i, bar in enumerate(last_bars):
        ax.text((bar.get_x() + bar.get_width()/(len(x_tick_labels)/2) - 0.1), 
                (bar.get_height() + previous_heigths[i] + 10), optional_annot[i], fontsize=35)
    plt.legend(fontsize=30, loc='center right', bbox_to_anchor=(1.3, 0.5))


    plt.title(title, fontsize=38)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.autoscale()
    plt.margins(0.02)
    plt.ylim(min_y, max_y)
    plt.savefig(path, bbox_inches='tight', dpi=600)

def stacked_bar2(lists, x_tick_labels, labels, title, xlabel, ylabel, colors, min_y, max_y, optional_annot, path, **kwargs):
    """
    Create a stacked bar chart from a list of lists ('lists').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    df = pd.DataFrame(lists, index=x_tick_labels, columns=labels)
    print(df)
    ax = df.plot.bar(stacked=True, figsize=(15, 12), color=colors, **kwargs)
    ax.set_xticklabels(x_tick_labels, rotation=45)
    
    # Add optional annotation above bars (ex: number of sno in each bar)
    n = len(labels)
    
    def chunker(list_, size):
        # Get chunks of n elements from list_ where n=size
        return (list_[pos:pos + size] for pos in range(0, len(list_), size))
    
    # Get the cumulative height of each bar until the before last stacked bar
    previous_heigths = [0] * len(x_tick_labels)
    for i, chunks_index in enumerate(chunker(list(range(0, n*len(x_tick_labels))[:-len(x_tick_labels)]), len(x_tick_labels))):  # from 0 to the before last bars
        for index_ in chunks_index:
            bar = list(ax.patches)[index_]
            previous_heigths[index_ - len(x_tick_labels) * i] += bar.get_height()
    
    # Add the cumulative height of previous bars to the last stacked bar of each stack
    last_bars = [bar_ for j, bar_ in enumerate(ax.patches) if j in list(range(0, n*len(x_tick_labels))[-len(x_tick_labels):])]  # last (i.e. topmost) bars
    for i, bar in enumerate(last_bars):
        ax.text((bar.get_x() + bar.get_width()/(len(x_tick_labels)/2) - 0.1), 
                (bar.get_height() + previous_heigths[i] + 20), optional_annot[i], fontsize=25)
    plt.legend(fontsize=30, loc='center right', bbox_to_anchor=(1, 0.9))

    #Add horizontal line
    ax.hlines(y=41, xmin=0, xmax=3,
            linestyles='dashed', colors='black')

    plt.title(title, fontsize=38)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.autoscale()
    plt.margins(0.02)
    plt.ylim(min_y, max_y)
    plt.savefig(path, bbox_inches='tight', dpi=600)

def percent_count(count_list):
    """
    Create from a list of lists a percentage list of lists.
    """
    percent_list = []

    for i, cat in enumerate(count_list):
        temp = []
        for j, crit in enumerate(cat):
            total = sum(cat)
            if total != 0:
                percent = crit/total * 100
                percent = round(percent, 2)
            else:
                percent = 0
            temp.append(percent)
        percent_list.append(temp)

    return percent_list
