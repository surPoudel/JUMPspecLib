import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from logFunctions import *

# # Add a Target-Decoy Column

def countcumulativeTarget(row):
    if row["Type"] == "Target":  #Define Target 
        value = 1
    else:
        value =0
    return value

def countcumulativeDecoy(row):
    if row["Type"] == "Decoy":  #Define Decoy
        value = 1
    else:
        value =0
    return value



def calcFDR(row):
    if row.cumsumTarget == 0:
        FDR = 100
    else:
        FDR = row.cumsumDecoy/row.cumsumTarget*100
    return FDR



def FDR_Target_Decoy(df,sortCol="JDscore"):
    reqd_cols = list(df.columns)
    df1 = df.copy()
    #define Target Decoy
    # df1["Target-Decoy"]=df1.Peptide_ID.apply(lambda row: "Decoy" if "Decoy" in row else "Target")
    # write_log("\nSorting matrix using {} and computing FDR\n".format(sortCol))
    if sortCol == "JDscore":
        df2=df1.sort_values([sortCol,"Type"],ascending=[False,False])   #for labeled dataset
    else:
        df2=df1.sort_values([sortCol,"Type"],ascending=[True,False])   #for labeled dataset
    df2["Target_value"] = df2.apply(countcumulativeTarget, axis=1)
    df2["Decoy_value"] = df2.apply(countcumulativeDecoy, axis=1)
    #df['SUM_C'].cumsum()
    df2["cumsumTarget"] = df2["Target_value"].cumsum()
    df2["cumsumDecoy"] = df2["Decoy_value"].cumsum()
    df2["FDR"] = df2.apply(calcFDR, axis=1)
    if "FDR" in reqd_cols:
        addedColsAll = reqd_cols
    else:
        addedColsAll = reqd_cols+["FDR"]
    return df2[addedColsAll]

def histogramPlot(matched_df, unmatched_df, xaxis, figname,label1, label2): #bins2 for xaxis label
    minv = np.min(unmatched_df[xaxis])
    maxv = np.max(matched_df[xaxis])
    bins = np.linspace(minv,maxv)

    plt.rcParams.update({'font.size': 10})
    fig,ax = plt.subplots(figsize=(4,2))
    plt.yticks(color="black")
    # size, scale = 1000, 10

    commutes2 = matched_df[xaxis]
    commutes2.plot.hist(grid=False, bins=bins, rwidth=0.9,
                   color='#F4F6F7',edgecolor='black', linewidth=1.0)
    commutes = unmatched_df[xaxis]
    commutes.plot.hist(grid=False, bins=bins, rwidth=0.9,
                   color='#808B96',edgecolor='black', linewidth=1.0)
    # bins_labels(bins2)
    # Hide grid lines
    # plt.grid(False)

    plt.title('')
    plt.xlabel(xaxis)
    plt.ylabel('Number of PSMs')
    #   plt.xticks(bins)
    # plt.grid(axis='y', alpha=0.75)
    plt.legend([label1, label2],loc="best")
    figurename = figname+".pdf"
    figurename1 = figname+".png"
    fig.savefig(figurename, bbox_inches="tight", dpi=600 )
    fig.savefig(figurename1, bbox_inches="tight", dpi=600 )


def scatterPlot(df, xaxis, yaxis, hueCol,xlabel,ylabel, figname):
    fig, ax = plt.subplots(figsize=(4,2))
    plt.style.use("ggplot")


    plt.rcParams['axes.edgecolor'] = "#010101"
    plt.rcParams['axes.facecolor'] = '#FFFFFF'


    plt.rcParams.update({'font.size': 6,'figure.max_open_warning': 0})
    ax = plt.axes(facecolor="w")


    sns.scatterplot(data=df, x=xaxis, y=yaxis,hue=hueCol,s =0.5, palette="deep", edgecolor = 'none')

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # plt.hlines(0, 50, 2200, color = "green", linestyles="dotted",linewidth=0.5)
    ax.set_xlabel(xlabel, color="black")
    ax.set_ylabel(ylabel, color="black")
    # ax.set_xlim(0,max(match_list)+100)
    ax.tick_params(axis='x', colors='black')
    ax.tick_params(axis='y', colors='black')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


    figurename = figname

    plt.savefig(figurename, bbox_inches="tight", dpi=600 )
