# %%
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 08:22:16 2018
@author: Jianbo Zhang
"""
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import random
import time
import csv
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
import os
import datetime
import warnings
# from scipy.signal import savgol_filter
import argparse


def smAlleleFreq(popStruc, sizeOfBulk, rep):
    '''
    An AA/Aa/aa individual carries 0%, 50%, and 100% of the alt (a) allele, respectively.
    The AA:Aa:aa ratios are 0.25 : 0.5 : 0.25, 0.5 : 0 : 0.5, and 0.5 : 0.5 : 0, repectively, in a F2 population, in a RIL population,
    and in a back crossed population if A/a does not affect the trait in the population (null hypothesis)
    '''
    freqL = []
    pop = [0.0, 0.5, 1.0]

    if popStruc == 'F2':
        prob = [0.25, 0.5, 0.25]
    elif popStruc == 'RIL':
        prob = [0.5, 0.0, 0.5]
    elif popStruc == 'BC':
        prob = [0.5, 0.5, 0.0]

    for __ in range(rep):
        altFreq = np.random.choice(pop, sizeOfBulk, p=prob).mean()
        freqL.append(altFreq)

    return sum(freqL)/len(freqL)


def stats(inPath, outPath):
    with open(inPath, 'r') as du, open(outPath, 'w', newline='') as outF:
        xie = csv.writer(outF)
        
        next(du)
        fbADI, sbADI = header.index(fbID+'.AD'), header.index(sbID+'.AD')
                
        header.extend([fbID+'.SI', sbID+'.SI', 'Delta.SI', 'DSI_CI', 'fisher_exact'])
        xie.writerow(header)

        for line in du:
            row = line.rstrip('\t\n').split('\t')
 
            try:
                fb_read = [int(row[fbADI].split(',')[0]), int(row[fbADI].split(',')[1])]
                sb_read = [int(row[sbADI].split(',')[0]), int(row[sbADI].split(',')[1])]
                fb_depth, sb_depth = fb_read[0]+fb_read[1], sb_read[0]+sb_read[1]
            except ValueError:          #ValueError: invalid literal for int() with base 10: 'NA'
                fb_read, sb_read = ['NA', 'NA'], ['NA', 'NA']
                fb_depth, sb_depth = 'NA', 'NA'

            try:
                fb_SI = fb_read[1]/fb_depth
                sb_SI = sb_read[1]/sb_depth
                delta_SI = sb_SI - fb_SI
                fe = fisher_exact([fb_read, sb_read])

                smDSI_List = []
                for __ in range(rep):
                    sm_fbAlt = np.random.binomial(fb_depth, fb_Freq)
                    sm_sbAlt = np.random.binomial(sb_depth, sb_Freq)

                    sm_fbSI = sm_fbAlt/fb_depth
                    sm_sbSI = sm_sbAlt/sb_depth
                    smDSI = sm_sbSI - sm_fbSI
                    smDSI_List.append(smDSI)

                ciDSI = np.percentile(smDSI_List, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])

            except (TypeError, ZeroDivisionError):
                fb_SI, sb_SI, delta_SI, ciDSI, fe  = 'NA', 'NA', 'NA', 'NA', 'NA'
                
            row.extend([fb_SI, sb_SI, delta_SI, ciDSI, fe])
            xie.writerow(row)
    print((time.time()-t0)/3600)

    with open(os.path.join(path, 'COMPLETE.txt'), 'w') as xie:
        xie.write('Statistic calculation is complete!')


def zeroSNP(li):
    # Replace 'divide by zero' with the nearnest value. Use 'empty' as a placeholder if it is the first element of the list
    if li != []:
        li.append(li[-1])   # Assign the previous value to the empty sliding window if the list is not empty
    else:
        li.append('empty')  # Assign 'empty' to the first sliding windows that is empty


def replaceZero(li):
    # Replace the 'empty' placeholders at the begining of the list with the nearnest non-empty value
    i = 0
    while li[i]=='empty':
        i += 1

    j = 0
    while j < i:
        li[j] = li[i]
        j += 1


def atoi(text):
    # Source: https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    Source: https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
    alist.sort(key=natural_keys) sorts in human order. http://nedbatchelder.com/blog/200712/human_sorting.html
    '''
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def smThresholds(DF):
    ratioLi = []
    for __ in range(rep):
        smSMPLR = DF.sample(sm_SmplSize, replace=True)
        ldrLi = zip(smSMPLR[fbID+'.LD'], smSMPLR[sbID+'.LD'])
        i = 0
        for ldr in ldrLi:
            fb_smAlt = np.random.binomial(ldr[0], fb_Freq)
            fb_smRead = [ldr[0]-fb_smAlt, fb_smAlt]
            sb_smAlt = np.random.binomial(ldr[1], sb_Freq)
            sb_smRead = [ldr[1]-sb_smAlt, sb_smAlt]
            
            __, fe_pvalue = fisher_exact([fb_smRead, sb_smRead])
            if fe_pvalue < smAlpha:
                i += 1

        ratioLi.append(i/sm_SmplSize)

    return np.percentile(ratioLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])


def bsaseqPlot(chrmIDL, datafr, datafrT):
    '''
    chrEnd: a list containing the endpoint of each chromosome
    wmL: list of warning message
    numOfSNPInSwL: a list containing the number of SNPs in each sliding window
    points: a dictionary with the chromosome ID as its keys; the value of each key is a list containing
            the chromosome ID, the ltaSNPs/totalSNPs in each sliding window, and the midpoint of the sliding window
    '''
    chrEnd, misc_info, wmL = [], [], []
    numOfSNPInSwL, numOfSNPOnChr, points = [], [], {}

    # if not os.path.exists(fltr):
    #     os.mkdir(fltr)

    # Analyze each chromsome separately
    i = 1
    for chrmID in chrmIDL:
        ch = datafr[datafr.CHROM==chrmID]
        chT = datafrT[datafrT.CHROM==chrmID]
        numOfSNPOnChr.append([len(ch.index), len(chT.index), len(ch.index)/len(chT.index)])

        # outFile = fltr + r'/' + chrmID + '.csv'
        # outChrLtaSNP = os.path.join(fltr, 'Chr' + chrmID + '.csv')
        # ch.to_csv(outChrLtaSNP, sep=',', encoding='utf-8', index=False)

        # Sliding window. swStr: the begining of the window; swEnd: the end of the window; icrs: incremental step
        swStr, swEnd, icrs = 1, 2000000, 10000
        # plotSP = swEnd/2
        # x, y, y5, y6 y8, y9 are lists, each sliding window represents a single data point
        x, y, y6, y8, y9, yT, yRatio = [], [], [], [], [], [], []
        while swEnd <= ch['POS'].max()+1:
            # A single sliding window - dataframe
            # swDF: ltaSNPs in a sliding window; swDFT: all SNPs in a sliding window
            swDF = ch[(ch.POS>=swStr) & (ch.POS<swEnd)]
            swDFT = chT[(chT.POS>=swStr) & (chT.POS<swEnd)]

            rowInSwDF = len(swDF.index)                 # number of ltaSNPs in a sliding window
            rowInSwDFT = len(swDFT.index)               # number of total SNPs in a sliding window

            x.append((swStr+swEnd)/2)                   # Append the midpoint of a sliding window to x
            y.append(rowInSwDF)                         # Append number of ltaSNPs in a sliding window to y
            yT.append(rowInSwDFT)                       # Append number of totalSNPs in a sliding window to yT

            # len(swDL) or len(swDLT) would be zero if no SNP in a sliding window
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    # y5.append(swDFT['G_S'].sum()/rowInSwDFT)
                    y6.append(swDFT['Delta.SI'].sum()/rowInSwDFT)
                    # y7.append(swDFT['GS_CI0995'].sum()/rowInSwDFT)
                    y8.append(swDFT['DSI_CI0005'].sum()/rowInSwDFT)
                    y9.append(swDFT['DSI_CI0995'].sum()/rowInSwDFT)
                except Warning as w:
                    wmL.append(['No ltaSNP', i, int((swStr+swEnd)/2), w])
                    # zeroSNP(y5)
                    zeroSNP(y6)
                    # zeroSNP(y7)
                    zeroSNP(y8)
                    zeroSNP(y9)

            try:
                yRatio.append(float(rowInSwDF/rowInSwDFT))  # Append the ratio of ltaSNP/totalSNP in a sliding window to yRatio
            except ZeroDivisionError as e:
                wmL.append(['No SNP', i, int((swStr+swEnd)/2), e])
                zeroSNP(yRatio)

            numOfSNPInSwL.append(rowInSwDF)

            if i not in points:
                points[i] = [['Chr'+str(i), yRatio[-1], int((swStr+swEnd)/2)]]
            else:
                points[i].append(['Chr'+str(i), yRatio[-1], int((swStr+swEnd)/2)])

            swStr += icrs
            swEnd += icrs

        if y6[-1]=='empty':
            print(f'No ltaSNP on Chr{i}')
            break

        # Replace the 'empty' values at the begining of the lists with nearest non-empty value
        yLists = [y6, y8, y9, yRatio]
        for yl in yLists:
            if 'empty' in yl:
                replaceZero(yl)

        pIndex = 0
        while points[i][pIndex][1] == 'empty':
            pIndex += 1

        j = 0
        while j < pIndex:
            points[i][j][1] = points[i][pIndex][1]
            j += 1

        # ax1[0,i-1].plot(x, y, c='k')
        # ax1[0,i-1].plot(x, yT, c='gray')

        if ldIndex == 0:
            axs[ldIndex,i-1].set_title('Chr'+str(i))
        # axs[0,i-1].plot(x, y5, c='k')
        # axs[0,i-1].plot(x, y7, c='r')

        # axs[2,i-1].set_xticks(np.arange(0, max(x), 10000000))
        # ticks = axs[2,i-1].get_xticks()*1e-6
        # axs[2,i-1].set_xticklabels(ticks.astype(int))

        axs[ldIndex,i-1].plot(x, y6, c='k')
        axs[ldIndex,i-1].plot(x, y8, c='r')
        axs[ldIndex,i-1].plot(x, y9, c='r')

        # axs[2,i-1].plot(x, yRatio, c='k')
        # sg_yRatio = savgol_filter(yRatio, 11, 3)
        # axs[2,i-1].plot(x, sg_yRatio, c='r')

        if ldIndex == len(ldUL)-1:
            axs[ldIndex,i-1].set_xticks(np.arange(0, max(x), 10000000))
            ticks = axs[ldIndex,i-1].get_xticks()*1e-6
            axs[ldIndex,i-1].set_xticklabels(ticks.astype(int))

        chrEnd.append((swStr+swEnd)/2)
        i += 1
    print(f'Preparation for plotting - complete - {time.time()-t0} seconds')

    wrnLog = os.path.join(results, 'wrnLog.csv')
    with open(wrnLog, 'w', newline='') as outF1:
        xie1 = csv.writer(outF1)
        xie1.writerow(['Type', 'Chr', 'Position', 'Warning Message'])
        xie1.writerows(wmL)

    numOfSNPOnChrFile = os.path.join(results, 'numOfSNPOnChrFile.csv')
    with open(numOfSNPOnChrFile, 'w', newline='') as outF2:
        xie2 = csv.writer(outF2)
        xie2.writerow(['Num of ltaSNPs', 'Num of totalSNPs', r'ltaSNPs/totalSNPs'])
        xie2.writerows(numOfSNPOnChr)
    return misc_info


t0 = time.time()
plt.rc('font', family='Arial', size=22)     # controls default text sizes
plt.rc('axes', titlesize=22)                # fontsize of the axes title
plt.rc('axes', labelsize=22)                # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)               # fontsize of the tick labels
plt.rc('ytick', labelsize=20)               # fontsize of the tick labels
plt.rc('legend', fontsize=20)               # legend fontsize
plt.rc('figure', titlesize=22)              # fontsize of the figure title
# plt.tick_params(labelsize=20)

# construct the argument parser and parse the arguments
# python PyBSASeq_MS.py -i snpPE_final.tsv -o BSASeqPE.csv -f 430 -s 385 -p F2
ap = argparse.ArgumentParser()
ap.add_argument('-i', '--input', required=False, help='file name of the GATK4-generated tsv file', default='snp_final.tsv')
ap.add_argument('-o', '--output', required=False, help='file name of the output csv file', default='BSASeq.csv')
ap.add_argument('-f', '--fbsize', type=int, required=False, help='number of individuals in the first bulk', default=430)
ap.add_argument('-s', '--sbsize', type=int, required=False, help='number of individuals in the second bulk', default=385)
ap.add_argument('-p', '--popstrct', required=False, choices=['F2','RIL','BC'], help='population structure', default='F2')
ap.add_argument('--alpha', type=float, required=False, help='p-value for fisher\'s exact test', default=0.01)
ap.add_argument('--smalpha', type=float, required=False, help='p-value for calculating threshold', default=0.1)
ap.add_argument('-r', '--replication', type=int, required=False, help='the number of replications for threshold calculation', default=10000)
ap.add_argument('--smplsize', type=int, required=False, help='simulation sample size for threshold calculation', default=300)
args = vars(ap.parse_args())

popStr = args['popstrct']
rep = args['replication']
fb_Size, sb_Size = args['fbsize'], args['sbsize']
sm_SmplSize = args['smplsize']
alpha, smAlpha = args['alpha'], args['smalpha']

# popStr = 'F2'
# rep = 10000
# fb_Size, sb_Size = 430, 385
# sm_SmplSize = 300

fb_Freq = smAlleleFreq(popStr, fb_Size, rep)
sb_Freq = smAlleleFreq(popStr, sb_Size, rep)

alpha, smAlpha, = 0.01, 0.1

path = os.getcwd()
currentDT = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
inFile, oiFile = os.path.join(path, args['input']), os.path.join(path,'snp_si.csv')
misc = []

with open(inFile, 'r') as du1:   
    header = next(du1).rstrip('\t\n').split('\t')
    bulks = []

    for ftrName in header:
        if ftrName.endswith('.AD'):
            bulks.append(ftrName.split('.')[0])

    fbID, sbID = bulks[0], bulks[1]

results = os.path.join(path, 'Results', currentDT)
if not os.path.exists(results):
    os.makedirs(results)

if os.path.isfile(os.path.join(path, 'COMPLETE.txt')) == False:
    stats(inFile, oiFile)

snpDF = pd.read_csv(oiFile, dtype={'CHROM': str})
columnL = snpDF.columns.values.tolist()
misc.append(['Header', columnL])

misc.extend([['Bulk ID', bulks], ['Shape - entire dataframe', snpDF.shape]])

snpDF.dropna(inplace=True)
misc.append(['Shape - dataframe without NA', snpDF.shape])

# Create and sort chromosome ID list
chrmL = list(set(snpDF['CHROM'].tolist()))
# Remove unordered, mt, and cp 'chromosomes' from the list
chrmIDL = []
for eml in chrmL:
    if eml.isdigit():
        chrmIDL.append(eml)

chrmIDL.sort(key=natural_keys)

numOfChrs = len(chrmIDL)
misc.extend([['Chromosome ID', chrmL], ['Cleaned chromosome ID',chrmIDL]])

# Calculate the total read depth of each SNP (Alt+Ref, or locus depth) in each bulk
fbAD = snpDF[fbID+'.AD'].str.split(',', expand=True)
snpDF[fbID+'.LD'] = fbAD[0].astype(int) + fbAD[1].astype(int)
sbAD = snpDF[sbID+'.AD'].str.split(',', expand=True)
snpDF[sbID+'.LD'] = sbAD[0].astype(int) + sbAD[1].astype(int)

snpDF['seqCvrg'] = snpDF[fbID+'.LD'] + snpDF[sbID+'.LD']

maxLD = max([snpDF[fbID+'.LD'].max(), snpDF[sbID+'.LD'].max()])
ldUL = [maxLD, 80, 60, 40]
# ldUL = [maxLD]

# Create a list of chromosome sizes
chrmSzL = []
for chrmID in chrmIDL:
    chrmSzL.append(snpDF[snpDF.CHROM==chrmID]['POS'].max())

misc.append(['Chromosome sizes', chrmSzL])

fig, axs = plt.subplots(nrows=4, ncols=numOfChrs, figsize=(20, 12), sharex='col', sharey='row', \
            gridspec_kw={'width_ratios': chrmSzL, 'height_ratios': [0.6,0.6,0.75,1]})

# Separate the 99%, 95%, and 90% confidence intervals and add them to the dataframe
snpDF['DSI_CI'] = snpDF['DSI_CI'].apply(lambda x: x[1:-1])
splitDSI_CI = snpDF['DSI_CI'].str.split(expand=True).astype(float)
snpDF['DSI_CI0005'] = splitDSI_CI[0]
snpDF['DSI_CI0995'] = splitDSI_CI[1]

snpDF['fisher_exact'] = snpDF['fisher_exact'].apply(lambda x: x[1:-1])
splitFE = snpDF['fisher_exact'].str.split(',', expand=True).astype(float)
snpDF['FE_P'] = splitFE[1]

ldIndex = 0
for ld in ldUL:
    miscLD = []
    miscLD.extend(misc)
    ldDir = os.path.join(results, str(ld))
    if not os.path.exists(ldDir):
        os.makedirs(ldDir)

    # fltrL001 = os.path.join(ldDir, 'l0' + str(alpha).split('.')[1])
    # fltrGE001 = os.path.join(ldDir, 'ge0' + str(alpha).split('.')[1])
    # nfltr = os.path.join(ldDir, 'total')
    # Filter out the entry in 'ALT' column with more than one base
    qualityDF = snpDF[(snpDF[fbID+'.GQ']>=20) & (snpDF[sbID+'.GQ']>=20) & (snpDF['ALT'].str.len()==1) \
            & (snpDF[fbID+'.LD']<=ld) & (snpDF[sbID+'.LD']<=ld)]
    miscLD.append(['Dataframe filtered with quality scores', qualityDF.shape])

    # thrshld = smThresholds(qualityDF)[0][1]

    # Filter out the entries with similar 'AD' ratio in both bulks
    fe = qualityDF[qualityDF['FE_P']<alpha]
    # nfe = qualityDF[qualityDF['G_P']>=alpha]
    #print(np.percentile(nfe['G_S'], [0.5,99.5,2.5,97.5]))

    miscLD.append([f'Average locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].mean(), fe[sbID+'.LD'].mean(), ld]])
    miscLD.append([f'Maximum locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].max(), fe[sbID+'.LD'].max(), ld]])
    miscLD.append([f'Minimum locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].min(), fe[sbID+'.LD'].min(), ld]])
    #miscLD.append([f'Mode locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].mode(), fe[sbID+'.LD'].mode(), ld]])
    miscLD.append([f'Median locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].median(), fe[sbID+'.LD'].median(), ld]])
    miscLD.append([f'DSI_CI0005 and DSI_CI0995', [fe.DSI_CI0005.mean(), fe.DSI_CI0995.mean(), ld]])
    miscLD.append([f'Dataframe filtered with Fisher-Exact - p<{alpha}', fe.shape])
    # miscLD.append([f'Dataframe filtered with Fisher-Exact - p>={alpha}', nfe.shape])
    print(f'Data manipulation - complete - {time.time()-t0} seconds')

    msc_i = bsaseqPlot(chrmIDL, fe, qualityDF)
    miscLD.extend(msc_i)

    miscLD.append(['Running time', [(time.time()-t0)/3600]])

    with open(os.path.join(ldDir, 'misc_info.csv'), 'w', newline='') as outF:
        xie = csv.writer(outF)
        xie.writerows(miscLD)
    
    ldIndex += 1

# miscLD = []
# miscLD.extend(misc)
# if not os.path.exists(results):
#     os.makedirs(results)

# # fltrL001 = os.path.join(results, 'l0' + str(alpha).split('.')[1])

# # Filter out the entry in 'ALT' column with more than one base
# # qualityDF = snpDF[(snpDF[fbID+'.GQ']>=20) & (snpDF[sbID+'.GQ']>=20) & (snpDF['ALT'].str.len()==1)]
# qualityDF = snpDF[(snpDF[fbID+'.GQ']>=99) & (snpDF[sbID+'.GQ']>=99) & (snpDF['ALT'].str.len()==1) \
#             & (snpDF['seqCvrg']<=400) & (snpDF['seqCvrg']>=100) & (snpDF[fbID+'.LD']>=40) & (snpDF[sbID+'.LD']>=40)]
# miscLD.append(['Dataframe filtered with quality scores', qualityDF.shape])

# thrshld = 0.12666666666666668
# # thrshld = smThresholds(qualityDF)[1]
# print(thrshld)

# # Filter out the entries with similar 'AD' ratio in both bulks
# fe = qualityDF[qualityDF['FE_P']<alpha]
# nfe = qualityDF[qualityDF['FE_P']>=alpha]
# # print(np.percentile(nfe['G_S'], [0.5,99.5,2.5,97.5]))

# miscLD.append([f'Average locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].mean(), fe[sbID+'.LD'].mean()]])
# miscLD.append([f'Maximum locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].max(), fe[sbID+'.LD'].max()]])
# miscLD.append([f'Minimum locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].min(), fe[sbID+'.LD'].min()]])
# # miscLD.append([f'Mode locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].mode(), fe[sbID+'.LD'].mode()]])
# miscLD.append([f'Median locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].median(), fe[sbID+'.LD'].median()]])
# miscLD.append([f'DSI_CI0005 and DSI_CI0995', [fe.DSI_CI0005.mean(), fe.DSI_CI0995.mean()]])
# miscLD.append([f'Dataframe filtered with Fisher-Exact - p<{alpha}', fe.shape])
# miscLD.append([f'Dataframe filtered with Fisher-Exact - p>={alpha}', nfe.shape])
# print(f'Data manipulation - complete - {time.time()-t0} seconds')


# fig.tight_layout(pad=0.2, rect=[0.002, 0.028, 1, 1])
fig.subplots_adjust(top=0.971, bottom=0.062, left=0.048, right=0.995, hspace=0.052, wspace=0.092)
fig.suptitle('Genomic position (Mb)', y=0.002, ha='center', va='bottom')
fig.text(0.001, 0.997, 'A', weight='bold', ha='left', va='top')
fig.text(0.001, 0.76, 'B', weight='bold', ha='left', va='bottom')
fig.text(0.001, 0.57, 'C', weight='bold', ha='left', va='bottom')
fig.text(0.001, 0.338, 'D', weight='bold', ha='left', va='bottom')
# fig.text(0.0, 0.115, 'E', weight="bold", ha='left', va='bottom')

fig.savefig(os.path.join(results, 'DSI.pdf'))
fig.savefig(os.path.join(results, 'DSI.png'), dpi=600)
