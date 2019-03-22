# -*- coding: utf-8 -*-
"""
Created on Wed May 16 12:05:53 2018

@author: Paulo
"""

import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import csv

path = ' '

#funções

def bandgap(a, b):
    return float(-b)/float(a)

def ph_e(x): #converte comprimento de onda (nm) em energia (eV)
    return 1240/x

def absobance(Trasmitance): #converte transmitância em absorbância
    return 2 - np.log10(Trasmitance)

def abs_coef(Abs): #calcula o coeficiente de aborção
    return (2.303*Abs)/10e-9


def tauc(abs, energy, thickness, r):  #r = 2 for direct allowed transitions;r = 3 for direct forbidden transitions;r = 1/2 for indirect allowed transitions;r = 3/2 for indirect forbidden transitions
     return (((abs*2.303)/thickness)*energy)**(r) #(ahv)^r  #thickness(nm)

def normaliza_0_1(data):
    return data/data.max()

def subzero(data):
    return data - data.min()


def find_files(path):
    files = []
    for i in os.listdir(path):
        if i.endswith(".txt"):
             files.append(path + str(i))
    files.sort()
    return files

def samplename(files, n):
    return files[n].split('/')[-1][:-4]

def create_data_dict(files):
    data={}
    for i in range(0, len(files)):

        data[samplename(files, i)]=pd.read_csv(files[i], sep='\t', header = None, decimal = ',')
        data[samplename(files, i)][3] = ph_e(data[samplename(files, i)][0])
        data[samplename(files, i)][4] = tauc(data[samplename(files, i)][1], data[samplename(files, i)][3], 400, (2))
        data[samplename(files, i)][2] = normaliza_0_1(data[samplename(files, i)][1])
        data[samplename(files, i)][5] = tauc(data[samplename(files, i)][2], data[samplename(files, i)][3], 400, (2))
        data[samplename(files, i)][6] = normaliza_0_1(subzero(data[samplename(files, i)][4]))
        data[samplename(files, i)][7] = normaliza_0_1(subzero(data[samplename(files, i)][1]))
    return data


def dics(data, limit, dif):
    dic = {}
    for z in data:

        a = []
        y=[]
        for i in range(0, len(data[z][6])-1):
            if data[z][6][i] <= 0.013:   #select linear part
                teste=(data[z][6][i] -data[z][6][i+1])
                if teste >= 5e-4:
                    y.append(data[z][6][i])
                    a.append(i)

        x=[]
        for i in range(0, len(data[z][3])):
            for j in a:
             if i == j:
                 x.append(data[z][3][i])

        dic[z] = np.polyfit(x,y,1)
    return dic


def plotes(data, files, xlim, ylim, type, figname, save=False, fit_interval=[1.6,1.8,20]): #type = "tauc" or  "abs"
    def plotpar(xlim, ylim, xlabel='', ylabel=''):
        plt.xlim(xlim) #650,850
        plt.ylim(ylim) #0,.4
        plt.xlabel(xlabel, size=20)
        plt.ylabel(ylabel, size=20)
        plt.tick_params(direction='in', which='both', length=5)
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)

    if type == "abs":
        fig = plt.figure(figsize=(10,8))
        for i in range(0, len(files)):
            plt.plot(data[samplename(files, i)][0], data[samplename(files, i)][7], label = samplename(files, i))
        plotpar(xlim, ylim, xlabel='Wavelength / nm', ylabel='Absorbance / a.u.')
        if save == True:
            plt.savefig(path  + type + "_" + figname + '.eps', format='eps')
        plt.legend(fontsize = 15)

    elif type == "tauc":
        dic = dics(data, 0.013, 5e-4)
        dic_fit ={}
        for z in data:
            dic_fit[z] = np.poly1d(dic[z])

        xfit = np.linspace(fit_interval[0], fit_interval[1], fit_interval[2])

        fig = plt.figure(figsize=(12,9))
        for i in range(0, len(files)):
            plt.plot(data[samplename(files, i)][3], data[samplename(files, i)][6], 'o', label = samplename(files, i))
        for z in range(0,len(dic_fit)):
            plt.plot(xfit,dic_fit[samplename(files, z)](xfit), '--', color='C'+str(z))
        plotpar(xlim, ylim, xlabel='Photon energy, hv / eV', ylabel=r'($\alpha$h$\nu$)$^{2}$ / (eV nm$^{-1}$)$^{2}$')
        if save == True:
            plt.savefig(path  + type + "_" + figname + '.eps', format='eps')
        plt.legend(fontsize = 15)



files = find_files(path)

data=create_data_dict(files)


abs_plot = plotes(data, files, [650,850], [0,.4], "abs", "teste")

tauc_plot = plotes(data, files, [1.55,1.9], [0,0.02], "tauc", "teste")

dic = dics(data, 0.013, 5e-4)

for i in dic:
    print(i + ':  ' + 'y = ' + str(dic[i][0]) + 'x  '  + str(dic[i][1]))


with open(path + 'eqs.csv','w') as data:
    for i in dic:
        data.write(i + ':  ' + 'y = ' + str(dic[i][0]) + 'x  '  + str(dic[i][1]) + '\n\n')
data.close()

bandagaps = {}
for i in dic:
    bandagaps[i]=bandgap(dic[i][0], dic[i][1])

with open(path + 'bandgaps.csv','w') as data:
    for i in dic:
        data.write(i + ':  ' + str(bandgap(dic[i][0], dic[i][1])) + ' eV' + '\n\n')
data.close()
