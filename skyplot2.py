# -*- coding: utf-8 -*-

"""
Created on Tue Mar  9 22:14:14 2021

@author: Maciek
updated by : Maria
"""
import numpy as np
from pylab import *
import matplotlib.pyplot as plt 
from matplotlib.pyplot import rc, rcParams, grid 
import matplotlib.patches as mpatches
import matplotlib.animation as animation

###wyglada lepiej niz to sie probuje robic xd

def plot_skyplot(sat_positions):
    # sat_positions - [PRN, el, az] w stopniach
    rc('grid', color='gray', linewidth=1, linestyle='--')
    fontsize = 20
    rc('xtick', labelsize = fontsize)
    rc('ytick', labelsize = fontsize)
    rc('font', size = fontsize)
    # define colors
    green   ='#467821'
    blue    ='#348ABD'
    red     ='#A60628'
    orange  ='#E24A33'
    purple  ='#7A68A6'
    # start ploting
    fig = plt.figure(figsize=(8, 6))
    fig.subplots_adjust(bottom=0.1,
                        top=0.85,
                        left=0.1,
                        right=0.74)
    ax = fig.add_subplot(polar=True)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticks(range(0, 90 + 10, 10))
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    #definiuję wykres
    scat = ax.scatter([], [], lw=2)
    #slownik annotacji, które potem są zmieniane ich pozycje
    anotations ={}
    #storzenie annotacji i zapis do słownika
    for (PRN, el, az) in sat_positions:
        a = len(el)
        anotations[PRN] = ax.annotate(PRN,
                    xy=(0, -70),                        #stworzenie gdzieś gdzie ich nie widać
                    bbox=dict(boxstyle="round", fc=green, alpha=0.5),
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    color='k')
    PG = 0
    gps = mpatches.Patch(color=green, label='{:02.0f}  GPS'.format(PG))
    legend = ax.legend(handles=[gps], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    def animate(i):
        cord = []
        PG = 0  # zliczanie widocznych satelitów GPS
        for (PRN, el, az) in sat_positions:
            #zmienianie pozycji satelitów i  ich podpisów
            cord.append([np.deg2rad(az[i]), 90 - el[i]])
            anotations[PRN].set_position((np.deg2rad(az[i]), 90-el[i]))
            anotations[PRN].xy = (np.deg2rad(az[i]), 90-el[i])
            if el[i] > 0:
                PG+=1
        htext = legend.get_texts()[0]
        htext.set_text('{:02.0f}  GPS'.format(PG))      #wyświetlanie widocznej liczby satelitów
        scat.set_offsets(cord)          #wyświetlanie punktów w ktorych znajduja sie satelity
    #wyywołanie animacji
    anim = animation.FuncAnimation(fig, animate, frames=a-1, interval=300)     # zmiana szybkości - interval - w milisekundach
    #animacja przedstawia jak zmienia się polozenie satelitow w zaleznosci od czasu -
    #dziala dla takiego okresu dla jakiego sa podane dane wejsciowe do funkji
    #anim.save('satelits.mp4')      #zapisywanie animacji
    plt.show()


if __name__ == '__main__':
    # sat_positions = np.ones((10,20,3))
    sat_positions = [['G04', np.hstack((np.nan,np.arange(30,90,step=0.1))) ,np.hstack((np.nan,np.arange(90,150,step=0.1)))],['G08', [55,56,57], [45,55,65]],['G09', [10], [270]]]
    plot_skyplot(sat_positions)
    
    