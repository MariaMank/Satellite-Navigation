from read_yuma import read_yuma
from date2tow import date2tow
from skyplot2 import plot_skyplot
from groundtrack_stud import latlon, groundtrack
import numpy as np
import matplotlib.pyplot as plt

yuma = read_yuma("almanac.txt")

### DANE DO WPISANIA ###
######DATA POCZĄTKOWA#########
data1 = [2022, 2, 22, 0, 0, 0]
########DATA KONCOWA###########
data2 = [2022, 2, 23, 0, 0, 0]
###SKOK - W JAKIM ODSTĘPIE CZASOWYM MAJĄ SIĘ PRZELICZAC WYNIKI#####
skok = 900 #sekundy co ile mają się zliczać wyniki w s - tutaj 900 s = 15 min
#####WYBÓR MIEJSCA DLA KTÓREGO LICZYMY WYNIKI #####
fi = 87 #stopnie
lb = 9
h = 100 #metry
########MASKA OBSERWACJI#####
maska = 10
######TABLICA NUMEROW INTERESUJACYCH NAS SATELITOW####
satelity = []


####OBLICZENIA###
week1, tow1 = date2tow(data1)
week2, tow2 = date2tow(data2)

a = 6378137 #metry
e2 = 0.00669438002290
fi = np.deg2rad(fi)
lb = np.deg2rad(lb)
R = np.array([[-np.sin(fi) * np.cos(lb), -np.sin(lb), np.cos(fi) * np.cos(lb)],
              [-np.sin(fi) * np.sin(lb), np.cos(lb), np.cos(fi) * np.sin(lb)],
              [np.cos(fi), 0, np.sin(fi)]])
result = {}               #wyniki azymuity etc poszczegolnych satelit,

#### funkcja zamieniające wspolrzedne geodezyjne na xyz
def geo2xyz(fi, lam, h, a, e2):
    N = a / (np.sqrt(1 - e2*np.sin(fi)**2))
    x = (N+h)*np.cos(fi)*np.cos(lam)
    y = (N+h)*np.cos(fi)*np.sin(lam)
    z = (N*(1-e2)+h)*np.sin(fi)
    return x, y, z

xr = geo2xyz(fi, lb, h, a , e2)         #xyz odbiornika


def toxyz (nav, week, tow):
    e = nav[2]
    toa = nav[3]
    i = nav[4]
    omega = nav[5]
    sa = nav[6]
    a = sa ** 2
    omega0 = nav[7]
    w = nav[8]
    M0 = nav[9]
    gpsweek = nav[12]
    q = 3.986005 * 10 ** (14)
    we = 7.2921151467 * 10 ** (-5)
    t = week * 7 * 86400 + tow
    toa2 = gpsweek * 7 * 86400 + toa
    tk = t - toa2
    n = np.sqrt(q/a**3)
    Mk = M0 + n*tk
    Ei = Mk
    while True:
        Ei1 = Mk + e*np.sin(Ei)
        if np.abs(Ei - Ei1) < 10 ** (-12):
            break
        Ei = Ei1
    Ek = Ei1
    vk = np.arctan2(np.sqrt(1 - e**2)*np.sin(Ek), np.cos(Ek) - e)
    fik = vk + w
    rk = a*(1 - e*np.cos(Ek))
    xk = rk*np.cos(fik)
    yk = rk*np.sin(fik)
    omegak = omega0 + (omega - we)*tk - we*toa
    Xk = xk*np.cos(omegak) - yk*np.cos(i)*np.sin(omegak)
    Yk = xk * np.sin(omegak) + yk * np.cos(i) * np.cos(omegak)
    Zk = yk*np.sin(i)
    return Xk, Yk, Zk

####funkcja odpowiedzialna za wyliczenia dla wszystkich satelitow
def wyliczenia(tow1, tow2, skok):
    parametryKon = np.zeros((0, 7))
    nry = []                        #tablica numerow wszystkich wystepujacych satelitow - kluczy slownika
    epoki = []                      #tablica godzin od poczatku obserwacji w ktorych zostaly zapisane dane
    i = 0
    tow0 = tow1
    tow3 = tow2
    #gdy zmienia się tydzień satelity zerują się sekundy i generuje to błędy
    #w zwiazku z tym
    if week1 != week2:
        i = week2 - week1
        tow2 = 604800
    while i>=0:
        for t in range(tow1, tow2, skok):       #co 15 min - 900 s
            A = np.zeros((0, 4))
            for nav in yuma:
                nr = nav[0]
                tow = t
                xs = toxyz(nav, week1, tow) #xyz satelity
                xsr = np.array([xs[0] - xr[0], xs[1] - xr[1], xs[2] - xr[2]]) #wektor od odbiornika do satelity geocentrycznie
                r = np.sqrt(xsr[0]**2 + xsr[1]**2 + xsr[2]**2)
                neu = R.transpose()@xsr.transpose()     #mnozenie macierzowe zamiast dot, R.T też powinno transponować
                #neu = np.dot(R.transpose(), xsr.transpose())
                Az = np.rad2deg(np.arctan2(neu[1],neu[0]))
                s = np.sqrt(neu[0] ** 2 + neu[1] ** 2 + neu[2] ** 2)
                el = np.rad2deg(np.arcsin(neu[2] / s))
                #zapis danych do słownika result
                try:
                    S = np.array(result[nr])
                except KeyError:
                    S = np.zeros((0, 9))
                    S = np.append(S, np.array([[t, Az, el, neu[0], neu[1], neu[2], xs[0], xs[1], xs[2]]]), axis=0)
                    nry.append(nr)
                else:
                    S = np.append(S, np.array([[t, Az, el, neu[0], neu[1], neu[2], xs[0], xs[1], xs[2]]]), axis=0)
                result[nr] = S
                if el > maska: #el tylko w zakresie od -90 do 90
                    A1 = np.array([[-(xs[0] - xr[0])/r, -(xs[1] - xr[1])/r, -(xs[2] - xr[2])/r, 1 ]])   #jeden wiersz macierzyA
                    A = np.append(A, A1, axis=0)
            Q = np.linalg.inv(A.T @ A)  # odwrotnoosc macierzy
            GDOP = np.sqrt(Q[0][0] + Q[1][1] + Q[2][2] + Q[3][3])
            Q1 = Q[:3, :3]
            Qneu = R.T @ Q1 @ R
            HDOP = np.sqrt(Qneu[0][0]+Qneu[1][1])
            VDOP = np.sqrt(Qneu[2][2])
            PDOP = np.sqrt(Q[0][0] + Q[1][1] + Q[2][2])
            TDOP = np.sqrt(Q[3][3])
            #PDOP1 = np.sqrt(HDOP**2 + VDOP**2)
            iloscSat = A.shape[0]       #liczba satelitow spelniajacych warunki - liczba wierszy macierzy A
            parametryKon = np.append(parametryKon, np.array([[t, GDOP, HDOP, VDOP, PDOP, TDOP,iloscSat]]), axis=0 )
            if i == (week2 - week1):
                epoki.append((t - tow1) / 3600)
            elif i == 0 and (week2 - week1) != 0:
                epoki.append((604800 - tow0 + (week2 - week1 -1)*604800 + t)/3600)
            else:
                epoki.append((604800 - tow0 + (week2 - week1 - i - 1) * 604800 + t) / 3600)
        tow1 = 0
        if i == 1:
            tow2 = tow3
        i-= 1
    return parametryKon, result, nry, epoki         #wynikiem funkcji jest tablica wspolczynnikow DOP dla każdej epoki - parametryKon,
    # słownik result - w ktorym sa wyniki z wszystkich epok dla kazdego satelity

def rysowanie():
    sat_positioins = []
    parametryKon, result, nry, epoki = wyliczenia(tow1, tow2, skok)
    liczbaSat = parametryKon[:, 6]
    GDOPY = parametryKon[:, 1]      #wyciecie tablic w ktorych ssa wyniki poszczegolnych wspolczynnikow dla kazdej epoki, azeby łatwo przedstawić je
    HDOPY = parametryKon[:, 2]
    VDOPY = parametryKon[:, 3]
    PDOPY = parametryKon[:, 4]
    TDOPY = parametryKon[:, 5]
    sat_pos_latlon = {}
    ilosc_godz = (data2[2]-data1[2])*24 + data2[3]-data1[3]
    fig2, ax = plt.subplots(figsize=[10, 5])
    ax.set_yticks(range(-90, 90 + 10, 30))
    ax.set_xticks(range(0, ilosc_godz + 2, int(ilosc_godz/10)))

    satelity1 = satelity        #przepisanie i sprawdzenie czy tablica satelit nie jest przypadkiem pusta
    if len(satelity1) == 0:     # w takim przypadku nadajemy jej wartośc wszystkich występujących satelit
        satelity1 = nry
    # petla odpowiedzialna za utworzenie tablicy do rysowania skyplot - skladajaca sie z numeru satelity i jego wszystkich elewacji i azymutow
    # a także tworzenia wykresu zmian elewacji satelitow
    for i in satelity1:
        #index = satelity.index(i)
        sat_positioins.append([i, result[i][:, 2], result[i][:, 1]])
        sat_pos_latlon[i] = np.array(latlon([result[i][:, 6], result[i][:, 7], result[i][:, 8] ]))
        ###wybrane co 6 satelita do wyświetlania jego elewacji, zeby cokolwiek bylo widac na wykresie - przy braniu wszystkich satelit
        #if index%6 == 0:
        ax.plot(epoki, result[i][:, 2], label = i)
        groundtrack(sat_pos_latlon[i], [fi, lb], maska)                 # nie mam pojecia czy to dziala bo mam straszne problemy przy insalacjach wszelkich bibliotek przestrzennych bo nie mam condy

    ax.set(ylim=(0, 95), title= 'Wykres elewacji wybranych satelitów w zależności od czasu')
    ax.set_xlabel('godz. (od początku obserwacji)')
    ax.legend(title='numery\nsatelitów')
    ax.grid(color='gray', axis='y')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[10,6], layout='tight')
    ax1.set_title('wykresy DOP')
    ax1.set_yticks(np.arange(0, 10, 0.5))
    ax1.set_xticks(range(0, ilosc_godz + 2, int(ilosc_godz/10)))
    ax1.plot(epoki, GDOPY, label= 'GDOP')
    ax1.plot(epoki, HDOPY, label= 'HDOP')
    ax1.plot(epoki, VDOPY, label = 'VDOP')
    ax1.plot(epoki, PDOPY, label = 'PDOP')
    ax1.plot(epoki, TDOPY, label= 'TDOP')
    ax1.set_xlabel('godz (od początku obserwacji)')
    ax1.set_ylabel('wartość współczynnika')
    ax1.grid(color='gray')
    ax1.legend()
    ax2.set_title('wykres widoczności satelit które spełniają konieczne warunki (el>maski)')
    ax2.set_yticks(range(0, 32, 2))
    ax2.set_xticks(range(0,ilosc_godz+2 , int(ilosc_godz/10)))
    ax2.stem(epoki, liczbaSat)
    ax2.set_xlabel('godz. (od początku obserwacji)')
    ax2.set_ylabel('ilość widocznych satelitów')
    ax2.grid(color='gray', axis='y')
    plot_skyplot(sat_positioins)

rysowanie()


