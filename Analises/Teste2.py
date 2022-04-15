import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from datetime import datetime, date, time
import pandas as pd
from gnss import *

rtkData = loadmat('20210529184255r.mat')
boatData = loadmat('00000124.mat')

latitude = rtkData['GPS']['Latitude'][0][0]
longitude = rtkData['GPS']['Longitude'][0][0]
gpsSamples = len(latitude)

boatLatitude = boatData['GPS'][:, 7]
boatLongitude = boatData['GPS'][:, 8]
boatSamples = len(boatLatitude)

lat = np.zeros(gpsSamples)
lon = np.zeros(gpsSamples)

boatLat = np.zeros(boatSamples)
boatLon = np.zeros(boatSamples)

# Converte arrays dentro do array para float
for i in range(0, gpsSamples):
    lat[i] = latitude[i].astype(float)
    lon[i] = longitude[i].astype(float)

for i in range(0, boatSamples):
    boatLat[i] = boatLatitude[i].astype(float)
    boatLon[i] = boatLongitude[i].astype(float)

# inicio do log do mission planner equivalente ao do adcp
startGps = 785
# frequencia do adcp é maior. Usa esse valor para equiparar os dados
freqGps = 5
tGps = np.zeros(gpsSamples)
for i in range(0, gpsSamples):
    tGps[i] = i * freqGps

startBoat = startGps
# valor do final do teste
stopBoat = 1300
tBoat = np.arange(startBoat, stopBoat, 1)

# conversao para enu
gnssData = gnss(lat, lon)
(x, y, z) = gnssData.selfGeodetic2ecef()
(x0, y0, z0) = gnssData.geodetic2ecef(gnssData.phiRad[0], gnssData.lbdRad[0], gnssData.h[0])
(e, n, u) = gnssData.ecef2enu(x, y, z, x0, y0, z0)

# conversao para enu
boatData = gnss(boatLat, boatLon)
(xBoat, yBoat, zBoat) = boatData.selfGeodetic2ecef()
(x0Boat, y0Boat, z0Boat) = gnssData.geodetic2ecef(gnssData.phiRad[0], gnssData.lbdRad[0], gnssData.h[0])
(eBoat, nBoat, uBoat) = boatData.ecef2enu(xBoat, yBoat, zBoat, x0Boat, y0Boat, z0Boat)

# devolve tempo para escala de 1 em 1 s
tG = tGps/freqGps
tB = (tBoat - startBoat)/freqGps
# pega somente o intervalo do barco em que o teste com o rtk ocorreu
eB = eBoat[startBoat:stopBoat]
nB = nBoat[startBoat:stopBoat]

# Downsample para deixar tamanho dos arrays do rtk e do barco iguais
eb = np.zeros(gpsSamples)
nb = np.zeros(gpsSamples)
c = 0
for i in range(len(eB)):
    if tB[i] % 1 == 0 and c < gpsSamples:
        eb[c] = eB[i]
        nb[c] = nB[i]
        c = c + 1

# remove ultimos elementos para tirar casas não varridas pelo downsample no fim do vetor
tG = tG[0:-5]
e = e[0:-5]
n = n[0:-5]
eb = eb[0:-5]
nb = nb[0:-5]

# Calculo do erro
erroE = np.sqrt((e-eb)**2)

sumPlan = 0
sumE = 0
sumN = 0
for i in range(len(e)):
    sumPlan = sumPlan + (eb[i] - e[i])**2 + (nb[i] - n[i])**2
    sumE = sumE + (eb[i] - e[i])**2
    sumN = sumN + (nb[i] - n[i])**2

reqmPlanimetrico = np.sqrt(sumPlan/len(e))
reqmE = np.sqrt(sumE/len(e))
reqmN = np.sqrt(sumN/len(e))

print("reqm: ", reqmPlanimetrico)

plt.figure()
plt.plot(tG, e, label = 'RTK')
plt.plot(tG, eb, label = 'Ponto simples')
plt.legend(loc='upper right')
plt.ylabel('e (m)')
plt.xlabel('Tempo (s)')
plt.grid()

plt.figure()
plt.plot(tG, n, label = 'RTK')
plt.plot(tG, nb, label = 'Ponto simples')
plt.legend(loc='upper right')
plt.ylabel('n (m)')
plt.xlabel('Tempo (s)')
plt.grid()

#plt.figure()
#plt.plot(e-eb)

#plt.figure()
#plt.plot(n-nb)

# plot do erro planimetrico
plt.figure()
plt.plot(erroE)
plt.ylabel('Erro de posicionamento simples - componente e (m)')
plt.xlabel('Tempo (s)')
plt.grid()

plt.figure()
plt.plot(e, n)

plt.show()