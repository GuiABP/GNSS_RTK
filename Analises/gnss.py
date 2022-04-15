import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

class gnss:

    # SIRGAS2000
    a = 6378137 # semi-eixo maior
    b = 6356752.298 # semi eixo menor 
    f = 1/298.25722356 # achatamento do elipsoide -> f = (a-b)/a
    e2 = 0.006694385 # primeira excentricidade ao quadrado -> e**2 = (a**2 - b**2)/(a**2) ou e**2 = 2*f - f**2
    e2l = (a**2 - b**2)/(b**2) # segunda excentricidade ao quadrado (e2l = e'**2)

    def __init__(self, phiDeg, lbdDeg):

        self.phiDeg = phiDeg
        self.lbdDeg = lbdDeg

        self.phiRad = np.deg2rad(self.phiDeg)
        self.lbdRad = np.deg2rad(self.lbdDeg)

        self.length = len(self.phiDeg)

        self.h = np.zeros(self.length)


    def geodetic2ecef(self, phiRad, lbdRad, h):
        if isinstance(phiRad, int) or isinstance(phiRad, float):
            N = self.a/((1 - self.e2*(np.sin(phiRad)) ** 2) ** 0.5)
            x = ((N + h)*np.cos(phiRad)*np.cos(lbdRad))
            y = ((N + h)*np.cos(phiRad)*np.sin(lbdRad))
            z = ((N*(1-self.e2)+h)*np.sin(phiRad))
        else:
            length = len(phiRad)
        
            x = np.zeros(length)
            y = np.zeros(length)
            z = np.zeros(length)
            
            for i in range(self.length):
                N = self.a/((1 - self.e2*(np.sin(phiRad[i])) ** 2) ** 0.5)
                x[i] = ((N + h[i])*np.cos(phiRad[i])*np.cos(lbdRad[i]))
                y[i] = ((N + h[i])*np.cos(phiRad[i])*np.sin(lbdRad[i]))
                z[i] = ((N*(1-self.e2)+h[i])*np.sin(phiRad[i]))
        
        return (x, y, z)

    def ecef2enu(self, x, y, z, x0, y0, z0):
        p = np.sqrt(x0**2 + y0**2) # parametro p - para simplificacao das expressoes de latitude e longitude

        mLbdRad = (np.arctan2(y0, x0)) # calculo da longitude em rad

        theta = np.arctan2((z0*self.a),(p*self.b)) # parametro theta - para simplificacao da expressao da latitude
        mPhiNum = z0 + self.e2l*self.b*np.sin(theta)**3 # numerador da expressao para calculo do valor medio da latitude
        mPhiDen = p - self.e2*self.a*np.cos(theta)**3 # numerador da expressao para calculo do valor medio latitude

        mPhiRad = np.arctan2(mPhiNum, mPhiDen) # expressao para calculo do valor medio da latitude

        mN = self.a/(1 - self.e2*(np.sin(mPhiRad) ** 2)) ** 0.5 # grande normal utilizando o valor medio da latitude

        phi0 = mPhiRad
        lbd0 = mLbdRad

        # conversao de coordenadas cartersianas geocentricas em coordenadas cartesianas locais (sistema geodesico local)
        if isinstance(x, int) or isinstance(x, float):
            east = -np.sin(lbd0)*(x - x0) + np.cos(lbd0)*(y - y0)
            north = -np.sin(phi0)*np.cos(lbd0)*(x - x0) - np.sin(phi0)*np.sin(lbd0)*(y - y0) + np.cos(phi0)*(z - z0)
            up = np.cos(phi0)*np.cos(lbd0)*(x - x0) + np.cos(phi0)*np.sin(lbd0)*(y - y0) + np.sin(phi0)*(z - z0)
        else:
            length = len(x)
            # inicializacao dos vetores e, n e u
            east = np.zeros(length)
            north = np.zeros(length)
            up = np.zeros(length)
            for i in range(length):
                east[i] = -np.sin(lbd0)*(x[i] - x0) + np.cos(lbd0)*(y[i] - y0)
                north[i] = -np.sin(phi0)*np.cos(lbd0)*(x[i] - x0) - np.sin(phi0)*np.sin(lbd0)*(y[i] - y0) + np.cos(phi0)*(z[i] - z0)
                up[i] = np.cos(phi0)*np.cos(lbd0)*(x[i] - x0) + np.cos(phi0)*np.sin(lbd0)*(y[i] - y0) + np.sin(phi0)*(z[i] - z0)

        return (east, north, up)


    def selfGeodetic2ecef(self):
        ecef = self.geodetic2ecef(self.phiRad, self.lbdRad, self.h)
        
        return ecef

    def histPlot(self, eDisc, nDisc, uDisc, eSdtDisc, nStdDisc, uStdDisc):
        bars = 20
        plt.figure(figsize=(5, 5))
        plt.subplot(311)
        count, bins, ignored = plt.hist(eDisc, bars, density=True)
        plt.plot(bins, 1/np.sqrt(2 * np.pi * eSdtDisc**2)*np.exp(-(bins-np.mean(eDisc))**2 / (2*eSdtDisc**2) ),
                linewidth=2, color='r')
        plt.legend(["Normal"], loc=2, prop={'size': 6})
        plt.ylabel('E')
        plt.grid()

        plt.subplot(312)
        count, bins, ignored = plt.hist(nDisc, bars, density=True)
        plt.plot(bins, 1/np.sqrt(2 * np.pi * nStdDisc**2)*np.exp(-(bins-np.mean(nDisc))**2 / (2*nStdDisc**2) ),
                linewidth=2, color='r')
        plt.legend(["Normal"], loc=2, prop={'size': 6})
        plt.ylabel('N')
        plt.grid() 

        plt.subplot(313)
        count, bins, ignored = plt.hist(uDisc, bars, density=True)
        plt.plot(bins, 1/np.sqrt(2 * np.pi * uStdDisc**2)*np.exp(-(bins-np.mean(uDisc))**2 / (2*uStdDisc**2) ),
                linewidth=2, color='r')
        plt.legend(["Normal"], loc=2, prop={'size': 6})
        plt.grid()
        plt.ylabel('U')
        plt.xlabel('Discrepância (m)')

    def enuPlot(self, e, n, u):
        plt.figure()

        plt.subplot(311)
        plt.plot(e, 'c')
        plt.ylabel('e (m)')
        plt.grid()

        plt.subplot(312)
        plt.plot(n, 'c')
        plt.ylabel('n (m)')
        plt.grid()

        plt.subplot(313)
        plt.plot(u, 'c')
        plt.xlabel('Tempo (s)')
        plt.ylabel('u (m)')
        plt.grid()
        
        plt.suptitle("Série temporal de posição")

    def discPlot(self, e, n, drm):
        plt.figure()
        # now make a circle with no fill, which is good for hi-lighting key results
        circle1 = plt.Circle((np.mean(e), np.mean(n)), drm, color='r', fill=False)
        circle2 = plt.Circle((np.mean(e), np.mean(n)), 2*drm, color='r', fill=False)
            
        ax = plt.gca()
        ax.cla() # clear things for fresh plot

        # change default range so that new circles will work
        #ax.set_xlim((-0.04, 0.04))
        #ax.set_ylim((-0.04, 0.04))
        # key data point that we are encircling
        ax.scatter(e, n, marker = ".")

        ax.set_aspect('equal', adjustable='box')
        ax.add_patch(circle1)
        ax.add_patch(circle2)

        plt.title('Discrepância planimétrica')
        plt.xlabel('Leste (m)')
        plt.ylabel('Norte (m)')

    def enuCinematicoPlot(self, e, n):
        plt.figure()

        plt.subplot(211)
        plt.plot(e, 'c')
        plt.ylabel('e (m)')
        plt.yticks([-15, -10, -5, 0, 5, 10, 15])
        plt.grid()

        plt.subplot(212)
        plt.plot(n, 'c')
        plt.xlabel('Tempo (s)')
        plt.ylabel('n (m)')
        plt.yticks([-15, -10, -5, 0, 5, 10, 15])
        plt.grid()

        #plt.suptitle("Série temporal do teste cinemático")

    def cinematicoPlot(self, e, n):
        plt.figure()
            
        ax = plt.gca()
        ax.cla() # clear things for fresh plot

        # change default range so that new circles will work
        #ax.set_xlim((-0.04, 0.04))
        #ax.set_ylim((-0.04, 0.04))
        # key data point that we are encircling
        ax.scatter(e, n, marker = ".")

        ax.set_aspect('equal', adjustable='box')

        #plt.grid()
        #plt.title('Trajetória percorrida no teste cinemático')
        plt.xlabel('Leste (m)')
        plt.ylabel('Norte (m)')

    def std(self, eDisc, nDisc, uDisc):
        eStdDisc = np.std(eDisc)
        nStdDisc = np.std(nDisc)
        uStdDisc = np.std(uDisc)
        return (eStdDisc, nStdDisc, uStdDisc)

    def drm(self, eStdDisc, nStdDisc):
        # distance root mean squared - expressa a acuracia 2D
        # referencia - GPS Position Accuracy Measures NOVATEL
        drm = np.sqrt(eStdDisc ** 2 + nStdDisc ** 2)
        return drm

    def accuracy(self, eDisc, nDisc, uDisc):
        # variavel para armazenar a soma dos quadrados das discrepancias
        sum_east_square = 0
        sum_north_square = 0
        sum_up_square = 0

        # somatorio dos quadrados das discrepancias
        for i in range(self.length):
            sum_east_square = sum_east_square + (eDisc[i]) ** 2
            sum_north_square = sum_north_square + (nDisc[i]) ** 2
            sum_up_square = sum_up_square + (uDisc[i]) ** 2

        # root mean square erro (RMSE) ou raiz do erro medio quadratico (REMQ)
        eReqm = np.sqrt((1/self.length)*sum_east_square)
        nReqm = np.sqrt((1/self.length)*sum_north_square)
        uReqm = np.sqrt((1/self.length)*sum_up_square)
        enAcc = np.sqrt((1/self.length)*(sum_east_square + sum_north_square))
        enuAcc = np.sqrt((1/self.length)*(sum_east_square + sum_north_square + sum_up_square))
        #print("eAcc: ", eReqm, " / nAcc: ", nReqm, " / uAcc: ", uReqm, " / enAcc: ", enAcc, " / enuAcc: ", enuAcc)
        return (eReqm, nReqm, uReqm, enAcc, enuAcc)

    def precision(self, e, n, u):
        eP = np.std(e)
        nP = np.std(n)
        uP = np.std(u)

        # variavel para armazenar a soma dos quadrados das discrepancias
        sum_east_square = 0
        sum_north_square = 0

        # somatorio dos quadrados das discrepancias
        for i in range(self.length):
            sum_east_square = sum_east_square + (e[i] - np.mean(e)) ** 2
            sum_north_square = sum_north_square + (n[i] - np.mean(n)) ** 2
        
        enP = np.sqrt((1/self.length)*(sum_east_square + sum_north_square))
        return (eP, nP, uP, enP)

    def printLatexAcuracia(self, name, e, n, u):
        [eReqm, nReqm, uReqm, enAcc, enuAcc] = self.accuracy(e, n, u)
        print("\\hline")
        print(name, " & ", "{:.3f}".format(eReqm), " & ", "{:.3f}".format(nReqm), " & ", "{:.3f}".format(uReqm), " & ", 
        "{:.3f}".format(enAcc), " & ", "{:.3f}".format(enuAcc), " \\\\")

    def printLatexPrecisao(self, name, e, n, u):
        [ep, np, up, enp] = self.precision(e, n, u)
        print("\\hline")
        print(name, " & ", "{:.3f}".format(ep), " & ", "{:.3f}".format(np), " & ", "{:.3f}".format(up), " & ", "{:.3f}".format(enp), " \\\\")

        

