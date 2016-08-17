__author__ = "Ben"

import numpy as np
# Plot tools
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import animation




### Unit
GHz = 10.0**9
MHz = 10.0**6
KHz = 10.0**3

def peak_width(array, center_idx=0, boundary_cond=3, span=0, mode='Down', reverse=False):
    """
    Get peak width by boundary condition
    Examples:
    FWHM(boundary_cond=0.5, mode='Multiple')
    3dB down(boundary_cond=3, mode='Down')

    @param array: 1-D array
    @param center_idx: Index of peak, if 0 auto peak
    @param boundary_cond: Boundary condition of edge
    @param span: Analysis span
    @param mode: Down or Multiple(boundary_cond < 1)
    @param reverse: True if peak is below
    @return: Index of edge
    """

    length = np.size(array)
    # Auto peak
    if center_idx == 0:
        center_idx = np.argmax(array, axis=0) if reverse is False else np.argmin(array, axis=0)

    # Select mode
    if mode == 'Down':
        height = array[center_idx] - boundary_cond * (1 - 2 * reverse)
    elif mode == 'Multiple':
        if boundary_cond > 1 or boundary_cond < 0: ValueError('boundary_cond should between 0 and 1')
        height = array[center_idx] * boundary_cond
    else:
        ValueError('mode = Down/Multiple')

    if span == 0:
        idx1 = np.abs(np.subtract.outer(array[0:center_idx:1], height)).argmin(0)
        idx2 = np.abs(np.subtract.outer(array[center_idx:length:1], height)).argmin(0)
    else:
        idx1 = np.abs(np.subtract.outer(array[center_idx - span / 2:center_idx:1], height)).argmin(0)
        idx2 = np.abs(np.subtract.outer(array[center_idx:center_idx + span / 2:1], height)).argmin(0)

    return idx1, idx2+center_idx







class Resonator:
    def __init__(self, LRCG_Ls, LRCG_Cs, cpw_pinw, cpw_gapw):
        self.Ls = LRCG_Ls
        self.Cs = LRCG_Cs
        self.pinw = cpw_pinw
        self.gapw = cpw_gapw

    def LengthCal(self, Freq, ShortEnd=False):

        x = np.linspace(0, 1, 100)
        if ShortEnd == True:
            P = 2
            y = np.sin(x * np.pi / 2 + np.pi / 2)
        else:
            P = 1
            y = np.cos(x * np.pi)
        Length = 1 / (np.sqrt(self.Cs * self.Ls) * 2 * Freq * 1000 * P)

        print 'For n=1 Frequency = %.2f GHz, Resonator Length = %.2f um ' % (Freq, Length)
        fig, ax = plt.subplots(figsize=(12, 2))
        cmap = mpl.cm.cool
        for i in np.linspace(1, -1, 20):
            plt.plot(x * Length, y * i, color=cmap(i / 2 + 0.5))
        plt.xticks(np.linspace(0, Length, 3))
        plt.xlim(0, Length)
        plt.ylim(-1, 1)
        plt.xlabel('Position is Resonator(um)')
        plt.ylabel('Potential (100%)')

    def Qfactor_simulator(self, Type='Hanger',  Qo=10000, Qe=10000, Swp_F=100):
        Freq = 5000
        if Type == 'Hanger':
            self.Qfactor_hanger(Freq, Qo, Qe, Swp_F, offset=True)
        if Type == 'TwoPort':
            print 'Not Finish Yet !'


    def Qfactor_hanger(self, r_freq=5, Unloaded_Q=438000, Couple_Q=43000, sweep_freq=1000, sweep_point = 10000, d_freq = 0,offset = False, autorange = False):
        """

        @param r_freq: Resonance Frequency (MHz)
        @param Unloaded_Q: Unloaded Q factor ; Internal Q factor
        @param Couple_Q: Couple Q factor ; External Q factor
        @param sweep_freq: Span of Simulate Frequency (MHz)
        @param sweep_point: Point of Simulate Frequency
        @param d_freq: asymmetry parameter (KHz)
        @param offset: move resonance frequency to middle if True
        @param autorange: autorange if True
        @return:
        """
        fig, ax = plt.subplots(figsize=(6, 4))
        r_freq = r_freq*MHz
        d_freq = d_freq*KHz
        sweep_freq = sweep_freq*1.0
        Loaded_Q = (1. / ((1. / Unloaded_Q) + (1. / Couple_Q)))
        autorange_goal = 0.1
        offset_center = r_freq if offset is True else 0

        resolution = 0
        while resolution <= autorange_goal:
            Sweep_freq = np.linspace(r_freq - sweep_freq * 10 ** 6 / 2, r_freq + sweep_freq * 10 ** 6 / 2, sweep_point)
            S21 = np.zeros(sweep_point)
            n = 0
            for i in Sweep_freq:
                s_freq = i
                a = (s_freq - (r_freq + d_freq))/(r_freq+d_freq)
                b = 2 * d_freq / r_freq
                S21[n] = (-2*Loaded_Q*Couple_Q+Couple_Q**2+(Loaded_Q**2)*(1+(Couple_Q**2)*((2*a+b)**2)))/((Couple_Q**2)*(1 +4*(Loaded_Q**2)*(a**2)))
                n += 1
            S21_dB = 10*np.log10(S21)

            #autorange
            if autorange is True:
                idx = peak_width(S21_dB, boundary_cond=0.8, mode='Multiple', reverse=True)
                bandwidth = abs(Sweep_freq[idx[1]] - Sweep_freq[idx[0]])
                resolution = bandwidth / (sweep_freq * 10 ** 6)
                sweep_freq = sweep_freq/10
            else: break

        #     S21_0 = np.min(S21)
        #     g = S21_0/(1-S21_0)
        #
        #
        # print S21_0
        # print 'Theory Loaded Q = ' + str(Loaded_Q)
        # print 'Theory g = ' + str(Unloaded_Q/Couple_Q)
        # print 'Measure g = ' + str(g)
        #
        # if (S21_dB.max() - S21_dB.min()) > 3:
        #     idx = peak_width(S21_dB, sweep_point / 2, 3, reverse=True)
        #     bandwidth = abs(Sweep_freq[idx[1]] - Sweep_freq[idx[0]])
        #     Cal_QL = r_freq / bandwidth
        #     print r_freq, bandwidth
        #     print 'Calculate Loaded Q = %.2f' % (Cal_QL)

        Sweep_freq = Sweep_freq-offset_center
        if abs(Sweep_freq[-1]-Sweep_freq[0])*10 >= GHz:
            Sweep_freq = Sweep_freq/GHz
            plt.xlabel('Frequency (GHz)')
        elif abs(Sweep_freq[-1]-Sweep_freq[0])*10 >= MHz:
            Sweep_freq = Sweep_freq/MHz
            plt.xlabel('Frequency (MHz)')
        elif abs(Sweep_freq[-1] - Sweep_freq[0])*10 >= KHz:
            Sweep_freq = Sweep_freq / KHz
            plt.xlabel('Frequency (KHz)')
        else:
            plt.xlabel('Frequency (Hz)')


        plt.plot(Sweep_freq, S21_dB, color='k')
        plt.title('QL = %s, g = %.2f' % (Loaded_Q,1))

        ax.get_xaxis().get_major_formatter().set_useOffset(False)

        plt.ylabel('S21 (dB)')
        plt.xlim((Sweep_freq[0], Sweep_freq[-1]))

        plt.xticks(np.linspace(Sweep_freq[0], Sweep_freq[-1], 5))
        plt.autoscale(axis='y')
        plt.show()


    def Qfactor_twoport(self, A = []):
        fig, ax = plt.subplots(figsize=(8, 8))
        x = np.linspace(-100,100,1000)
        y=np.zeros(1000)
        n=0
        for i in x:
            y[n] = A[0] + A[1]*i/(1+((A[2]*i*i)**2/(A[3]*i*i*i)**2))
            n+=1

        plt.plot(x,y)
        plt.show()










if __name__ == "__main__":
    cpw_pinw = 10.0
    cpw_gapw = 4.8
    LRCG_Ls = 397.161 * 10 ** -9
    LRCG_Cs = 158.297 * 10 ** -12

    Res = Resonator(LRCG_Ls,LRCG_Cs,cpw_pinw,cpw_gapw)
    Res.Qfactor_hanger(Unloaded_Q=430000,Couple_Q=43000,sweep_freq=10,autorange=True,offset=True)
