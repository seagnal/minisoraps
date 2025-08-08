#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A Python script for passive sonar performance prediction in 1D and constant soundspeed.

Copyright (C) 2025  Johann Baudy

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

mini-SORAPS is a free, lightweight version of SORAPS for public use. It simulates acoustic wave propagation in a 1D marine environment, aiding predictions for sonar and underwater communication systems..
"""

__author__ = "Johann Baudy"
__version__ = "1.0.0"
__license__ = "GPL"

import sys
import argparse
import numpy as np
import math
from PyQt5.QtWidgets import QAction, QScrollArea, QApplication, QMainWindow, QDockWidget, QVBoxLayout, QWidget, QSpinBox, QDoubleSpinBox, QComboBox, QPushButton, QLabel, QSplashScreen, QTextBrowser
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QTimer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Fonction fictive pour l'absorption

def absorption_fisher_simmons(f, T=25, D=100):
    T_kel = 273.1 + T
    P = D / 10.0

    A1 = 1.03e-8 + 2.36e-10 * T - 5.22e-12 * T**2
    P1 = 1
    f1 = 1.32e3 * T_kel * np.exp(-1700 / T_kel)
    Boric = (A1 * P1 * f1 * f**2) / (f**2 + f1**2)

    A2 = 5.62e-8 + 7.52e-10 * T
    P2 = 1 - 10.3e-4 * P + 3.7e-7 * P**2
    f2 = 1.55e7 * T_kel * np.exp(-3052 / T_kel)
    MgSO4 = (A2 * P2 * f2 * f**2) / (f**2 + f2**2)

    A3 = (55.9 - 2.37 * T + 4.77e-2 * T**2 - 3.48e-4 * T**3) * 1e-15
    P3 = 1 - 3.84e-4 * P + 7.57e-8 * P**2
    H2O = A3 * P3 * f**2

    Alpha = (Boric + MgSO4 + H2O) * 20 * np.log10(np.e)
    return Alpha


def calculate_nl(f_sea_state, f_noise, f_margin_shallow_water=0):
    if f_sea_state < 0.5:
        nl = 44 - 17 * np.log10(f_noise / 1000.0) + f_margin_shallow_water
    elif f_sea_state < 1.5:
        nl = 55 - 17 * np.log10(f_noise / 1000.0) + f_margin_shallow_water
    elif f_sea_state < 2.5:
        nl = 61.5 - 17 * np.log10(f_noise / 1000.0) + f_margin_shallow_water
    elif f_sea_state < 3.5:
        nl = 64.5 - 17 * np.log10(f_noise / 1000.0) + f_margin_shallow_water
    elif f_sea_state < 4.5:
        nl = 66.5 - 17 * np.log10(f_noise / 1000.0) + f_margin_shallow_water
    elif f_sea_state < 5.5:
        nl = 68.5 - 17 * np.log10(f_noise / 1000.0) + f_margin_shallow_water
    else:
        nl = 70 - 17 * np.log10(f_noise / 1000.0) + f_margin_shallow_water
    return nl


class MainWindow(QMainWindow):
    def __init__(self, args):
        super().__init__()

        self.init = False
        self.args = args
        self.initUI()

    def initUI(self):
        self.setWindowTitle('mini SORAPS by SEAGNAL')
        self.setGeometry(100, 100, 800, 600)

        # Menu Bar
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        saveAction = QAction('Save', self)
        saveAction.triggered.connect(self.save_plot)
        fileMenu.addAction(saveAction)

        # Dock pour les spinboxes
        dock_config = QDockWidget('Configuration', self)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_config)

        config_widget = QWidget()
        config_layout = QVBoxLayout()

        self.spinboxes = {
            'Fmin': QDoubleSpinBox(),
            'Fmax': QDoubleSpinBox(),
            'SL': QDoubleSpinBox(),
            'SH': QDoubleSpinBox(),
            'B': QDoubleSpinBox(),
            'DI': QDoubleSpinBox(),
            'T': QDoubleSpinBox(),
            'Dmin': QDoubleSpinBox(),
            'Dmax': QDoubleSpinBox(),
            'SeaState': QDoubleSpinBox(),
            'Depth': QDoubleSpinBox(),
            'Temp': QDoubleSpinBox(),
            'DT': QDoubleSpinBox(),
            'PFA': QDoubleSpinBox()
        }

        self.spinboxes['Fmin'].setRange(100, 1000000)
        self.spinboxes['Fmin'].setSuffix(' Hz')
        self.spinboxes['Fmin'].setValue(self.args.Fmin)

        self.spinboxes['Fmax'].setRange(100, 1000000)
        self.spinboxes['Fmax'].setSuffix(' Hz')
        self.spinboxes['Fmax'].setValue(self.args.Fmax)

        self.spinboxes['SL'].setRange(10, 230)
        self.spinboxes['SL'].setSuffix(' dB')
        self.spinboxes['SL'].setValue(self.args.SL)

        self.spinboxes['SH'].setRange(-220, -150)
        self.spinboxes['SH'].setSuffix(' dBV/uPa')
        self.spinboxes['SH'].setValue(self.args.SH)

        self.spinboxes['B'].setRange(0.01, 50000)
        self.spinboxes['B'].setSuffix(' Hz')
        self.spinboxes['B'].setValue(self.args.B)

        self.spinboxes['DI'].setRange(0, 50)
        self.spinboxes['DI'].setSuffix(' dB')
        self.spinboxes['DI'].setValue(self.args.DI)

        self.spinboxes['T'].setRange(0.001, 50)
        self.spinboxes['T'].setSuffix(' s')
        self.spinboxes['T'].setValue(self.args.T)

        self.spinboxes['Dmin'].setRange(1, 100000)
        self.spinboxes['Dmin'].setSuffix(' m')
        self.spinboxes['Dmin'].setValue(self.args.Dmin)

        self.spinboxes['Dmax'].setRange(self.args.Dmin, 100000)
        self.spinboxes['Dmax'].setSuffix(' m')
        self.spinboxes['Dmax'].setValue(self.args.Dmax)

        self.spinboxes['DT'].setRange(0, 50)
        self.spinboxes['DT'].setSuffix(' dB')
        self.spinboxes['DT'].setValue(self.args.DT)
        
        self.spinboxes['PFA'].setRange(-9, -2)
        self.spinboxes['PFA'].setSuffix(' 10**')
        self.spinboxes['PFA'].setValue(self.args.PFA)

        self.spinboxes['SeaState'].setRange(0, 6)
        self.spinboxes['SeaState'].setSuffix('')
        self.spinboxes['SeaState'].setValue(self.args.SeaState)

        self.spinboxes['Depth'].setRange(0, 10000)
        self.spinboxes['Depth'].setSuffix(' m')
        self.spinboxes['Depth'].setValue(self.args.Depth)

        self.spinboxes['Temp'].setRange(0, 60)
        self.spinboxes['Temp'].setSuffix(' degC')
        self.spinboxes['Temp'].setValue(self.args.Temp)

        for name, spinbox in self.spinboxes.items():
            config_layout.addWidget(QLabel(name))
            config_layout.addWidget(spinbox)
            spinbox.valueChanged.connect(self.update_plots)

        self.combo_type = QComboBox()
        self.combo_type.addItems(['Incoherent', 'Coherent'])
        self.combo_type.currentTextChanged.connect(self.update_plots)
        config_layout.addWidget(QLabel('TypeTraitement'))
        self.combo_type.setCurrentText(self.args.TypeTraitement)
        config_layout.addWidget(self.combo_type)

        config_widget.setLayout(config_layout)


        # Ajouter une barre de défilement verticale
        scroll_area = QScrollArea()
        scroll_area.setWidget(config_widget)
        scroll_area.setWidgetResizable(True)
        dock_config.setWidget(scroll_area)
        scroll_area.setMaximumWidth(200)


        # Dock pour les plots ES
        dock_es = QDockWidget('Plots ES', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_es)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout()

        self.figure_es = Figure()
        self.canvas_es = FigureCanvas(self.figure_es)
        plot_layout.addWidget(self.canvas_es)

        plot_widget.setLayout(plot_layout)
        dock_es.setWidget(plot_widget)
        
        # Dock pour les plots POD
        dock_pod = QDockWidget('Plots POD', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_pod)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout()

        self.figure_pod = Figure()
        self.canvas_pod = FigureCanvas(self.figure_pod)
        plot_layout.addWidget(self.canvas_pod)

        plot_widget.setLayout(plot_layout)
        dock_pod.setWidget(plot_widget)

        # Dock pour le plot CP
        dock_cp = QDockWidget('Plot CP', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_cp)

        cp_widget = QWidget()
        cp_layout = QVBoxLayout()

        self.figure_cp = Figure()
        self.canvas_cp = FigureCanvas(self.figure_cp)
        cp_layout.addWidget(self.canvas_cp)

        cp_widget.setLayout(cp_layout)
        dock_cp.setWidget(cp_widget)

        # Dock pour le plot SH
        dock_v = QDockWidget('Plot Tension', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_v)

        sh_widget = QWidget()
        sh_layout = QVBoxLayout()

        self.figure_sh = Figure()
        self.canvas_sh = FigureCanvas(self.figure_sh)
        sh_layout.addWidget(self.canvas_sh)

        sh_widget.setLayout(sh_layout)
        dock_v.setWidget(sh_widget)


        # Dock pour le bar graph de N en fonction de f
        dock_bargraph_N = QDockWidget('Noise', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_bargraph_N)

        bargraph_N_widget = QWidget()
        bargraph_N_layout = QVBoxLayout()

        self.figure_bargraph_N = Figure()
        self.canvas_bargraph_N = FigureCanvas(self.figure_bargraph_N)
        bargraph_N_layout.addWidget(self.canvas_bargraph_N)

        bargraph_N_widget.setLayout(bargraph_N_layout)
        dock_bargraph_N.setWidget(bargraph_N_widget)

        # Dock pour le bar graph de N en fonction de f
        dock_bargraph_A = QDockWidget('Absorption', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_bargraph_A)

        bargraph_A_widget = QWidget()
        bargraph_A_layout = QVBoxLayout()

        self.figure_bargraph_A= Figure()
        self.canvas_bargraph_A = FigureCanvas(self.figure_bargraph_A)
        bargraph_A_layout.addWidget(self.canvas_bargraph_A)

        bargraph_A_widget.setLayout(bargraph_A_layout)
        dock_bargraph_A.setWidget(bargraph_A_widget)

        # Create a QDockWidget
        dock_about = QDockWidget("About SORAPS", self)
        dock_about.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)

        # Create a QTextBrowser to display rich text with a hyperlink
        text_browser = QTextBrowser()
        text_browser.setOpenExternalLinks(True)  # Allow opening external links
        text_browser.setHtml("""
        <p><b>SORAPS (SOnar RAnge Prediction Software)</b> is a powerful software for the performance prediction of systems based on underwater acoustic waves propagation, such as SONAR or digital communication systems.</p>
        <p>This software is interfaced to 3D GHOM (Geography, Hydrography, Oceanography & Meteorology) databases to enable finer estimations of the characteristics of the acoustic channels between a reference point of the water column and many surrounding others, including full 4D simulations (3D+time) and 3D outputs.</p>
        <p>The <b>SORAPS</b> software is designed to predict the performance of sonar and underwater communication systems using 3D GHOM databases and detailed system specifications. These predictions can be obtained in real-time through actual measurements of parameters or in a predictive manner using statistical environmental data from specialized databases and predictive models.</p>

        <p>For more information, on the full featured Software, visit the <a href="http://products.seagnal.fr/soraps">SORAPS product page</a>.</p>
        
        <p>or write us a <a href="mailto:contact@seagnal.fr?subject=SORAPS">mail</a></p>
        
        <p>Author:Johann Baudy.</p>
        <p>Company:<a href="http://www.seagnal.fr">SEAGNAL</a>.</p>
        """)

        dock_about.setWidget(text_browser)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_about)

        self.tabifyDockWidget(dock_about, dock_bargraph_N)
        self.tabifyDockWidget(dock_bargraph_N, dock_bargraph_A)
        self.tabifyDockWidget(dock_bargraph_A, dock_v)
        self.tabifyDockWidget(dock_v, dock_cp)
        self.tabifyDockWidget(dock_cp, dock_pod)
        self.tabifyDockWidget(dock_pod, dock_es)

        self.init = True
        self.update_plots()

        if self.args.save:
            self.save_plot()
            sys.exit(0)

    def update_plots(self):
        if not self.init:
            return
        #print('--- START UPDATE PLOTS')
        Fmin = self.spinboxes['Fmin'].value()
        Fmax = self.spinboxes['Fmax'].value()
        SL = self.spinboxes['SL'].value()
        B = self.spinboxes['B'].value()
        DI = self.spinboxes['DI'].value()
        T = self.spinboxes['T'].value()
        Dmin = self.spinboxes['Dmin'].value()
        Dmax = self.spinboxes['Dmax'].value()
        DT = self.spinboxes['DT'].value()
        SH = self.spinboxes['SH'].value()
        Depth = self.spinboxes['Depth'].value()
        Temp = self.spinboxes['Temp'].value()
        Pfa = self.spinboxes['PFA'].value()
        SeaState = self.spinboxes['SeaState'].value()

        d = np.arange(Dmin, Dmax + 1, 1)
        f = np.transpose(np.array([Fmin, (Fmin + Fmax) / 2, Fmax]))
        np_erf = np.vectorize(math.erf)
        Pfa_ = math.pow(10, Pfa)


        # Calculate the noise level (N) based on the sea state and frequency
        # SeaState: Represents the condition of the sea, affecting ambient noise levels
        # f: Frequency at which the noise level is calculated
        N = calculate_nl(SeaState, f)

        # Calculate the absorption coefficient (alpha) using the Fisher-Simmons model
        # This model estimates sound absorption in seawater based on frequency, depth, and temperature
        # f: Frequency of the sound wave
        # D: Depth of the water, set to the variable Depth
        # T: Temperature of the water, set to the variable Temp
        alpha = absorption_fisher_simmons(f, D=Depth, T=Temp)

        # CP: Computation of Propagation Loss
        # This line calculates the propagation loss (CP) which includes spherical spreading loss and absorption loss.
        # 20 * np.log10(d[:, np.newaxis]) represents spherical spreading loss.
        # alpha * d[:, np.newaxis] represents absorption loss, where alpha is the absorption coefficient.
        CP = 20 * np.log10(d[:, np.newaxis]) + alpha * d[:, np.newaxis]

        # V: Voltage or signal level after propagation
        # This line calculates the voltage or signal level (V) at the receiver.
        # SL is the Source Level, CP is the Propagation Loss, and SH is the Signal Excess or Hydrophone Sensitivity.
        V = 10.0**((SL - CP + SH)/20)

        # ES: Signal Excess for Coherent and Incoherent Processing
        # Depending on the type of processing (Coherent or Incoherent), the Signal Excess (ES) is calculated differently.
        if self.combo_type.currentText() == 'Coherent':
            # For Coherent processing, the Signal Excess includes terms for Bandwidth (B) and Integration Time (T).
            ES = SL - CP + 10.0 * np.log10(B) + 10.0 * np.log10(T) - (N - DI) - 10.0 * np.log10(B) - DT
            E_N = (N - DI) - 10.0 * np.log10(T)
        else:# self.combo_type.currentText() == 'Incoherent':
            # For Incoherent processing, the Signal Excess calculation is slightly different.
            ES = SL - CP + 10.0 * np.log10(T) - (N-DI) - 10.0 * np.log10(B) - DT
            E_N = (N - DI) - 10.0 * np.log10(T) + 10.0 * np.log10(B)

        # Convert ES and E_N from dB to linear scale
        E_S = 10.**((SL - CP)/10);
        E_N = 10.**(E_N/10);


        # Calculate detection threshold (dt_s) based on the noise level and probability of false alarm
        dt_s= E_N *np.sqrt(-np.log(Pfa_));
        
        # This uses the error function to compute the probability of detection based on the detection threshold, signal excess, and noise level.
        POD = 0.5 - np_erf((dt_s - E_S) / (np.sqrt(2) * E_N)) / 2

        # Plot POD
        self.figure_pod.clear()
        ax = self.figure_pod.add_subplot(111)
        ax.plot(d, POD[:, 0], label='Fmin', color='blue')
        ax.plot(d, POD[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax.plot(d, POD[:, 2], label='Fmax', color='red')
        ax.grid()
        ax.legend()
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('POD (dB)')
        ax.set_title(f'POD {self.combo_type.currentText()}')
        ax.set_ylim(-0.1, 1.1)
        self.canvas_pod.draw()

        # Plot ES
        self.figure_es.clear()
        ax = self.figure_es.add_subplot(111)
        ax.plot(d, ES[:, 0], label='Fmin', color='blue')
        ax.plot(d, ES[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax.plot(d, ES[:, 2], label='Fmax', color='red')
        ax.grid()
        ax.legend()
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('ES (dB)')
        ax.set_title(f'ES {self.combo_type.currentText()}')
        self.canvas_es.draw()

        # Plot CP
        self.figure_cp.clear()
        ax_cp = self.figure_cp.add_subplot(111)
        ax_cp.plot(d, CP[:, 0], label='Fmin', color='blue')
        ax_cp.plot(d, CP[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax_cp.plot(d, CP[:, 2], label='Fmax', color='red')
        ax_cp.legend()
        ax_cp.grid()
        ax_cp.set_xlabel('Distance (m)')
        ax_cp.set_ylabel('Champ de perte (dB)')
        ax_cp.set_title(f'Champ de perte')
        self.canvas_cp.draw()

        # Plot SH
        self.figure_sh.clear()
        ax_sh = self.figure_sh.add_subplot(111)
        ax_sh.semilogy(d, V[:,0], label='Fmin', color='blue')
        ax_sh.semilogy(d, V[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax_sh.semilogy(d, V[:, 2], label='Fmax', color='red')
        ax_sh.legend()
        ax_sh.grid()
        ax_sh.set_xlabel('Distance (m)')
        ax_sh.set_ylabel('Vrms')
        ax_sh.set_title(f'Voltage at sensor')
        self.canvas_sh.draw()

        # Bar Graph N vs f
        self.figure_bargraph_N.clear()
        ax_bargraph_N = self.figure_bargraph_N.add_subplot(111)
        ax_bargraph_N.bar(['Fmin','Fmoy','Fmax'], N, color='blue')
        ax_bargraph_N.set_xlabel('Frequency (Hz)')
        ax_bargraph_N.set_ylabel('N (dBuPa/sqrt(Hz))')
        ax_bargraph_N.set_title(f'Noise over Bandwidth - SeaState:{SeaState}')

        ax_bargraph_N.grid(axis='y')
        self.canvas_bargraph_N.draw()

        # Bar Graph N vs f
        self.figure_bargraph_A.clear()
        ax_bargraph_A = self.figure_bargraph_A.add_subplot(111)
        ax_bargraph_A.bar(['Fmin','Fmoy','Fmax'], alpha, color='blue')
        ax_bargraph_A.set_xlabel('Frequency (Hz)')
        ax_bargraph_A.set_ylabel('absorption db/m')
        ax_bargraph_A.set_title(f'Absorption over Bandwidth - T:{Temp} - D:{Depth}')
        ax_bargraph_A.grid(axis='y')
        self.canvas_bargraph_A.draw()



        #print('--- END UPDATE PLOTS')

    def save_plot(self):
        print('--- START SAVE PLOTS')
        self.figure_pod.savefig('plot_pod.png')
        self.figure_es.savefig('plot_es.png')
        self.figure_cp.savefig('plot_cp.png')
        self.figure_sh.savefig('plot_sh.png')
        self.figure_bargraph_N.savefig('plot_bargraph_N.png')
        self.figure_bargraph_A.savefig('plot_bargraph_A.png')

        # Open a file for writing
        with open('variables.txt', 'w') as file:
            for name, spinbox in self.spinboxes.items():
                # Extract the unit from the spinbox suffix
                unit = spinbox.suffix().strip()
                value = spinbox.value()
                file.write(f"{name}: {value} {unit}\n")

        print('--- END UPDATE PLOTS')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Mini-Soraps')
    parser.add_argument('--Fmin', type=float, default=10000, help='Fmin value (Hz)')
    parser.add_argument('--Fmax', type=float, default=20000, help='Fmax value (Hz)')
    parser.add_argument('--SL', type=float, default=150, help='SL value (dBuPa)')
    parser.add_argument('--SH', type=float, default=-200, help='SH value (dBV/uPa)')
    parser.add_argument('--B', type=float, default=500, help='B value (Hz)')
    parser.add_argument('--DI', type=float, default=0, help='DI value (dB)')
    parser.add_argument('--T', type=float, default=0.01, help='T value (s)')
    parser.add_argument('--Dmin', type=float, default=10, help='Dmin value (m)')
    parser.add_argument('--Dmax', type=float, default=3000, help='Dmax value (m)')
    parser.add_argument('--DT', type=float, default=0, help='DT value (dB)')
    parser.add_argument('--PFA', type=float, default=-7, help='PFA value (10**)')
    parser.add_argument('--SeaState', type=int, default=1, help='SeaState value (0-6)')
    parser.add_argument('--Temp', type=float, default=15, help='Temp value (degC)')
    parser.add_argument('--Depth', type=float, default=100, help='Depth value (m)')
    parser.add_argument('--TypeTraitement', type=str, choices=['Incoherent', 'Coherent'], default='Coherent', help='Type of treatment')
    parser.add_argument('--save', action='store_true', default=False, help='Save and Exit')
    return parser.parse_args()


def signal_handler(sig, frame):
    print('Control-C detected. Exiting gracefully...')
    QApplication.quit()

import signal
signal.signal(signal.SIGINT, signal_handler)

if __name__ == '__main__':
    args = parse_arguments()
    app = QApplication(sys.argv)

    # Create a splash screen
    splash_pix = QPixmap('seagnal.png')  # Make sure to have a logo.png file in your working directory
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()

    # Simulate some loading time
    QTimer.singleShot(3000, splash.close)  # Close the splash screen after 3 seconds

    mainWin = MainWindow(args)
    mainWin.showMaximized()
    sys.exit(app.exec_())
