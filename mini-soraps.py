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
from PyQt5.QtWidgets import QStatusBar, QAction, QScrollArea, QApplication, QMainWindow, QDockWidget, QVBoxLayout, QGridLayout, QWidget, QSpinBox, QDoubleSpinBox, QComboBox, QPushButton, QLabel, QSplashScreen, QTextBrowser
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QTimer
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor

linestyle_tuple = {
     'loosely dotted':        (0, (1, 10)),
     'dotted':                (0, (1, 5)),
    'densely dotted':        (0, (1, 1)),

     'long dash with offset': (5, (10, 3)),
     'loosely dashed':        (0, (5, 10)),
     'dashed':                (0, (5, 5)),
     'densely dashed':        (0, (5, 1)),

     'loosely dashdotted':    (0, (3, 10, 1, 10)),
     'dashdotted':            (0, (3, 5, 1, 5)),
     'densely dashdotted':    (0, (3, 1, 1, 1)),

     'dashdotdotted':         (0, (3, 5, 1, 5, 1, 5)),
     'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
     'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1))
}

def inverse_erf(y, max_iter=10, tol=1e-12):
    """
    Approximation de l'inverse de la fonction d'erreur.
    
    Parameters
    ----------
    y : array_like
        Valeurs dans l'intervalle [-1, 1] (ou plus strictement [-1+ε, 1-ε] pour éviter l'extremité).
    max_iter : int, optional
        Nombre maximal d'itérations de Newton.
    tol : float, optional
        Tolérance de convergence (abs(Δx) < tol).
    
    Returns
    -------
    x : ndarray
        Valeur de x telle que erf(x) ≈ y.
    """
    np_erf = np.vectorize(math.erf)
    y = np.asarray(y, dtype=float)

    # Cas particulier : y==±1 => x→±∞, on limite artificiellement
    y = np.clip(y, -0.9999999999999999, 0.9999999999999999)

    # Initialisation : approximation simple
    x = np.sqrt(np.pi)/2 * y   # x ≈ (√π/2) y

    for _ in range(max_iter):
        # Calculer erf(x) et sa dérivée
        erf_x = np_erf(x)
        deriv = 2/np.sqrt(np.pi) * np.exp(-x**2)

        # Correction de Newton
        delta = (erf_x - y) / deriv
        x -= delta

        # Vérifier convergence
        if np.all(np.abs(delta) < tol):
            break

    return x

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


def circular_arc_length(R, theta_deg):
    """
    Calculate the arc length of a circle defined by a radius and an angular opening.

    Args:
        R (float or np.ndarray): Radius of the circle (in meters or any length unit).
        theta_deg (float or np.ndarray): Angular opening in degrees.

    Returns:
        float or np.ndarray: Arc length (in meters or any length unit).
    """
    arc_length = 2 * np.pi * R * (theta_deg / 360)
    return arc_length


def spherical_cap_area_vectorized(R, delta_phi_deg, delta_theta_deg):
    """
    Calculate the area of spherical caps for vectorized inputs.

    Args:
        R (np.ndarray): Array of sphere radii (in meters or any length unit).
        delta_phi_deg (np.ndarray): Array of azimuth (bearing) openings in degrees.
        delta_theta_deg (np.ndarray): Array of elevation (site) openings in degrees.

    Returns:
        np.ndarray: Array of spherical cap areas (in m² or any area unit).
    """
    # Convert angles from degrees to radians
    delta_phi_rad = np.radians(delta_phi_deg)
    delta_theta_rad = np.radians(delta_theta_deg)

    # Approximation for small angles: solid angle Omega ≈ delta_phi * delta_theta
    omega = delta_phi_rad * delta_theta_rad

    # Area of the spherical caps
    area = R**2 * omega

    return area

def distance_and_grazing_angle(horizontal_distance, height, angle_in_degrees=True):
    """
    Calculate the Euclidean distance and the grazing angle (angle between the sound ray and the boundary).

    Args:
        horizontal_distance (float or np.ndarray): Horizontal distance (in meters or any length unit).
        height (float or np.ndarray): Height (in meters or any length unit).
        angle_in_degrees (bool): If True, returns the grazing angle in degrees. Otherwise, in radians. Default: True.

    Returns:
        tuple: (euclidean_distance, grazing_angle)
            - euclidean_distance (float or np.ndarray): Euclidean distance (in meters or any length unit).
            - grazing_angle (float or np.ndarray): Grazing angle between the sound ray and the boundary (in degrees or radians).
    """
    # Euclidean distance
    euclidean_distance = np.sqrt(horizontal_distance**2 + height**2)

    # Angle between the vertical and the hypotenuse (using arctan2)
    angle_from_vertical_rad = np.arctan2(horizontal_distance, height)

    # Grazing angle = 90° - angle_from_vertical
    grazing_angle_rad = (np.pi / 2) - angle_from_vertical_rad
    grazing_angle = np.degrees(grazing_angle_rad) if angle_in_degrees else grazing_angle_rad

    return euclidean_distance, grazing_angle


import numpy as np
from typing import Union

# --------------------------------------------------------------------
# Beaufort‑to‑wind‑speed lookup (m/s).  Values are the mid‑points
# of the standard wind‑speed ranges for each Beaufort number.
# --------------------------------------------------------------------
BEAUFORT_WIND_SPEEDS = np.array([
    0.15,   # 0 – <0.3  m/s   (no wind)
    1.0,   # 1 – 0.3–1.5 m/s
    2.4,   # 2 – 1.5–3.3 m/s
    4.2,   # 3 – 3.3–5.5 m/s
    6.2,   # 4 – 5.5–7.9 m/s
    9.0,   # 5 – 7.9–10.7 m/s
    12.2,  # 6 – 10.7–13.8 m/s
    15.8,  # 7 – 13.8–17.1 m/s
    19.4,  # 8 – 17.1–20.8 m/s
    23.5,  # 9 – 20.8–24.9 m/s
    28.0,  # 10 – 24.9–29.5 m/s
    33.5,  # 11 – 29.5–34.5 m/s
    35.0   # 12 – >34.5  m/s   (take 35 m/s as a typical upper limit)
])

def beaufort_to_wind(beaufort: Union[int, float, np.ndarray]) -> np.ndarray:
    """
    Convert a Beaufort number (or array of numbers) to a wind speed in m/s.
    Values outside the 0‑12 range are clipped to the nearest valid entry.
    """
    bf = np.asarray(beaufort).astype(int)
    # Clip to 0…12 so that indexing is safe
    bf = np.clip(bf, 0, BEAUFORT_WIND_SPEEDS.size - 1)
    return BEAUFORT_WIND_SPEEDS[bf]

# --------------------------------------------------------------------
# Chapman‑Harris noise model – now takes sea state instead of VVent
# --------------------------------------------------------------------
def surface_scattering_strength(
    angle_rasance: Union[float, np.ndarray, list, tuple],
    f: Union[np.ndarray, list, tuple],
    sea_state: Union[int, float, np.ndarray, list, tuple]
) -> np.ndarray:
    """
    Compute noise level (dB) from the Chapman‑Harris formula,
    but the wind speed is inferred from the sea state (Beaufort).

    Parameters
    ----------
    sea_state : array_like
        Sea‑state number(s) on the Beaufort scale (0–12).  Scalars or
        arrays of any shape are accepted.
    f : array_like
        Frequency vector (row).  Scalars are treated as 1‑element arrays.
    angle_rasance : array_like or scalar
        Rake angle(s) in degrees.  Scalars are broadcast to all combinations.

    Returns
    -------
    bruit : np.ndarray
        Noise level(s) in dB for each (sea_state, f) pair.
    """
    # 1. Convert sea state → wind speed
    VVent = beaufort_to_wind(sea_state)          # (n,) or (n,1)

    # 2. Ensure correct shapes for broadcasting
    VVent = np.asarray(VVent).reshape(-1, 1)      # (n,1)
    f = np.asarray(f).reshape(1, -1)              # (1,m)
    angle = np.asarray(angle_rasance)            # scalar or (p,)

    # 3. Compute Beta
    beta = 158.0 * (VVent * np.power(f, 1.0/3.0)) ** (-0.58)

    # 4. Compute Scatering
    s = (
        2.6
        - 42.4 * np.log10(beta)
        + 3.3 * beta * np.log10(np.abs(angle[:,np.newaxis]) / 30.0)
    )
    #(s.shape)
    return s


def bottom_scattering_strength(
    angle_deg: Union[float, np.ndarray],
    f: Union[float, np.ndarray],
    sediment_type: str = "Sand",
    flag_rasance: int = 1
) -> np.ndarray:
    """
    Mackinney‑Anderson loss (dB) for a given rake angle and frequency.

    Parameters
    ----------
    angle_deg : float or array_like
        Rake angle(s) in degrees.
    f : float or array_like
        Frequency(ies) in hertz (Hz).
    sediment_type : str, default "Sand"
        One of: "Rock", "Sand", "Mud".
    flag_rasance : int, default 1
        Flag selecting the formula branch (1 or 0).

    Returns
    -------
    pertes : np.ndarray
        Loss in dB, broadcasted to the shape of the inputs.
    """
    # ----------------------------------------------------------------------
    # Mapping from string → integer (original code used 1‑3)
    SEDIMENT_MAP = {
        "sand": 1,
        "mud": 2,
        "rock": 3
    }

    # ------------------------------------------------------------------
    # Convert the sediment string to the integer code used by the model
    # ------------------------------------------------------------------
    sed_key = sediment_type.strip().lower()
    if sed_key not in SEDIMENT_MAP:
        raise ValueError(
            f"Unknown sediment_type '{sediment_type}'. "
            f"Choose one of {list(SEDIMENT_MAP.keys())}."
        )
    type_fond = SEDIMENT_MAP[sed_key]

    # ------------------------------------------------------------------
    # Prepare arrays
    # ------------------------------------------------------------------
    angle_deg = np.asarray(angle_deg, dtype=float)
    f = np.asarray(f, dtype=float)

    f_khz = f / 1e3                     # convert to kHz

    # Helper trigonometric functions (degrees → radians)
    tan_deg = np.tan(np.radians(angle_deg))
    sin_deg = np.sin(np.radians(angle_deg))
    cos_deg = np.cos(np.radians(angle_deg))

    # ------------------------------------------------------------------
    # Main computation
    # ------------------------------------------------------------------
    if flag_rasance:
        # ----- Branch 1 (flag_rasance == 1) ----------------------------
        B = 1 + 125 * np.exp(
            -2.64 * (type_fond - 1.75) ** 2
            - 50 / (type_fond * tan_deg ** 2)
        )
        C = B * np.power(
            sin_deg + 0.19,
            (cos_deg ** 16)* f_khz[:,np.newaxis]
        )
        bs = 10 * np.log10(
            2.53 * C * f_khz[:,np.newaxis] ** (3.2 - 0.8 * type_fond)
            * 10 ** (2.8 * type_fond - 12)
            + 10 ** (-4.5)
        )
    else:
        # ----- Branch 2 (flag_rasance == 0) ----------------------------
        # Adjust the angle as in the original code
        angle_deg = 90.0 - angle_deg
        tan_deg = np.tan(np.radians(angle_deg))
        sin_deg = np.sin(np.radians(angle_deg))
        cos_deg = np.cos(np.radians(angle_deg))

        B = 1 + 125 * np.exp(
            -2.64 * (type_fond - 1.75) ** 2
            - (50 * tan_deg ** 2) / type_fond
        )
        C = B * np.power(
            cos_deg + 0.19,
            f_khz * (sin_deg ** 16)
        )
        bs = 10 * np.log10(
            2.53 * C * f_khz ** (3.2 - 0.8 * type_fond)
            * 10 ** (2.8 * type_fond - 12)
            + 10 ** (-4.5)
        )
    bs = bs.T
    #print(bs.shape)
    return bs

def volume_scattering_strength(
    frequency_hz,
    turbidity="Moderate",
    scatterer_density=1.0,
):
    """
    Estimate the volume scattering strength (in dB/m) based on frequency and water turbidity.

    Args:
        frequency_hz (float or np.ndarray): Frequency in Hz.
        turbidity (str): Turbidity level ("Clear", "Moderate", "Turbid", or "VeryTurbid").
                        Default: "Moderate".
        scatterer_density (float): Relative density of scatterers (1.0 = typical). Default: 1.0.

    Returns:
        float or np.ndarray: Volume scattering strength in dB/m.

    Raises:
        ValueError: If turbidity is not one of the supported levels.
    """
    # Define empirical models for each turbidity level
    turbidity_models = {
        "Clear": {
            # Low scattering due to minimal suspended particles
            "formula": lambda f: -90*np.ones(f.shape),# + 5 * np.log10(f) + 10 * np.log10(scatterer_density),
        },
        "Moderate": {
            # Typical scattering for coastal or open ocean waters
            "formula": lambda f: -70*np.ones(f.shape),# + 15 * np.log10(f) + 10 * np.log10(scatterer_density),
        },
        "Turbid": {
            # High scattering due to suspended sediments or high biological activity
            "formula": lambda f: -50*np.ones(f.shape),# + 25 * np.log10(f) + 10 * np.log10(scatterer_density),
        },
        "VeryTurbid": {
            # Very high scattering due to extreme turbidity (e.g., near river deltas or after storms)
            "formula": lambda f: -30*np.ones(f.shape),# + 35 * np.log10(f) + 10 * np.log10(scatterer_density),
        }
    }

    # Check if turbidity is valid
    if turbidity not in turbidity_models:
        raise ValueError("turbidity must be 'Clear', 'Moderate', 'Turbid', or 'VeryTurbid'.")

    # Get the formula for the specified turbidity level
    model = turbidity_models[turbidity]
    volume_scattering_db = model["formula"](frequency_hz)

    return volume_scattering_db


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
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
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
        config_layout = QGridLayout()

        self.spinboxes = {
            'Fmin': QDoubleSpinBox(),
            'Fmax': QDoubleSpinBox(),
            'SL': QDoubleSpinBox(),
            'SH': QDoubleSpinBox(),
            'Ne': QDoubleSpinBox(),
            'B': QDoubleSpinBox(),
            'DI': QDoubleSpinBox(),
            'T': QDoubleSpinBox(),
            'Dmin': QDoubleSpinBox(),
            'Dmax': QDoubleSpinBox(),
            'SeaState': QDoubleSpinBox(),
            'Depth': QDoubleSpinBox(),
            'Temp': QDoubleSpinBox(),
            'DT': QDoubleSpinBox(),
            'PFA': QDoubleSpinBox(),
            'Bottom': QDoubleSpinBox(),
            'BearingAperture': QDoubleSpinBox(),
            'ElevationAperture': QDoubleSpinBox(),
            'TS': QDoubleSpinBox(),
            'C': QDoubleSpinBox(),
            'NbCell': QDoubleSpinBox(),
        }

        self.spinboxes['Fmin'].setRange(100, 1000000)
        self.spinboxes['Fmin'].setSuffix(' Hz')
        self.spinboxes['Fmin'].setValue(self.args.Fmin)
        self.spinboxes['Fmin'].setToolTip("Minimum frequency of the sonar bandwidth (Hz). Controls the lowest frequency component of the sonar signal.")
        
        self.spinboxes['Fmax'].setRange(100, 1000000)
        self.spinboxes['Fmax'].setSuffix(' Hz')
        self.spinboxes['Fmax'].setValue(self.args.Fmax)
        self.spinboxes['Fmax'].setToolTip("Maximum frequency of the sonar bandwidth (Hz). Controls the highest frequency component of the sonar signal.")
        

        self.spinboxes['SL'].setRange(10, 230)
        self.spinboxes['SL'].setSuffix(' dB')
        self.spinboxes['SL'].setValue(self.args.SL)
        self.spinboxes['SL'].setToolTip("Source level (SL) of the sonar signal (dBuPa). Represents the acoustic power emitted by the sonar transducer.")
        

        self.spinboxes['SH'].setRange(-220, -150)
        self.spinboxes['SH'].setSuffix(' dBV/uPa')
        self.spinboxes['SH'].setValue(self.args.SH)
        self.spinboxes['SH'].setToolTip("Sensor sensitivity (SH) (dBV/uPa). Indicates the voltage output of the hydrophone for a given acoustic pressure input.")
        

        self.spinboxes['Ne'].setRange(1, 200)
        self.spinboxes['Ne'].setSuffix(' nV/sqrt(Hz)')
        self.spinboxes['Ne'].setValue(self.args.Ne)
        self.spinboxes['Ne'].setToolTip("Electronic noise (Ne) (nV/sqrt(Hz)). Represents the noise floor of the sonar receiver system.")
        

        self.spinboxes['B'].setRange(0.01, 50000)
        self.spinboxes['B'].setSuffix(' Hz')
        self.spinboxes['B'].setValue(self.args.B)
        self.spinboxes['B'].setToolTip("Bandwidth (B) of the sonar system (Hz). Determines the frequency range over which the sonar operates.")
        

        self.spinboxes['DI'].setRange(0, 50)
        self.spinboxes['DI'].setSuffix(' dB')
        self.spinboxes['DI'].setValue(self.args.DI)
        self.spinboxes['DI'].setToolTip("Detection integration time (DI) (dB). Represents the system's ability to integrate signals over time for detection.")
        

        self.spinboxes['T'].setRange(0.001, 50)
        self.spinboxes['T'].setSuffix(' s')
        self.spinboxes['T'].setValue(self.args.T)
        self.spinboxes['T'].setToolTip("Integration time (T) (s). The time over which the sonar integrates the received signal for processing.")
        

        self.spinboxes['Dmin'].setRange(1, 100000)
        self.spinboxes['Dmin'].setSuffix(' m')
        self.spinboxes['Dmin'].setValue(self.args.Dmin)
        self.spinboxes['Dmin'].setToolTip("Minimum operating range (Dmin) (m). The closest distance at which the sonar can effectively operate.")
        

        self.spinboxes['Dmax'].setRange(self.args.Dmin, 100000)
        self.spinboxes['Dmax'].setSuffix(' m')
        self.spinboxes['Dmax'].setValue(self.args.Dmax)
        self.spinboxes['Dmax'].setToolTip("Maximum operating range (Dmax) (m). The farthest distance at which the sonar can effectively operate.")
        

        self.spinboxes['DT'].setRange(0, 50)
        self.spinboxes['DT'].setSuffix(' dB')
        self.spinboxes['DT'].setValue(self.args.DT)
        self.spinboxes['DT'].setToolTip("Detection threshold (DT) (dB). The minimum signal level required for detection above the noise floor.")
        
        
        self.spinboxes['PFA'].setRange(-9, -2)
        self.spinboxes['PFA'].setSuffix(' 10**')
        self.spinboxes['PFA'].setValue(self.args.PFA)
        self.spinboxes['PFA'].setToolTip("Probability of false alarm (PFA) (10**). The likelihood of incorrectly detecting a signal when no actual target is present.")
        

        self.spinboxes['SeaState'].setRange(0, 6)
        self.spinboxes['SeaState'].setSuffix('')
        self.spinboxes['SeaState'].setValue(self.args.SeaState)
        self.spinboxes['SeaState'].setToolTip("Sea state (0-6). Indicates the wave height and sea conditions (0 = calm, 6 = very rough).")
        

        self.spinboxes['Depth'].setRange(0, 10000)
        self.spinboxes['Depth'].setSuffix(' m')
        self.spinboxes['Depth'].setValue(self.args.Depth)
        self.spinboxes['Depth'].setToolTip("Sonar operating depth (m). The depth at which the sonar is positioned in the water column.")
        

        self.spinboxes['Temp'].setRange(0, 60)
        self.spinboxes['Temp'].setSuffix(' degC')
        self.spinboxes['Temp'].setValue(self.args.Temp)
        self.spinboxes['Temp'].setToolTip("Water temperature (°C). Affects sound velocity and absorption in seawater.")
        

        self.spinboxes['BearingAperture'].setRange(1, 360)
        self.spinboxes['BearingAperture'].setSuffix(' deg')
        self.spinboxes['BearingAperture'].setValue(self.args.BearingAperture)
        self.spinboxes['BearingAperture'].setToolTip("Bearing aperture (deg). The angular width of the active sonar beam in the horizontal plane (0-360 degrees).")

        self.spinboxes['ElevationAperture'].setRange(1, 360)
        self.spinboxes['ElevationAperture'].setSuffix(' deg')
        self.spinboxes['ElevationAperture'].setValue(self.args.ElevationAperture)
        self.spinboxes['ElevationAperture'].setToolTip("Elevation aperture (deg). The angular width of the active sonar beam in the vertical plane (0-360 degrees).")
        

        self.spinboxes['TS'].setRange(-50, 50)
        self.spinboxes['TS'].setSuffix(' dB')
        self.spinboxes['TS'].setValue(self.args.TS)
        self.spinboxes['TS'].setToolTip("Target strength (TS) (dB). A measure of how much sound a target reflects back to the active sonar.")
        

        self.spinboxes['NbCell'].setRange(1, 100)
        self.spinboxes['NbCell'].setSuffix(' dB')
        self.spinboxes['NbCell'].setValue(self.args.NbCell)
        self.spinboxes['NbCell'].setToolTip("Number of cells integrated (NbCell). The number of accuracy celles integrated for detection.")
        

        self.spinboxes['Bottom'].setRange(0, 10000)
        self.spinboxes['Bottom'].setSuffix(' m')
        self.spinboxes['Bottom'].setValue(self.args.Bottom)
        self.spinboxes['Bottom'].setToolTip("Seabed depth (m). The depth of the ocean floor at the sonar location. Used for active sonar reverberation")
        

        self.spinboxes['C'].setRange(1300, 1600)
        self.spinboxes['C'].setSuffix(' m/s')
        self.spinboxes['C'].setValue(self.args.C)
        self.spinboxes['C'].setToolTip("Sound velocity (C) (m/s). The speed of sound in seawater, which varies with temperature, salinity, and pressure. Used for active sonar resolution")
        

        row = 0
        self.sonar_type = QComboBox()
        self.sonar_type.addItems(['Passive', 'Active'])
        self.sonar_type.setToolTip("Select sonar type: Passive (listens to ambient sounds) or Active (emits and receives sonar signals).")
        self.sonar_type.currentTextChanged.connect(self.update_plots)
        config_layout.addWidget(QLabel('SonarType'), row, 0)
        self.sonar_type.setCurrentText(self.args.SonarType)
        config_layout.addWidget(self.sonar_type, row, 1)

        row += 1
        self.combo_type = QComboBox()
        self.combo_type.addItems(['Incoherent', 'Coherent'])
        self.combo_type.setToolTip("Select processing type: Incoherent (sums power) or Coherent (sums complex signals).")
        self.combo_type.currentTextChanged.connect(self.update_plots)
        config_layout.addWidget(QLabel('ProcessingType'), row, 0)
        self.combo_type.setCurrentText(self.args.ProcessingType)
        config_layout.addWidget(self.combo_type, row, 1)

        for name, spinbox in self.spinboxes.items():
            row += 1
            config_layout.addWidget(QLabel(name), row, 0)
            config_layout.addWidget(spinbox, row, 1)
            spinbox.valueChanged.connect(self.update_plots)

        row += 1
        self.seabed_type = QComboBox()
        self.seabed_type.addItems(['Rock', 'Sand', 'Mud'])
        self.seabed_type.setToolTip("Select seabed type: Rock (hard surface), Sand (sediment), or Mud (soft sediment). Affects active sonar reverberation.")
        self.seabed_type.currentTextChanged.connect(self.update_plots)
        config_layout.addWidget(QLabel('SeabedType'), row, 0)
        self.seabed_type.setCurrentText(self.args.SeabedType)
        config_layout.addWidget(self.seabed_type, row, 1)

        row += 1
        self.turbidity = QComboBox()
        self.turbidity.addItems(["Clear", "Moderate", "Turbid", "VeryTurbid"])
        self.turbidity.currentTextChanged.connect(self.update_plots)
        self.turbidity.setToolTip("Select water turbidity level: Clear (low particles) to VeryTurbid (high particles). Affects active sonar signal attenuation.")
        config_layout.addWidget(QLabel('Turbidity'), row, 0)
        self.turbidity.setCurrentText(self.args.Turbidity)
        config_layout.addWidget(self.turbidity, row, 1)

        config_widget.setLayout(config_layout)


        # Ajouter une barre de défilement verticale
        scroll_area = QScrollArea()
        scroll_area.setWidget(config_widget)
        scroll_area.setWidgetResizable(True)
        dock_config.setWidget(scroll_area)
        scroll_area.setMaximumWidth(300)

        figsize=(6, 6)

        # Dock pour les plots ES
        dock_es = QDockWidget('Plots ES', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_es)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout()

        self.figure_es = Figure(figsize=figsize)
        self.canvas_es = FigureCanvas(self.figure_es)
        self.canvas_es.mpl_connect("motion_notify_event", self.on_mouse_move)

        plot_layout.addWidget(self.canvas_es)

        plot_widget.setLayout(plot_layout)
        dock_es.setWidget(plot_widget)
        
        # Dock pour les plots POD
        dock_pod = QDockWidget('Plots POD', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_pod)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout()

        self.figure_pod = Figure(figsize=figsize)
        self.canvas_pod = FigureCanvas(self.figure_pod)
        self.canvas_pod.mpl_connect("motion_notify_event", self.on_mouse_move)
        plot_layout.addWidget(self.canvas_pod)

        plot_widget.setLayout(plot_layout)
        dock_pod.setWidget(plot_widget)

        # Dock pour le plot CP
        dock_cp = QDockWidget('Plot CP', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_cp)

        cp_widget = QWidget()
        cp_layout = QVBoxLayout()

        self.figure_cp = Figure(figsize=figsize)
        self.canvas_cp = FigureCanvas(self.figure_cp)
        self.canvas_cp.mpl_connect("motion_notify_event", self.on_mouse_move)
        cp_layout.addWidget(self.canvas_cp)

        cp_widget.setLayout(cp_layout)
        dock_cp.setWidget(cp_widget)

        # Dock pour le plot Noises
        dock_noise = QDockWidget('Plot Noises', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_noise)

        noise_widget = QWidget()
        noise_layout = QVBoxLayout()

        self.figure_noise = Figure(figsize=figsize)
        self.canvas_noise = FigureCanvas(self.figure_noise)
        self.canvas_noise.mpl_connect("motion_notify_event", self.on_mouse_move)
        noise_layout.addWidget(self.canvas_noise)

        noise_widget.setLayout(noise_layout)
        dock_noise.setWidget(noise_widget)

        # Dock pour le plot SH
        dock_v = QDockWidget('Plot Tension', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_v)

        sh_widget = QWidget()
        sh_layout = QVBoxLayout()

        self.figure_volt = Figure(figsize=figsize)
        self.canvas_volt = FigureCanvas(self.figure_volt)
        self.canvas_volt.mpl_connect("motion_notify_event", self.on_mouse_move)
        sh_layout.addWidget(self.canvas_volt)

        sh_widget.setLayout(sh_layout)
        dock_v.setWidget(sh_widget)


        # Dock pour le bar graph de N en fonction de f
        dock_bargraph_N = QDockWidget('Noise', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_bargraph_N)

        bargraph_N_widget = QWidget()
        bargraph_N_layout = QVBoxLayout()

        self.figure_bargraph_N = Figure(figsize=figsize)
        self.canvas_bargraph_N = FigureCanvas(self.figure_bargraph_N)
        bargraph_N_layout.addWidget(self.canvas_bargraph_N)

        bargraph_N_widget.setLayout(bargraph_N_layout)
        dock_bargraph_N.setWidget(bargraph_N_widget)

        # Dock pour le bar graph de N en fonction de f
        dock_bargraph_A = QDockWidget('Absorption', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_bargraph_A)

        bargraph_A_widget = QWidget()
        bargraph_A_layout = QVBoxLayout()

        self.figure_bargraph_A= Figure(figsize=figsize)
        self.canvas_bargraph_A = FigureCanvas(self.figure_bargraph_A)
        bargraph_A_layout.addWidget(self.canvas_bargraph_A)

        bargraph_A_widget.setLayout(bargraph_A_layout)
        dock_bargraph_A.setWidget(bargraph_A_widget)

        # Dock pour le plot Scattering en fonction de l'angle
        dock_scattering_surface = QDockWidget('Scattering Surface', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_scattering_surface)

        scattering_surface_widget = QWidget()
        scattering_surface_layout = QVBoxLayout()

        self.figure_scattering_surface = Figure(figsize=figsize)
        self.canvas_scattering_surface = FigureCanvas(self.figure_scattering_surface)
        self.canvas_scattering_surface.mpl_connect("motion_notify_event", self.on_mouse_move)
        scattering_surface_layout.addWidget(self.canvas_scattering_surface)

        scattering_surface_widget.setLayout(scattering_surface_layout)
        dock_scattering_surface.setWidget(scattering_surface_widget)

        # Dock pour le plot Scattering en fonction de l'angle
        dock_scattering_bottom = QDockWidget('Scattering bottom', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_scattering_bottom)

        scattering_bottom_widget = QWidget()
        scattering_bottom_layout = QVBoxLayout()

        self.figure_scattering_bottom = Figure(figsize=figsize)
        self.canvas_scattering_bottom = FigureCanvas(self.figure_scattering_bottom)
        self.canvas_scattering_bottom.mpl_connect("motion_notify_event", self.on_mouse_move)
        scattering_bottom_layout.addWidget(self.canvas_scattering_bottom)

        scattering_bottom_widget.setLayout(scattering_bottom_layout)
        dock_scattering_bottom.setWidget(scattering_bottom_widget)


        # Dock pour le bar graph de N en fonction de f
        dock_bargraph_VS = QDockWidget('Volume Scattering', self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_bargraph_VS)

        bargraph_VS_widget = QWidget()
        bargraph_VS_layout = QVBoxLayout()

        self.figure_bargraph_VS= Figure(figsize=figsize)
        self.canvas_bargraph_VS = FigureCanvas(self.figure_bargraph_VS)
        bargraph_VS_layout.addWidget(self.canvas_bargraph_VS)

        bargraph_VS_widget.setLayout(bargraph_VS_layout)
        dock_bargraph_VS.setWidget(bargraph_VS_widget)

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

        <p>For more information, on the full featured Software, visit the <a href="http://www.seagnal.fr">SORAPS product page</a>.</p>
        
        <p>or write us a <a href="mailto:contact@seagnal.fr?subject=SORAPS">mail</a></p>
        
        <p>Author:Johann Baudy.</p>
        <p>Company:<a href="http://www.seagnal.fr">SEAGNAL</a>.</p>
        """)

        dock_about.setWidget(text_browser)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_about)

        self.tabifyDockWidget(dock_about, dock_bargraph_N)
        self.tabifyDockWidget(dock_bargraph_N, dock_bargraph_A)
        self.tabifyDockWidget(dock_bargraph_A, dock_bargraph_VS)
        self.tabifyDockWidget(dock_bargraph_VS, dock_scattering_surface)
        self.tabifyDockWidget(dock_scattering_surface, dock_scattering_bottom)
        self.tabifyDockWidget(dock_scattering_bottom, dock_v)
        self.tabifyDockWidget(dock_v, dock_noise)
        self.tabifyDockWidget(dock_noise, dock_cp)
        self.tabifyDockWidget(dock_cp, dock_pod)
        self.tabifyDockWidget(dock_pod, dock_es)

        self.dock_bargraph_N = dock_bargraph_N
        self.dock_bargraph_A = dock_bargraph_A
        self.dock_bargraph_VS = dock_bargraph_VS
        self.dock_scattering_surface = dock_scattering_surface
        self.dock_scattering_bottom = dock_scattering_bottom
        self.dock_v = dock_v
        self.dock_noise = dock_noise
        self.dock_cp = dock_cp
        self.dock_pod = dock_pod
        self.dock_es = dock_es

        self.init = True
        self.update_plots()

        if self.args.save:
            self.save_plot()
            sys.exit(0)


    def on_mouse_move(self, event):
        if event.inaxes:
            # Get the mouse coordinates in data units
            x, y = event.xdata, event.ydata
            print(f"x: {x:.2f}, y: {y:.2f}")
            # Update the status bar
            self.status_bar.showMessage(f"x: {x:.2f}, y: {y:.2f}")

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
        Ne = 20*np.log10(self.spinboxes['Ne'].value()*1e-9)
        Depth = self.spinboxes['Depth'].value()
        Temp = self.spinboxes['Temp'].value()
        Pfa = self.spinboxes['PFA'].value()
        SeaState = self.spinboxes['SeaState'].value()
        BearingAperture = self.spinboxes['BearingAperture'].value()
        ElevationAperture = self.spinboxes['ElevationAperture'].value()
        ProcessingType = self.combo_type.currentText()
        SonarType = self.sonar_type.currentText()
        Turbidity = self.turbidity.currentText()
        SeaBedType = self.seabed_type.currentText()
        C = self.spinboxes['C'].value()
        Bottom = self.spinboxes['Bottom'].value()
        TS = self.spinboxes['TS'].value()
        K = self.spinboxes['NbCell'].value()

        if SonarType == 'Active':
            is_active = True
            self.spinboxes['NbCell'].setEnabled(True)
        else:
            is_active = False
            if ProcessingType == "Coherent":
                self.spinboxes['NbCell'].setEnabled(True)
            else:
                self.spinboxes['NbCell'].setEnabled(False)

        self.seabed_type.setEnabled(is_active)
        self.turbidity.setEnabled(is_active)
        self.spinboxes['C'].setEnabled(is_active)
        self.spinboxes['TS'].setEnabled(is_active)
        self.spinboxes['Bottom'].setEnabled(is_active)
        self.spinboxes['BearingAperture'].setEnabled(is_active)
        self.spinboxes['ElevationAperture'].setEnabled(is_active)


        d = np.arange(Dmin, Dmax + 1, 1)
        f = np.transpose(np.array([Fmin, (Fmin + Fmax) / 2, Fmax]))
        np_erf = np.vectorize(math.erf)
        Pfa_ = math.pow(10, Pfa)


        # Calculate the noise level (Na) based on the sea state and frequency
        # SeaState: Represents the condition of the sea, affecting ambient noise levels
        # f: Frequency at which the noise level is calculated
        Na = calculate_nl(SeaState, f)

        # Calculate the absorption coefficient (alpha) using the Fisher-Simmons model
        # This model estimates sound absorption in seawater based on frequency, depth, and temperature
        # f: Frequency of the sound wave
        # D: Depth of the water, set to the variable Depth
        # T: Temperature of the water, set to the variable Temp
        alpha = absorption_fisher_simmons(f, D=Depth, T=Temp)

        # CP: Computation of    
        # This line calculates the propagation loss (CP) which includes spherical spreading loss and absorption loss.
        # 20 * np.log10(d[:, np.newaxis]) represents spherical spreading loss.
        # alpha * d[:, np.newaxis] represents absorption loss, where alpha is the absorption coefficient.
        

        #print(SonarType)
        if is_active:
            CP = 2.0*(20.0 * np.log10(d[:, np.newaxis]) + alpha * d[:, np.newaxis])

            Resolution = C/(2*B)
            distance_surface, grazing_angle_deg_suface = distance_and_grazing_angle(d, Depth)
            distance_seabed, grazing_angle_deg_seabed = distance_and_grazing_angle(d, Bottom-Depth)

            A_surface = Resolution*circular_arc_length(
                d, 
                BearingAperture)

            A_seabed = Resolution*circular_arc_length(
                d, 
                BearingAperture)

            V = Resolution*spherical_cap_area_vectorized(
                d, 
                BearingAperture, 
                ElevationAperture)

            volume_scattering = volume_scattering_strength(f,turbidity=Turbidity)
            seabed_scattering = bottom_scattering_strength(
                grazing_angle_deg_seabed,
                f,
                sediment_type=SeaBedType,
            )
            surface_scattering = surface_scattering_strength(
                grazing_angle_deg_suface,
                f, 
                sea_state=SeaState)
            print('V',volume_scattering)
            print('SS',surface_scattering)
            print('BS',seabed_scattering)
            Nr_seabed = SL - CP + 10*np.log10(A_seabed[:,np.newaxis])+(seabed_scattering)
            Nr_seabed[d<(Bottom-Depth),:] = 0
            Nr_surface =  SL - CP + 10*np.log10(A_surface[:,np.newaxis])+(surface_scattering)
            Nr_surface[d<Depth,:] = 0
            Nv =  SL - CP + 10*np.log10(V[:,np.newaxis])+(volume_scattering)
            print('V',Nv)
            print('SS',Nr_surface)
            print('BS',Nr_seabed)
            
           
            #Nt = 10*np.log10(np.pow(10.0,(Ne - SH - 10.0 * np.log10(T))/10) + np.pow(10.0,(Na - DI - 10.0 * np.log10(T))/10) + np.pow(10.0,(Nr_seabed)/10) + np.pow(10.0,(Nr_surface)/10) + np.pow(10.0,(Nv)/10))- 10*np.log10((2*K)/(K+1))
            Nte = Ne - SH - DI - 10.0 * np.log10(T) - 10*np.log10((2*K)/(K+1))
            Nte = Nte*np.ones(CP.shape)

            Nta = Na - DI - 10.0 * np.log10(T) - 10*np.log10((2*K)/(K+1))
            Nta = Nta[np.newaxis,:]*np.ones(CP.shape)

            Ntrv = Nv - 10*np.log10((2*K)/(K+1))
            Ntrss = Nr_surface - 10*np.log10((2*K)/(K+1))
            Ntrbs = Nr_seabed - 10*np.log10((2*K)/(K+1))

            Nt = 10.0*np.log10(np.pow(10.0,Nta/10.0) + np.pow(10.0,Nte/10.0) + np.pow(10.0,Ntrv/10.0) + np.pow(10.0,Ntrss/10.0) + np.pow(10.0,Ntrbs/10.0))
            print(Nt[0,0])
            print(Nta[0,0], Na, DI)
            print(Nte[0,0], Ne, SH)
            print(Ntrv[0,0])
            print(Ntrss[0,0])
            print(Ntrbs[0,0])
            
           

            if not ProcessingType == "Coherent":
                print("Incoherent processing is not supported")


        else:
            CP = 20 * np.log10(d[:, np.newaxis]) + alpha * d[:, np.newaxis]
            
            if ProcessingType == "Coherent":

                Nta = Na - DI + 10.0 * np.log10(B) - 10.0 * np.log10(T*B) - 10*np.log10((2*K)/(K+1))
                Nte = Ne - SH + 10.0 * np.log10(B) - 10.0 * np.log10(T*B) - 10*np.log10((2*K)/(K+1))
        
            else:
                K = B*T
                Nta = Na - DI + 10.0 * np.log10(B) -                    0 - 10*np.log10((2*K)/(K+1))
                Nte = Ne - SH + 10.0 * np.log10(B) -                    0 - 10*np.log10((2*K)/(K+1))
        
            Nt = 10*np.log10(pow(10.0,Nta/10) + np.pow(10.0,Nte/10))

        # V: Voltage or signal level after propagation
        # This line calculates the voltage or signal level (V) at the receiver.
        # SL is the Source Level, CP is the Propagation Loss, and SH is the Signal Excess or Hydrophone Sensitivity.
        V = np.pow(10,(SL - CP + SH)/20)

        # ES: Signal Excess for Coherent and Incoherent Processing
        # Depending on the type of processing (Coherent or Incoherent), the Signal Excess (ES) is calculated differently.
        E_N_post = 10.**((Nt)/10);
        Gain=(2*K)/(K+1)
        E_N = E_N_post*Gain;
        E_Smin=np.pow(10,(DT/10))*E_N_post;
        
        # Calculate detection threshold (dt_s) based on the noise level and probability of false alarm of quadratic detector
        if K == 1:
            dt_s= E_N *(-np.log(Pfa_));
        else:
            dt_s=(np.sqrt(2)*E_N*inverse_erf(1-2*Pfa_)+E_N*np.sqrt(K))/np.sqrt(K);

        E_S = np.pow(10,(SL - CP + TS)/10);
        
        DT_db=(10*np.log10(E_Smin/E_N_post));
        
        SNR_db =(10*np.log10(E_S/E_N_post));

        ES_db = SNR_db - DT_db;
        
        # Probability of detection
        if K == 1:
            S_s = np.sqrt(E_N);
            POD_o=1-(np_erf((np.sqrt(dt_s*dt_s)+E_S)/(np.sqrt(2)*S_s*S_s))+np_erf((np.sqrt(dt_s*dt_s)-E_S)/(np.sqrt(2)*S_s*S_s)))/2;
        else:
            S_n = E_N;
            POD_o=1/2-np_erf((np.sqrt(K)*(dt_s-E_S))/(np.sqrt(2)*S_n))/2;
        
        # Plot POD
        self.figure_pod.clear()
        ax_pod = self.figure_pod.add_subplot(111)
        cursor_pod = Cursor(ax_pod, color='black', linewidth=1)
        ax_pod.plot(d, POD_o[:, 0], label='Fmin', color='blue')
        ax_pod.plot(d, POD_o[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax_pod.plot(d, POD_o[:, 2], label='Fmax', color='red')
        ax_pod.grid()
        ax_pod.legend()
        ax_pod.set_xlabel('Distance (m)')
        ax_pod.set_ylabel('POD (dB)')
        ax_pod.set_title(f'POD {self.combo_type.currentText()}')
        ax_pod.set_ylim(-0.1, 1.1)
        self.canvas_pod.draw()

        # Plot ES
        self.figure_es.clear()
        ax_es = self.figure_es.add_subplot(111)
        cursor_es = Cursor(ax_es, color='black', linewidth=1)
        ax_es.plot(d, ES_db[:, 0], label='Fmin', color='blue')
        ax_es.plot(d, ES_db[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax_es.plot(d, ES_db[:, 2], label='Fmax', color='red')
        ax_es.grid()
        ax_es.legend()
        ax_es.set_xlabel('Distance (m)')
        ax_es.set_ylabel('ES (dB)')
        ax_es.set_title(f'ES {self.combo_type.currentText()}')
        self.canvas_es.draw()

        # Plot CP
        self.figure_cp.clear()
        ax_cp = self.figure_cp.add_subplot(111)
        cursor_cp = Cursor(ax_cp, color='black', linewidth=1)
        ax_cp.plot(d, CP[:, 0], label='Fmin', color='blue')
        ax_cp.plot(d, CP[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax_cp.plot(d, CP[:, 2], label='Fmax', color='red')
        ax_cp.legend()
        ax_cp.grid()
        ax_cp.set_xlabel('Distance (m)')
        ax_cp.set_ylabel('Field of loss (dB)')
        ax_cp.set_title(f'Field of loss')
        self.canvas_cp.draw()

        # Plot N
        #print(Nt.shape)
        #print(Nta.shape)
        #Ntap = Nta[np.newaxis, :]*np.ones(Nt.shape)
        #Ntep = Nte*np.ones(Nt.shape)
        #print(Ntap.shape)
        #print(Ntrv.shape)
        #print(Ntrss.shape)
        #print(Ntrbs.shape)
        self.figure_noise.clear()
        ax_noise = self.figure_noise.add_subplot(111)

        if is_active:
            #cursor_noise = Cursor(ax_noise, color='black', linewidth=2)
            ax_noise.plot(d, Nt[:, 0], label='Total Fmin', color='blue')
            ax_noise.plot(d, Nt[:, 1], label='Total Fmoy', color='blue', linestyle='dotted')
            ax_noise.plot(d, Nt[:, 2], label='Total Fmax', color='blue', linestyle='dashdot')

            ax_noise.plot(d, Nte[:, 1], label='Electronic', color='gray')

            ax_noise.plot(d, Nta[:, 0], label='Ambiant Fmin', color='red')
            ax_noise.plot(d, Nta[:, 1], label='Ambiant Fmoy', color='red', linestyle='dotted')
            ax_noise.plot(d, Nta[:, 2], label='Ambiant Fmax', color='red', linestyle='dashdot')


            ax_noise.plot(d, Ntrv[:, 0], label='Volume Fmin', color='green')
            ax_noise.plot(d, Ntrv[:, 1], label='Volume Fmoy', color='green', linestyle='dotted')
            ax_noise.plot(d, Ntrv[:, 2], label='Volume Fmax', color='green', linestyle='dashdot')

            ax_noise.plot(d, Ntrss[:, 0], label='Surface Fmin', color='orange')
            ax_noise.plot(d, Ntrss[:, 1], label='Surface Fmoy', color='orange', linestyle='dotted')
            ax_noise.plot(d, Ntrss[:, 2], label='Surface Fmax', color='orange', linestyle='dashdot')

            ax_noise.plot(d, Ntrbs[:, 0], label='Bottom Fmin', color='purple')
            ax_noise.plot(d, Ntrbs[:, 1], label='Bottom Fmoy', color='purple', linestyle='dotted')
            ax_noise.plot(d, Ntrbs[:, 2], label='Bottom Fmax', color='purple', linestyle='dashdot')

            ax_noise.legend()
            ax_noise.grid()
            ax_noise.set_xlabel('Distance (m)')
            ax_noise.set_ylabel('Noise (dBuPa)')
            ax_noise.set_title(f'Noise')
        else:
            ax_noise.bar(['Total Fmin','Total Fmoy','Total Fmax', 'Ambiant Fmin','Ambiant Fmoy','Ambiant Fmax', 'Electronic'], [Nt[0], Nt[1], Nt[2], Nta[0], Nta[1], Nta[2], Nte], color=['blue','blue','blue', 'red', 'red', 'red', 'gray'])
            #print(Nta.shape)
            #print(Nt.shape)
            #print(Nte.shape)
            #ax_noise.bar(['Fmin','Fmoy','Fmax'], Nta, color='red')
            ax_noise.set_xlabel('Frequency (Hz)')
            ax_noise.set_ylabel('Noise Total (dBuPa/sqrt(Hz))')
            ax_noise.set_title(f'Total Noise over Bandwidth - SeaState:{SeaState} - Nte:{Nte}')
            ax_noise.grid(axis='y')

        self.canvas_noise.draw()

        # Plot SH
        self.figure_volt.clear()
        ax_volt = self.figure_volt.add_subplot(111)
        #cursor_volt = Cursor(ax_volt, color='black', linewidth=1)
        ax_volt.semilogy(d, V[:,0], label='Fmin', color='blue')
        ax_volt.semilogy(d, V[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax_volt.semilogy(d, V[:, 2], label='Fmax', color='red')
        ax_volt.legend()
        ax_volt.grid()
        ax_volt.set_xlabel('Distance (m)')
        ax_volt.set_ylabel('Vrms')
        ax_volt.set_title(f'Voltage at sensor')
        self.canvas_volt.draw()

        # Bar Graph N vs f
        self.figure_bargraph_N.clear()
        ax_bargraph_N = self.figure_bargraph_N.add_subplot(111)
        #cursor_bargraph_N = Cursor(ax_bargraph_N, color='black', linewidth=2)
        ax_bargraph_N.bar(['Fmin','Fmoy','Fmax'], Na, color=['blue','gray','red'])
        ax_bargraph_N.set_xlabel('Frequency (Hz)')
        ax_bargraph_N.set_ylabel('N (dBuPa/sqrt(Hz))')
        ax_bargraph_N.set_title(f'Ambiant Noise over Bandwidth - SeaState:{SeaState}')

        ax_bargraph_N.grid(axis='y')
        self.canvas_bargraph_N.draw()

        VS = volume_scattering_strength(f,turbidity=Turbidity)
        self.figure_bargraph_VS.clear()
        ax_bargraph_VS = self.figure_bargraph_VS.add_subplot(111)
        #cursor_bargraph_VS = Cursor(ax_bargraph_VS, color='black', linewidth=2)
        ax_bargraph_VS.bar(['Fmin','Fmoy','Fmax'], VS, color=['blue','gray','red'])
        ax_bargraph_VS.set_xlabel('Frequency (Hz)')
        ax_bargraph_VS.set_ylabel('(dB)')
        ax_bargraph_VS.set_title(f'Volume Scattering over Bandwidth - SeaType:{SeaBedType}')

        ax_bargraph_VS.grid(axis='y')
        self.canvas_bargraph_VS.draw()

        # Bar Graph N vs f
        self.figure_bargraph_A.clear()
        ax_bargraph_A = self.figure_bargraph_A.add_subplot(111)
        #cursor_bargraph_A = Cursor(ax_bargraph_A, color='black', linewidth=2)
        ax_bargraph_A.bar(['Fmin','Fmoy','Fmax'], alpha, color=['blue','gray','red'])
        ax_bargraph_A.set_xlabel('Frequency (Hz)')
        ax_bargraph_A.set_ylabel('absorption db/m')
        ax_bargraph_A.set_title(f'Absorption over Bandwidth - T:{Temp} - D:{Depth}')
        ax_bargraph_A.grid(axis='y')
        self.canvas_bargraph_A.draw()

        # Plot Scattering
        angles_deg  = np.linspace(0, 90, 90)
        SS = surface_scattering_strength(angles_deg,f,SeaState)
        
        self.figure_scattering_surface.clear()
        ax_scattering_surface = self.figure_scattering_surface.add_subplot(111)
        #cursor_scattering_surface = Cursor(ax_scattering_surface, color='black', linewidth=2)
        ax_scattering_surface.plot(angles_deg, SS[:, 0], label='Fmin', color='blue')
        ax_scattering_surface.plot(angles_deg, SS[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax_scattering_surface.plot(angles_deg, SS[:, 2], label='Fmax', color='red')
        ax_scattering_surface.legend()
        ax_scattering_surface.grid()
        ax_scattering_surface.set_xlabel('Angle (deg)')
        ax_scattering_surface.set_ylabel('dB')
        ax_scattering_surface.set_title(f'Surface scattering')
        self.canvas_scattering_surface.draw()

        BS = bottom_scattering_strength(angles_deg,f,SeaBedType)
        #print(BS)
        self.figure_scattering_bottom.clear()
        ax_scattering_bottom = self.figure_scattering_bottom.add_subplot(111)
        #cursor_scattering_bottom = Cursor(ax_scattering_bottom, color='black', linewidth=2)
        ax_scattering_bottom.plot(angles_deg, BS[:, 0], label='Fmin', color='blue')
        ax_scattering_bottom.plot(angles_deg, BS[:, 1], label='Fmoy', color='gray', linestyle='dotted')
        ax_scattering_bottom.plot(angles_deg, BS[:, 2], label='Fmax', color='red')
        ax_scattering_bottom.legend()
        ax_scattering_bottom.grid()
        ax_scattering_bottom.set_xlabel('Angle (deg)')
        ax_scattering_bottom.set_ylabel('dB')
        ax_scattering_bottom.set_title(f'Bottom scattering')
        self.canvas_scattering_bottom.draw()

        #print('--- END UPDATE PLOTS')

    def save_plot(self):
        print('--- START SAVE PLOTS')
        dpi = 150

        #figure_size = (6, 4)  # Width: 6 inches, Height: 4 inches
        self.dock_pod.show()
        self.dock_pod.raise_()
        self.figure_pod.savefig('plot_pod.png', dpi=dpi)
        
        self.dock_es.show()
        self.dock_es.raise_()
        self.figure_es.savefig('plot_es.png', dpi=dpi)

        self.dock_cp.show()
        self.dock_cp.raise_()
        self.figure_cp.savefig('plot_cp.png', dpi=dpi)

        self.dock_v.show()
        self.dock_v.raise_()
        self.figure_volt.savefig('plot_voltage.png', dpi=dpi)

        self.dock_noise.show()
        self.dock_noise.raise_()
        self.figure_noise.savefig('plot_noise.png', dpi=dpi)

        self.dock_scattering_bottom.show()
        self.dock_scattering_bottom.raise_()
        self.figure_scattering_bottom.savefig('plot_scattering_bottom.png', dpi=dpi)

        self.dock_scattering_surface.show()
        self.dock_scattering_surface.raise_()
        self.figure_scattering_surface.savefig('plot_scattering_surface.png', dpi=dpi)

        self.dock_bargraph_N.show()
        self.dock_bargraph_N.raise_()
        self.figure_bargraph_N.savefig('plot_bargraph_N.png', dpi=dpi)

        self.dock_bargraph_A.show()
        self.dock_bargraph_A.raise_()
        self.figure_bargraph_A.savefig('plot_bargraph_A.png', dpi=dpi)

        self.dock_bargraph_VS.show()
        self.dock_bargraph_VS.raise_()
        self.figure_bargraph_VS.savefig('plot_bargraph_VS.png', dpi=dpi)

        # Open a file for writing
        with open('variables.txt', 'w') as file:
            file.write(f"ProcessingType: {self.combo_type.currentText()}\n")
            file.write(f"SonarType: {self.sonar_type.currentText()}\n")
            file.write(f"Turbidity: {self.turbidity.currentText()}\n")
            file.write(f"SeaBedType: {self.seabed_type.currentText()}\n")

            for name, spinbox in self.spinboxes.items():
                # Extract the unit from the spinbox suffix
                unit = spinbox.suffix().strip()
                value = spinbox.value()
                file.write(f"{name}: {value} {unit}\n")

        print('--- END UPDATE PLOTS')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Mini-Soraps')
    parser.add_argument('--SonarType', type=str, choices=["Active", "Passive"], default='Active', help='Type of sonar')
    parser.add_argument('--ProcessingType', type=str, choices=['Incoherent', 'Coherent'], default='Coherent', help='Type of treatment')
    parser.add_argument('--Fmin', type=float, default=10000, help='Fmin value (Hz)')
    parser.add_argument('--Fmax', type=float, default=20000, help='Fmax value (Hz)')
    parser.add_argument('--SL', type=float, default=150, help='SL value (dBuPa)')
    parser.add_argument('--SH', type=float, default=-180, help='SH value (dBV/uPa)')
    parser.add_argument('--Ne', type=float, default=10, help='Electronic noise (nV/sqrt(Hz))')
    parser.add_argument('--B', type=float, default=500, help='B value (Hz)')
    parser.add_argument('--DI', type=float, default=0, help='DI value (dB)')
    parser.add_argument('--T', type=float, default=0.01, help='T value (s)')
    parser.add_argument('--Dmin', type=float, default=10, help='Dmin value (m)')
    parser.add_argument('--Dmax', type=float, default=3000, help='Dmax value (m)')
    parser.add_argument('--DT', type=float, default=0, help='DT value (dB)')
    parser.add_argument('--PFA', type=float, default=-7, help='PFA value (10**)')
    parser.add_argument('--NbCell', type=float, default=1, help='Number of cells intergrated')
    parser.add_argument('--SeaState', type=int, default=1, help='SeaState value (0-6)')
    parser.add_argument('--Temp', type=float, default=15, help='Temp value (degC)')
    parser.add_argument('--Depth', type=float, default=100, help='Depth of sonar (m)')
    parser.add_argument('--C', type=float, default=1500, help='Sound Velocity (m/s) - Active only')
    parser.add_argument('--Bottom', type=float, default=1000, help='Seabed Depth value (m) - Active only')
    parser.add_argument('--BearingAperture', type=float, default=360, help='Bearing Aperture (deg) - Active only')
    parser.add_argument('--ElevationAperture', type=float, default=180, help='Elevation Aperture (deg) - Active only')
    parser.add_argument('--TS', type=float, default=0, help='Target Strength value (dB) - Active only')
    parser.add_argument('--Turbidity', type=str, choices=["Clear", "Moderate", "Turbid", "VeryTurbid"], default='Clear', help='Type of tubidity - Active only')
    parser.add_argument('--SeabedType', type=str, choices=["Rock", "Sand", "Mud"], default='Rock', help='Type of seabed - Active only')
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
