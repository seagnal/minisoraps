# minisoraps
mini-SORAPS is a free, lightweight version of SORAPS for public use. It simulates acoustic wave propagation in a 1D marine environment, aiding predictions for sonar and underwater communication systems.

# Installation Procedure

## Prerequisites

Before you begin, ensure you have the following installed on your system:

- Python (3.6 or higher)
- pip (Python package installer)

## Step 1: Set Up a Virtual Environment

It is recommended to use a virtual environment to manage dependencies.

```bash
python -m venv myenv
source myenv/bin/activate  # On Windows use `myenv\Scripts\activate`
```

## Step 2: Install Required Packages

Install the necessary Python packages using pip:

```bash
pip install numpy PyQt5 matplotlib
```

## Troubleshooting

- If you encounter any issues during installation, ensure that your pip version is up-to-date by running `pip install --upgrade pip`.
- Make sure you have the necessary permissions to install packages on your system.
- If you are behind a proxy, configure pip to use the proxy settings.

# Running The Application

Ensure your main application script includes the necessary imports and run it using Python:

```bash
python mini-soraps.py
```
or with arguments:
```bash
python mini-soraps.py [-h] [--Fmin FMIN] [--Fmax FMAX] [--SL SL] [--SH SH] [--B B] [--DI DI] [--T T] [--Dmin DMIN] [--Dmax DMAX] [--DT DT] [--PFA PFA] [--SeaState SEASTATE] [--Temp TEMP] [--Depth DEPTH] [--TypeTraitement {Incoherent,Coherent}] [--save]

```

# Additional Resources

- [SORAPS SEAGNAL Website] (http://products.seagnal.fr)

