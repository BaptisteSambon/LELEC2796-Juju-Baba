# LELEC2796-Juju-Baba

Welcome to the repository of  [Baptiste](https://github.com/BaptisteSambon)  & [Justin](https://github.com/Just1Wmls) for the project of the course [LELEC2796 - Wireless Communications](https://uclouvain.be/en-cours-2023-lelec2796) given at [UClouvain](https://uclouvain.be/fr/index.html) ! 

[![made-with-python](https://img.shields.io/badge/python-%2314354C.svg?&style=for-the-badge&logo=python&logoColor=white)](https://www.youtube.com/watch?v=dQw4w9WgXcQ)

[![made-with-sionna](https://img.shields.io/badge/Sionna_library-orange)](https://nvlabs.github.io/sionna/)

This Readme contains information about the organisation of the repository, the codes and how to run them.

**Keywords** : NOMA, Non Orthogonal Multiple Access, CDMA, Code Division Multiple Access, BER, SNR, Spectral Efficiency, Users Scalability, LELEC2796, Wireless Communications. 

## Introduction

This repository contains the codes for the project of the course [LELEC2796 - Wireless Communications](https://uclouvain.be/en-cours-2023-lelec2796) at [UClouvain](https://uclouvain.be/fr/index.html). The goal of this project is to compare the performance of NOMA and CDMA in a wireless communication system. The comparison is based on the BER, the spectral efficiency and the users scalability. The codes are written in Python and use the Sionna library. The results are plotted using Matplotlib.

## Table of contents

- [Introduction](#introduction)
- [Table of contents](#table-of-contents)
- [How to run](#how-to-run)
- [Scenario](#scenario)
- [Available codes](#available-codes)
- [Note](#note)
- [References](#references)
- [Authors](#authors)

## How to run

Clone the project

```bash
  git clone https://github.com/BaptisteSambon/LELEC2796-Juju-Baba.git 
```

Go to the project directory

```bash
  cd LELEC2796-Juju-Baba
```

### Requirements

To run the codes, you need to install the following libraries:

- [Sionna](https://nvlabs.github.io/sionna/)
- [Matplotlib](https://matplotlib.org/)
- [Numpy](https://numpy.org/)
- [Scipy](https://www.scipy.org/)
- [Pandas](https://pandas.pydata.org/)
- [Random](https://docs.python.org/3/library/random.html)
- [Plotly](https://plotly.com/python/)
- [Tensorflow](https://www.tensorflow.org/)
- [Termcolor](https://pypi.org/project/termcolor/)
- [Seaborn](https://seaborn.pydata.org/)

Sionna is a library developed by the [NVIDIA](https://www.nvidia.com/en-us/) team. It is a library for simulating wireless communication systems. It is based on the [Tensorflow](https://www.tensorflow.org/) library. It is a very powerful library but it is not easy to install. To install it, you need to follow the instructions on the [Sionna website](https://nvlabs.github.io/sionna/installation.html).


> [!TIP]\
> If an error occurs the first time you run the import cell in the notebooks, don't panic it's totally fine. Despite a "try expect" method, Sionna is not easy to import. Thus you often need to run the codes two times before seeing the plotted results.

> [!NOTE]\
> The codes were written in Python 3.8.5. They were tested on Windows 10. They should work on any operating system and any version of Python 3. However, we cannot guarantee that they will work on your computer. If you encounter any problem, please contact us.

> [!IMPORTANT]\
> The variables and parameters used for the simulations are defined in the first cell of the notebooks. You can change them to see the impact on the results. 

## Scenario

A specific scenario is defined with shared assumptions for fair evaluation. It involves downlink communication from one Base Station (BS) to several User Equipments (UEs) with fixed initial power, using Binary Phase-Shift Keying (BPSK) modulation. Perfect Channel State Information (CSI) is assumed, with a flat fading Line of Sight (LoS) configuration and non-zero Inter-Cell and Inter-User Interference. Minimum Mean Square Error equalization (MMSE) is applied at the receiver. 

The inter-cell interference (ICI) is generated using the [cells.ipynb](./cells.ipynb) notebook by considering a grid of cells. The number of cells, the density of cells, the obstacles, etc. can be changed. The ICI is then used in the [NOMA.ipynb](./NOMA.ipynb) and [CDMA.ipynb](./CDMA.ipynb) notebooks to simulate the NOMA and CDMA systems. The results are plotted in the [NOMAvsCDMA.ipynb](./NOMAvsCDMA.ipynb) notebook.

## Available codes

The codes in this repository are divided into 4 Jupyter notebooks and 3 python files. 

The notebooks are:

- #### [NOMA.ipynb](./NOMA.ipynb)
Contains code for NOMA and plots the BER, the spectral efficiency and the users scalability for different values of the SNR. 

- #### [CDMA.ipynb](./CDMA.ipynb)
Contains code for CDMA and plots the BER, the spectral efficiency and the users scalability for different values of the SNR.

- #### [NOMAvsCDMA.ipynb](./NOMAvsCDMA.ipynb)
Contains code for plots comparing the performance (BER, spectral efficiency, users scalability) of the NOMA and CDMA systems for different values of the SNR, different number of users. 

- #### [cells.ipynb](./cells.ipynb)
Contains the code used to generate the inter-cells interference with different scenarios investigated (density of cells, obstacles, etc.).

The python files are:

- #### [utils.py](./utils.py)
Contains all the functions used in the notebooks. Each function is documented.

- #### [cells.py](./cells.py)
Contains the class used to generate the points (coordinates) for the obstacles in [cells.ipynb](./cells.ipynb) and used in [obstacles.py](./obstacles.py).

- #### [obstacles.py](./obstacles.py)
Contains the class used to generate the circular obstacles in [cells.ipynb](./cells.ipynb).


### Note

Folders and files with the prefix "old" are not used in the notebooks. They were used for previous versions of the codes. They are kept for reference.

In the Folder [data](./data), you can find the data used to plot the results in the [NOMAvsCDMA.ipynb](./NOMAvsCDMA.ipynb) notebook and the interference arrays generated in the [cells.ipynb](./cells.ipynb) notebook (representing the values of interference obtained by simulations and with several iterations and random initialisations of the BSs coordinates).

## References

- . 


### Links

- [LELEC2796 - Wireless Communications](https://uclouvain.be/en-cours-2023-lelec2796). 
- [LELEC2796 - Wireless Communications - Moodle](https://moodle.uclouvain.be/course/view.php?id=1465).

## Authors
- [@Baptiste Sambon](https://github.com/BaptisteSambon) 
- [@Justin Weemaels](https://github.com/Just1Wmls)