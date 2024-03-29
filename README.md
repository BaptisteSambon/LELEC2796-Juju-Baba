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

A specific scenario is defined with shared assumptions for fair evaluation. It involves downlink communication from one Base Station (BS) to several User Equipments (UEs) with fixed initial power, using Binary Phase-Shift Keying (BPSK) modulation. Perfect Channel State Information (CSI) and perfect synchronization are assumed, with a flat fading Line of Sight (LoS) configuration and non-zero Inter-Cell and Inter-User Interference. Minimum Mean Square Error equalization (MMSE) is applied at the receiver. 

The Inter-Cell Interferences (ICI) are generated using the [cells.ipynb](./cells.ipynb) notebook by considering a grid of cells. The number of cells, the density of cells, the obstacles, etc. can be changed. The ICI are then used in the [NOMA.ipynb](./NOMA.ipynb) and [CDMA.ipynb](./CDMA.ipynb) notebooks to simulate the NOMA and CDMA systems. The results are plotted in the [NOMAvsCDMA.ipynb](./NOMAvsCDMA.ipynb) notebook.

## Available codes

The codes in this repository are divided into 4 Jupyter notebooks and 3 python files. 

The **notebooks** are:

- #### [NOMA.ipynb](./NOMA.ipynb)
Contains code for NOMA and plots the BER, the spectral efficiency and the users scalability for different values of the SNR. 

- #### [CDMA.ipynb](./CDMA.ipynb)
Contains code for CDMA and plots the BER, the spectral efficiency and the users scalability for different values of the SNR.

- #### [NOMAvsCDMA.ipynb](./NOMAvsCDMA.ipynb)
Contains code for plots comparing the performance (BER, spectral efficiency, users scalability) of the NOMA and CDMA systems for different values of the SNR, different number of users. 

- #### [cells.ipynb](./cells.ipynb)
Contains the code used to generate the ICI with different scenarios investigated (density of cells, obstacles, etc.).

The **python files** are:

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

1. Iswarya, N., and L. S. Jayashree. 2021. ‘A Survey on Successive Interference Cancellation Schemes in Non-Orthogonal Multiple Access for Future Radio Access’. *Wireless Personal Communications*, 120(2), 1057–78. 

2. Dai, L., Wang, B., Yuan, Y., Han, S., Chih-Lin, I., & Wang, Z. (2015). Non-orthogonal multiple access for 5G: Solutions, challenges, opportunities, and future research trends. *IEEE Communications Magazine*, 53, 74–81. 

3. A. Benjebbour, Y. Saito, Y. Kishiyama, A. Li, A. Harada, and T. Nakamura, ”Concept and practical considerations of non-orthogonal multiple access (NOMA) for future radio access,” *Proc. IEEE ISPACS*, Nov. 2013. 

4. Y. Saito, A. Benjebbour, Y. Kishiyama, and T. Nakamura, “System-level performance evaluation of downlink nonorthogonal multiple access (NOMA),” *Proc. IEEE PIMRC 2013*, Sept. 2013. 

5. A. Li, A. Benjebbour, and A. Harada, "Performance evaluation of non-orthogonal multiple access combined with opportunistic beamforming," *Proc. IEEE VTC Spring 2014*, May 2014. 

6. A. Benjebbour, A. Li, Y. Kishiyama, H. Jiang, and T. Nakamura, “System-level performance of downlink NOMA combined with SU-MIMO for future LTE enhancements,” *IEEE Globecom*, Dec. 2014. 

7. K. Higuchi and A. Benjebbour, “Non-orthogonal multiple access (NOMA) with successive interference cancellation for future radio access,” *IEICE Trans. on Commun.*, Jan. 2015. 

8. Wang, Hong, Zhaoyang Zhang, and Xiaoming Chen. 2017. ‘Energy-Efficient Power Allocation for Non-Orthogonal Multiple Access with Imperfect Successive Interference Cancellation’. In *2017 9th International Conference on Wireless Communications and Signal Processing (WCSP)*, 1–6. 

9. Yang, Zhaohui, Wei Xu, Cunhua Pan, Yijin Pan, and Ming Chen. 2017. ‘On the Optimality of Power Allocation for NOMA Downlinks With Individual QoS Constraints’. *IEEE Communications Letters*, 21(7), 1649–52. 

10. Wang, Chin-Liang, Jyun-Yu Chen, and Yi-Jhen Chen. 2016. ‘Power Allocation for a Downlink Non-Orthogonal Multiple Access System’. *IEEE Wireless Communications Letters*, 5(5), 532. 

11. López, Carlos Alberto Rodríguez, and Vitalio Alfonso Reguera. 2021. ‘Simple Fair Power Allocation for NOMA-Based Visible Light Communication Systems’. *ArXiv*. 

12. L. Hanzo, L-L. Yang, E-L. Kuan and K. Yen. 2003. 'Single- and Multi-Carriers DS-CDMA: Multi-User Detection, Space-Time Spreading, Synchronisation, Networking and Standards'. Chapter 2 'CDMA Overview'. 

13. Nguyen Thi-Mai-Trang. "Génération des codes CDMA dans l'UMTS". "Network and Performance Analysis Team". Posted in 2017. Accessed in november 2023. 

14. Pei Xiao, Erik Ström. Added to IEEE Xplore in March 2004. 'Delay Estimation and Data Detection in Long-Code DS\_CDMA System'. DOI:10.1109/APCC.2003.1274425. 



## Links

- [LELEC2796 - Wireless Communications](https://uclouvain.be/en-cours-2023-lelec2796). 
- [LELEC2796 - Wireless Communications - Moodle](https://moodle.uclouvain.be/course/view.php?id=1465).

## Authors
- [@Baptiste Sambon](https://github.com/BaptisteSambon) 
- [@Justin Weemaels](https://github.com/Just1Wmls)