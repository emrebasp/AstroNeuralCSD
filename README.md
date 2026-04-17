# AstroNeuralCSD
An astro-neural field model with application to cortical spreading depression 


CSD propagation simulations based on an astro-neural field model presented in the article

[1] "An astro-neural field model with application to cortical spreading depression"

by E. Baspinar, D. Avitabile, C. Nouveau, M. Desroches, F. Campillo, M. Mantegazza

This set of Matlab codes provide the numerical framework for the simulations of the
model presented in the aforementioned article. We refer to the article for details 
regarding the model and the experimental cases which are simulated via the Matlab
codes. The codes were prepared in Matlab R2023a (64-bit).

The code was designed and edited by Daniele Avitabile, VU Amsterdam, Netherlands (d.avitabile@vu.nl) and by Emre Baspinar, INRIA, MathNeuro Team, France (emre.baspinar@inria.fr).

We thank to Fabien Campillo, INRIA, MathNeuro Team, for his advices in preparation of this set of Matlab codes.

USAGE:

This folder contains a set of m-files for performing simulations regarding propagation of cortical spreading depression (CSD), as well as the subfolders which are used by these m-files.

Run main.m to run the simulations providing the data and the figures corresponding to the results presented in Figures 3-7 of the article. Use blockRun.m to produce space-time plots where parameters gamma_an, gamma_aa and gamma_a are swept separately.

> main.m: Runs the m.files performing the simulations which give the results presented in [1]. An index between 1-7 should be chosen in main.m to choose which m.file to run.

> ISO_or_GBZ.m: It performs the simulations regarding the experimental conditions (C1) and (C2) (see [1]). These experiments are related to the effects of isoguvacine (ISO) and gabazine (GBZ) on CSD propagation speed. This m.file can be independently run to reproduce Figure 3a in the article instead of using main.m with the corresponding index.

> ChR2.m: It performs the simulations regarding the experimental condition (C3) (see the article). This experiment explores the effects of channelrhodopsin (ChR2) activation on CSD propagation speed. This m.file can be independently run to reproduce Figure 3b in the article instead of using main.m with the corresponding index.

> GBZ_plus_ChR2.m: It simulates the experimental condition (C4). It is based on simultaneous application of GBZ and activation of ChR2. This m.file can be independently run to reproduce Figure 3c in the article instead of using main.m with the corresponding index.

> sweep_gamma_aa.m: It performs the simulations providing the space-time plots of the model state variables where gamma_aa is varied in a given range. This m.file can be independently run to reproduce Figure 4 in the article instead of using main.m with the corresponding index.

> sweep_gamma_a.m: It performs the simulations providing the space-time plots of the model state variables where gamma_a is varied in a given range. This m.file can be independently run to reproduce Figure 5 in the article instead of using main.m with the corresponding index.

> beyond_propagation.m: It performs the simulations providing the space-time plots of the model state variables corresponding to the simulation scenario where astrocyte increased activity extends beyond the CSD propagation. This m.file can be independently run to reproduce Figure 6 in the article instead of using main.m with the corresponding index.

> astrocyte_stops.m: It performs the simulations providing the space-time plots of the model state variables corresponding to the simulation scenario where astrocyte increased activation stops CSD propagation. This m.file can be independently run to reproduce Figure 7 in the article instead of using main.m with the corresponding index.

> blockRun.m: It performs simulations where gamma_an, gamma_aa and gamma_an are varied separately. It provides the space-time plots of the state variables and the corresponding transfer functions.

> loadPlot.m: It can be employed for generating the plots by using the saved data in "Data" folder.

> control_Parameters_ICs_0.m: It can be generate and save the control parameters, as well as the default initial conditions, in the “Parameter” folder as .mat files. It need not be run as long as the control values are not changed.


Default simulation parameters correspond to the results presented in the article. For a fast inspection of the model, the number of discretization nodes can be reduced (typically not below 500 nodes for Lx = 50).


FOLDERS:

> Data: It contains the data of the simulations performed by the m-files given above, except for csd.m and csd_Control.m. Here csd.m is to play with the parameters, to manipulate the model to see how the changes made by the user affect the results. The m-file csd_Control.m corresponds to only one parameter set, it is very quick to run the simulation. Naming of the csv files in Data provides the dates and the names of the simulations performed:
	- GBZ_or_ISO_DATE_TIME.csv                     <— csd_GBZ_or_ISO.m
	- ChR2_DATE_TIME.csv                           <— csd_ChR2.m
	- GBZ_Plus_ChR2_DATE_TIME.csv                  <— csd_GBZ_Plus_ChR2.m
You can directly import and use them in Matlab to manipulate the figures given 	in Figure 5.4 and Figure 5.6 in the master thesis.

> Parameters: It contains two subfolders, “modelParameters” and “simulationParameters”. In “modelParameters”, you find the .mat file which contains the model parameters corresponding to the control case. In “simulationsParameters”, you find the .mat file containing the simulation parameters which were used to run the m-files producing Figures 5.4 and 5.6 in the master thesis.

> Plots: It contains the saved figures.

> PlotFunctions: It contains the functions which are used for plotting purposes.

> SimulationFunctions: It contains the functions which are crucial for running the m-files for the model simulations.


## License
See LICENSE files in the repository.
