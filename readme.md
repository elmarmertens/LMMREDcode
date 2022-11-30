# Code for "Indeterminacy and Imperfect Information" by Thomas Lubik, Christian Matthes, and Elmar Mertens

Review of Economic Dynamics (in press): [https://doi.org/10.1016/j.red.2022.09.003](https://doi.org/10.1016/j.red.2022.09.003)


All codes were executed using Matlab R2022b on a Windows 10 desktop. The files listed below can be run in any order. We use the default Matlab random number generator. The seed is set at the top of each file using the rng command. For the equilibrium specific impulse responses (figures 4 and 9) the seed matters - the specific plots produced by the corresponding figures should therefore be thought of as illustrative only.

The code requires access to Matlab's control system toolbox. Should the toolbox be unavailable, please set useDare=false in each script. 

## List of files to run (names of files correspond to figures in the main text):

1. figure_1.m in folder Figures_1_and_2 (982 seconds runtime)
2. figure_2.m in folder Figures_1_and_2 (1021 seconds runtime)
3. figures_3_and_5.m in folder Figures_3_4_5_8_9_10 (367 seconds runtime)  
4. figure_4.m in folder Figures_3_4_5_8_9_10 (18 seconds runtime)
5. figure_6.m in folder Figures_6_and_7 (21 seconds runtime)
6. figure_7.m in folder Figures_6_and_7 (22 seconds runtime)
7. figures_8_and_10.m in folder Figures_3_4_5_8_9_10 (681 seconds runtime)
8. figure_9.m in folder Figures_3_4_5_8_9_10 (5 seconds runtime)
