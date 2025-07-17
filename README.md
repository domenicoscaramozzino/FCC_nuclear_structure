# FCC nuclear structure
This Matlab code contains an adaptation of the Face-Centered-Cubic (FCC) lattice model proposed by Norman D. Cook ("Models of the Atomic Nucleus: Unification Through a Lattice of Nucleons", 2nd Edition, 2010) for visualizing and modeling the structures of atomic nuclei. The code contains an interactive interface to insert input data.

First, you will be asked which nucleus you want to analyze, specifying the number of protons (Z) and neutrons (N). For example, in case of Fe56 you type Z = 26 and N = 30. Following the default build-up proposed by Cook, the code will place protons and neutrons in the proper FCC positions and calculate basic features such as the Binding Energy (BE).

![FCC](images/Fe56_FCC.png)

Then the code will carry out a randomization of positions of the nucleons lying on the external surface. You will be asked, for both protons and neutrons, which are the minimum and maximum n-values that define the borders of the shell layer where nucleon positions can change randomly. Then the code will generate thousands of additional FCC structures, that have the same spin as the default build-up and do not exhibit disconnected nucleons. For example, inserting n = 1 and n = 4 as minimum and maximum values for both protons and neutrons for Fe56, one obtains 43,431 feasible structures, with BE values ranging from ~417 MeV to ~477 MeV.

![FCC](images/BEs_Fe56.png)

The code will also provide a picture of the randomized FCC structure with the highest BE (in this case, ~477 MeV).

![FCC](images/Optimal_BE_Fe56_FCC.png)

Then you can carry out the fission analysis of the selected nucleus. The code will carry out fission calculations by fracturing each of the feasible FCC structures along the main 21 crystallographic planes of the FCC lattice and generate a list of fragments, associated to Z, N, A numbers, Coulomb repulsion and energy values needed to break the first-neighboring bonds. Note: for very small nuclei (i.e. He4) you might need to reduce the number of crystallographic planes, since some of them might lie outside the overall surface of the nucleus and will not generate separate fragments. The code will then plot a histogram distribution of fission energies, e.g. for Fe56 these are found to range from ~15 MeV to ~147 MeV.

![FCC](images/Fission_energies_Fe56.png)

Then, you can select a energy range (in MeV) to carry out a statistical analysis 
