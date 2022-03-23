# The-Coupled-Crystal-Plasticity-and-Stochastic-Cluster-Dynamics-model-with-Qianran-Yu-
The CP/SCD simulation package contains a computer code written in C++. It couples the 0-dimensional SCD method with a general crystal plasticity formulation. We have used the code to study mechanical properties change of W and Fe metals during irradiation.

****Developed by****

Qianran Yu

Research projects on:
1. Irradiation hardening of self-ion irradiating tungsten material
2. Irradiation creep/swelling of DEMO neutron irradiated Fe material

Please send an email to the following address for more information:

Qianran Yu (yuqianran0709@gmail.com)

Jaime Marian (jmarian@ucla.edu)

****Introduction****

Mechanical properties of materials is subjected to significant degradation under irradiation condition. In CP/SCD model, we coupled a dislocation-mediated crystal plasticity (CP) formulation with the one-point stochastic cluster dynamics method (SCD). The SCD method simulates the microstructure evolutions during irradiation, and the CP model deals with the stress-strain correlations. This code has been used to study irradiation hardening under uniaxial loading condition and high dose irradiation creep behaviors. Please see the following journal papers for details:

[1] Qianran Yu, Sabyasachi Chatterjee, Kenneth J. Roche, Giacomo Po and Jaime Marian, “Coupling crystal plasticity and stochastic cluster dynamics models of irradiation damage in tungsten,” Modeling and Simulation in Materials Science and Engineering. 29, 055021 (2021) (link: https://iopscience.iop.org/article/10.1088/1361-651X/ac01ba)

[2] Qianran Yu, Sabyasachi Chatterjee, Giacomo Po, Jaime Marian, “Simulations of irradiation creep in ferritic materials under DEMO first-wall operation con- dition,” (processing).

****How to use****

Before use: 
1. Please install gcc/g++ for the newest version (c++11 or newer)
2. Please either unzip eigen library and put it into src directory or install eigen  
