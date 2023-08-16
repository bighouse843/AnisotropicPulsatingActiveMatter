# AnisotropicPulsatingActiveMatter

In this repository, you can find the scripts of my master's thesis project on the role of anisotropy in pulsating active matter. As four models have been studied, there are four different folders. In each folder, there is a C++ script for running the simulations and a Julia script to visualize the evolution of the simulated system. Furthermore, there are a .c and a .h file to execute the .cpp file. 
To execute the .cpp file run the command:
`g++ -o out Param_order_ord.cpp`
Then to launch the simulation run for example
`./out test 50 50 1.8 1 1 1 1 10 0.0001 100 0 1 1 1234 0.5 1 1 1 0.1 0.001 0.1 0.1 0.1`
The inputs are: output Lx Ly rho T mu_r mu_phi mu_theta omega dt totaltime initial_time interval_record interval_snap seed radius amplitude_r epsilon amplitude_theta lambda dist drho dx dphi dtheta

