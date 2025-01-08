# free-brownian-motion
Brownian dynamics simulation of nanoparticle diffusion in a solvent - (a C++ code)
# Brownian dynamics simulation of nanoparticle diffusion in water
#### Video Demo:  https://www.youtube.com/watch?v=zcFOznMuqDo
#### Description:
Brownian motion is the random, erratic movement of microscopic particles suspended in a fluid, resulting from collisions with the fluid's molecules. First observed by Robert Brown in 1827, this phenomenon is a key example of stochastic processes in physics. It plays a fundamental role in explaining diffusion, thermal motion, and the behavior of particles at the microscopic scale.
(More information: https://en.wikipedia.org/wiki/Brownian_motion)

Brownian dynamics simulation is a computational technique used to study the motion of particles influenced by Brownian motion and other forces in a fluid. By integrating the Langevin equation, these simulations account for stochastic thermal forces and deterministic interactions, such as hydrodynamic and interparticle forces. Brownian dynamics is widely used in fields like soft matter physics, biophysics, and nanotechnology to model phenomena such as diffusion, self-assembly, and the behavior of colloidal systems. Its efficiency in capturing the essential dynamics without explicitly simulating the solvent makes it a powerful tool for studying complex systems.
(more information: https://en.wikipedia.org/wiki/Brownian_dynamics)

The Stokes-Einstein relation is a fundamental equation used to estimate the diffusivity of nanoparticles (NPs) in a solvent. It links the diffusion coefficient D of a particle to the thermal energy k_B*T, the particle's hydrodynamic radius r, and the solvent's viscosity 𝜂, through the expression D=(k_B*T)/(6*𝜋*𝜂*𝑟). This relation assumes that the particle is spherical, the fluid is homogeneous, and the motion is dominated by thermal fluctuations. Widely applied in nanotechnology, biophysics, and colloidal science, the Stokes-Einstein relation provides a convenient way to predict how environmental factors like temperature and viscosity affect the mobility of nanoparticles in a solvent.
(more information: https://en.wikipedia.org/wiki/Einstein_relation_(kinetic_theory) )

In our simulation, we release N particles within a periodic unit cell with a side length b=200nm. The particles undergo random Brownian motion, and their trajectories are tracked by monitoring the center of mass (COM) over time. Using these trajectories, we calculate the mean square displacement (MSD) as a function of time, which provides insights into the particles' diffusive behavior. From the MSD data, we derive the diffusivity of the particles and compare these results with theoretical predictions. This approach allows us to validate the simulation and explore the influence of various factors on the system's dynamics.
(more on MSD: https://en.wikipedia.org/wiki/Mean_squared_displacement)

When the code is run, all the details of simulation is extracted in txt files which can be used to visualize trajectories, plot MSD vs time. You may add more interactions (electrostaitc) to the system. All parameters of simulation can be found here at the end of simulation: parameters.txt. You may change simulation parameters such as time step, periodic box size, number of particles... Adding more particles may result in better accuracy.
This code can be extended to reproduce this paper:
Hansing, Johann, et al. "Nanoparticle filtering in charged hydrogels: Effects of particle size, charge asymmetry and salt concentration." The European Physical Journal E 39 (2016): 1-13.
https://link.springer.com/article/10.1140/epje/i2016-16053-2
In other words, this code is the starting point of the paper by Hansing et al. 2016

When you run the code, you need to enter: T(temperature), mu= viscosity, d= diameter of NP to initiate the simulation
You find some background information when run the code, as well!
Use this code to delete all text files: rm -f *.txt (when you rerun the code, old txt files are replated with updated ones)
Final results (diffusivity, MSD,...) can be found in two files: Results.txt and Internal_TranslationRotation.txt. We used a more advanced formula for calculation of MSD in the second file.
Any question about this can be sent to : mrokhforouz92@gmail.com

