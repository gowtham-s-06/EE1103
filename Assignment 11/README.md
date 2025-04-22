Switching of a magnetic particle is commonly used in data storage. The switching dynamics are typically described by the Landau-Lifshitz-Gilbert Equation (LLGE).   
  
{\displaystyle {\frac {\mathrm {d} \mathbf {M} }{\mathrm {d} t}}=-\gamma \mathbf {M} \times \mathbf {H_{\mathrm {eff} }} -\lambda \mathbf {M} \times \left(\mathbf {M} \times \mathbf {H_{\mathrm {eff} }} \right)}  
  
For a spherical particle, the switching can be simplified using polar coordinates. Use equations (8) and (18),  from Mallinson's paper, to solve for the trajectory theta(t) and phi(t) for  M on the unit sphere, and compare your solutions against what you obtain using cartesian coordinates.  
  
Start with small angle theta_start with respect to z, and apply a field along -z  
For the same timestep delta_t, use  Heun's method to estimate theta(t) and phi(t)   
Or, you can use the RK45 code available here after suitably modifying it for your problem  
Assume H >> Hk so you can ignore Hk, and alpha < 0.2.  
Stop your simulation when theta > theta_stop  
Repeat the above step, but in cartesian coordinates i.e. solve for dmx/dt, dmy/dt and dmz/dt. [You must normalize M after every time step to continue to stay on the sphere.]  
Calculate the R^2 difference between the two methods  
Estimate the total time (in nanoseconds) it takes to switch between the two states  
Inputs: theta_start theta_stop alpha delta_t  
Output: alpha R^2 switch_time_ns  
Plot the trajectory (either solution) on the unit sphere and label your axis  
Filenames:  rollnum_llg.c, rollnum_llg.jpg  
  
Optional:   
  
Profile your code (using separate functions allows you to study how the code performs)  
Plot the switching time (between theta_stop = 170deg and theta_start=10deg) versus alpha  
Read and understand the physics of Mallinson's paper. Repeat the simulation for Hk ~ 0.1 H  
