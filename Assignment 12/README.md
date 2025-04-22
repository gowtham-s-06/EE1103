Simulate the onset of ferromagnetism and antiferromagnetism using a 2D Ising model and the Metropolis-Hastings algorithm.  This is a simple and fun exercise as you begin to see structure in what starts off as a random lattice of spins  
  
Try to use pointers and structures in your code  
Test your code with small lattices first  
Take as input the lattice row length (N), the constant coupling between nearest neighbours (J) and a normalized temperature (T)  
Start with a random distribution of spins, +1 and -1 in a lattice of size NxN with periodic boundary conditions  
Decide to retain, or flip, the value of a spin by calculating the energy  over nearest neighbours [-J*s_i*sum(s_j)]  
If the energy change delta_E <= 0 then flip  
If the energy change delta_E > 0, then flip with a probability P = exp(-delta_E/T)  
Keep 2 copies of the lattice (current and new), and update the new lattice with the new spin configuration  
After you have taken a decision on all the spins in the current lattice, copy the new lattice over to the old one and start a new iteration  
End your simulation when the change in total energy (totalE) over all spins in the lattice does not change significantly (define a convergence criterion totdelE) and report the number of iterations used (maxiter) and totdelE  
Plot total_M versus T and identify the critical temperature Tc (use a text label). Avoid screenshots. Learn to export the image from gnuplot  
Outputs: rollnum N J T maxiter totdelE total_M  
Files: rollnum_ising.c and rollnum_ising.jpg (or .png)  
Optional: Create a movie to show how the lattice is updated with iterations and eventually reaches a ferromagnetic or antiferromagnetic configuration. Share the movie as a link in the Comments while uploading your submission.  
