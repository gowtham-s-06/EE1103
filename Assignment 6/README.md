The Mandelbrot set is the set of complex numbers c for which the function fc(z) = z2 + c does not diverge when iterated from z = 0, i.e, for which the sequence fc(0), fc(fc(0)) etc, remains bounded in absolute value.  
  
The Mandelbrot set is the set of values of c in the complex plane for which the orbit of 0 under iteration of the quadratic map Zn+1 = Zn2 + c remains bounded. That is, a complex number c is part of the Mandelbrot set if, when starting with Z0 = 0 and applying the iteration repeatedly, the absolute value of Zn remains bounded however large n gets. Below given is the initial image of a Mandelbrot set zoom sequence. Black points correspond to numbers that are outside the set.  
  
Use structures to define complex numbers and their operations  
Generate a Mandelbrot set zoom sequence, with black points corresponding to numbers that are outside the set.  
Usage: ./mandelbrot <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.jpg>  
Upload: rollnum_mandelbrot.c  and rollnum_mandel.jpg  
