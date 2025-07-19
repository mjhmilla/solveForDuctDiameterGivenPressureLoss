function soln = evaluatePressureLoss(d,mdot,rho,L,nu,k)

soln.A = pi*(d*0.5)*(d*0.5);
soln.v = (mdot/rho)/soln.A;
soln.Re = soln.v*d / nu;

soln.f = evaluateFrictionFactor(d,k,soln.Re);
soln.deltaP = (soln.f*L/d)*rho*(soln.v*soln.v)/2;

