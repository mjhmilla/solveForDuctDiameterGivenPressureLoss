function soln = evaluatePressureLoss(d,mdot,rho,L,nu,k,useFrictionApproximation)

% 1. Given mdot guess, evaluate f
% 2. Iterate over mdot equation until the change in mdot is small
% 3. Evaluate delta p

% mdot =  soln.A * sqrt( (ductStruct(i).deltaPTarget*2*rho)...
%       /(ductParams.f/ductParams.d) );
mdotP = 1;
mdotDelta = 1;
i=1;
for i=1:1:2
    soln.A = pi*(d*0.5)*(d*0.5);
    soln.mdot_guess=mdot;
    soln.v = (mdot/rho)/soln.A;
    soln.Re = soln.v*d / nu;
    soln.f = evaluateFrictionFactor(d,k,soln.Re,useFrictionApproximation);
    soln.deltaP = (soln.f*L/d)*rho*(soln.v*soln.v)/2;
    
    mdotP = mdot;
    soln.mdot = soln.A * sqrt(soln.deltaP * 2 * rho)/(soln.f *L/ d);
    soln.mdot_error = abs(mdot-mdotP);
end