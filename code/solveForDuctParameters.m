function ductParams = solveForDuctParameters(...
                          ductName,ductNumber, ...
                          mdot_guess,L,k,rho,nu,deltaPTarget,...
                          dmin,dmax,...
                          useFrictionApproximation,...
                          maxIterations, numericalTolerance,...
                          fileName, appendToFile)


soln = evaluatePressureLoss(dmin,mdot_guess,rho,L,nu,k,useFrictionApproximation);

assert(soln.deltaP > deltaPTarget, ...
'Pick a smaller starting diameter');

soln = evaluatePressureLoss(dmax,mdot_guess,rho,L,nu,k,useFrictionApproximation);

assert(soln.deltaP < deltaPTarget, ...
'Pick a bigger ending diameter');

%Now use the bisection method to solve our root

dBest           = 0.5*(dmin+dmax);
soln            = evaluatePressureLoss(dBest,mdot_guess,rho,L,nu,k,useFrictionApproximation);
deltaPErrorBest = abs(deltaPTarget-soln.deltaP);

dDelta= (dmax-dmin)*0.25;


for i=1:1:maxIterations
  dL          = dBest - dDelta;
  solnL       = evaluatePressureLoss(dL,mdot_guess,rho,L,nu,k,useFrictionApproximation);
  deltaPError = abs(deltaPTarget - solnL.deltaP);

  if(deltaPError < deltaPErrorBest)
    soln            = solnL;
    dBest           = dL;
    deltaPErrorBest = deltaPError;
  else
    dR          = dBest + dDelta;
    solnR       = evaluatePressureLoss(dR,mdot_guess,rho,L,nu,k,useFrictionApproximation);
    deltaPError = abs(deltaPTarget - solnR.deltaP);
    if(deltaPError < deltaPErrorBest)
      soln            = solnR;
      dBest           = dR;
      deltaPErrorBest = deltaPError;
    end
  end

  dDelta = dDelta*0.5;

  if(deltaPErrorBest < numericalTolerance)
    break;
  end

end



%%
% Do not touch the code below-------------------------
%%

disp(ductName);
disp(sprintf('\t%1.2f\tPa\tTarget Pressure',deltaPTarget));
disp(sprintf('\t%1.3e\tPa\tSolution Error',deltaPErrorBest));

disp(sprintf('\t%1.6f\tm^3\tarea',      soln.A));
disp(sprintf('\t%1.6f\tm/2\tvelocity',  soln.v));
disp(sprintf('\t%1.6f\t\tRe',           soln.Re));
disp(sprintf('\t%1.6f\t\tfriction factor',soln.f));
disp(sprintf('\t%1.3f\t\tm\thydraulic diameter',dBest));


if(appendToFile==1)
  fid =fopen(fileName,'a');
else
  fid =fopen(fileName,'w');
end

fprintf(fid,'%s\n',ductName);
fprintf(fid,'\t%1.6f\tkg/s\t%s\n', soln.mdot,       'mdot');
fprintf(fid,'\t%1.6f\tkg/s\t%s\n', soln.mdot_error, 'mdot_error');
fprintf(fid,'\t%1.6f\tm\t%s\n',     L,            'L');
fprintf(fid,'\t%1.6f\tm\t%s\n',     k,            'k');
fprintf(fid,'\t%1.6f\tm/s\t%s\n',   soln.v,       'v');
fprintf(fid,'\t%1.6f\tm^2\t%s\n',   soln.A,       'A');
fprintf(fid,'\t%1.0f\t\t\t%s\n',      soln.Re,        'Re');
fprintf(fid,'\t%i\t\t\t%s\n',      (soln.Re>4000),  'Re>4000');
fprintf(fid,'\t%1.6f\t\t%s\n',      soln.f,           'f');
fprintf(fid,'\t%1.6f\tPa\t%s\n',    deltaPTarget,   'delta_P_Target');
fprintf(fid,'\t%1.2e\tPa\t%s\n',    deltaPErrorBest,'delta_P_Error');
fprintf(fid,'\t%1.6f\tm\t%s\n\n',   dBest,          'd');
fclose(fid);

ductParams = soln;
ductParams.deltaPTarget = deltaPTarget;
ductParams.deltaPError  = deltaPError;
ductParams.d            = dBest;



here=1;

