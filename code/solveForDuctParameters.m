function ductParams = solveForDuctParameters(...
                          ductName,ductNumber, ...
                          mdot,L,k,rho,nu,deltaPTarget,...
                          dmin,dmax,...
                          maxIterations, numericalTolerance,...
                          fileName, appendToFile,...
                          flag_generatePlot,...
                          plotFolder,...
                          usingOctave)


soln = evaluatePressureLoss(dmin,mdot,rho,L,nu,k);

assert(soln.deltaP > deltaPTarget, ...
'Pick a smaller starting diameter');

soln = evaluatePressureLoss(dmax,mdot,rho,L,nu,k);

assert(soln.deltaP < deltaPTarget, ...
'Pick a bigger ending diameter');

%Now use the bisection method to solve our root

dBest           = 0.5*(dmin+dmax);
soln            = evaluatePressureLoss(dBest,mdot,rho,L,nu,k);
deltaPErrorBest = abs(deltaPTarget-soln.deltaP);

dDelta= (dmax-dmin)*0.25;


for i=1:1:maxIterations
  dL          = dBest - dDelta;
  solnL       = evaluatePressureLoss(dL,mdot,rho,L,nu,k);
  deltaPError = abs(deltaPTarget - solnL.deltaP);

  if(deltaPError < deltaPErrorBest)
    soln            = solnL;
    dBest           = dL;
    deltaPErrorBest = deltaPError;
  else
    dR          = dBest + dDelta;
    solnR       = evaluatePressureLoss(dR,mdot,rho,L,nu,k);
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
disp(sprintf('\t%1.3f m\thydraulic diameter',dBest));


if(appendToFile==1)
  fid =fopen(fileName,'a');
else
  fid =fopen(fileName,'w');
end

fprintf(fid,'%s\n',ductName);
fprintf(fid,'\t%1.6f\tkg/s\t%s\n', mdot,          'mdot');
fprintf(fid,'\t%1.6f\tm\t\t%s\n',     L,            'L');
fprintf(fid,'\t%1.6f\tm\t\t%s\n',     k,            'k');
fprintf(fid,'\t%1.6f\tm/s\t\t%s\n',   soln.v,       'v');
fprintf(fid,'\t%1.6f\tm^2\t\t%s\n',   soln.A,       'A');
fprintf(fid,'\t%1.0f\t\t%s\n',      soln.Re,      'Re');
fprintf(fid,'\t%i\t\t\t%s\n',      (soln.Re>4000),      'Re>4000');
fprintf(fid,'\t%1.6f\t%s\n',      soln.f,       'f');
fprintf(fid,'\t%1.6f\tPa\t%s\n',    deltaPTarget, 'delta_P_Target');
fprintf(fid,'\t%1.2e\tPa\t%s\n',    deltaPErrorBest, 'delta_P_Error');
fprintf(fid,'\t%1.6f\tm\t%s\n\n',   dBest,        'd');
fclose(fid);

ductParams = soln;
ductParams.deltaPTarget = deltaPTarget;
ductParams.deltaPError  = deltaPError;
ductParams.d            = dBest;

if(flag_generatePlot==1)
  n=101;
  dSeries = dmin + ([0:(1/(n-1)):1]' .* (dmax-dmin) ); %Initial guess for the hydralic diameter

  solnSeries = struct('A',zeros(n,1),'v',zeros(n,1),...
                    'Re',zeros(n,1),'f',zeros(n,1),'deltaP',zeros(n,1));

  for i=1:1:length(dSeries)
      solnTemp = evaluatePressureLoss(dSeries(i,1),mdot,rho,L,nu,k);

      solnSeries.A(i,1)   = solnTemp.A;
      solnSeries.v(i,1)   = solnTemp.v;
      solnSeries.Re(i,1)  = solnTemp.Re;
      solnSeries.f(i,1)   = solnTemp.f;
      solnSeries.deltaP(i,1) = solnTemp.deltaP;


      if(i==1)
          assert(solnSeries.deltaP(i,1) > deltaPTarget, ...
          'Pick a smaller starting diameter');
      end
      if(i==length(dSeries))
          assert(solnSeries.deltaP(i,1) < deltaPTarget, ...
          'Pick a bigger ending diameter');
      end

  end



  fig=figure;
    %yyaxis left;
    subplot(2,2,1)
      plot(dSeries,solnSeries.deltaP);
      hold on;
      plot(dBest,deltaPTarget,'o');
      hold on;
      xlabel('Hydraulic Diameter (m)');
      ylabel('Pressure Loss (Pa)');
      title('Numerically Solved Duct Pressure Loss vs Hydraulic Diameter');
      ylim([0,deltaPTarget*1.5]);
      box off;

    subplot(2,2,2)
      plot(dSeries,solnSeries.v);
      hold on;
      plot(dBest,soln.v,'o');
      hold on;
      xlabel('Hydraulic Diameter (m)');
      ylabel('Velocity (m/s)');
      title('Fluid velocity');
      ylim([0,soln.v*1.5]);
      box off;

    subplot(2,2,3)
      plot(dSeries,solnSeries.Re);
      hold on;
      plot(dBest,soln.Re,'o');
      hold on;
      title('Reynolds number of fluid');
      xlabel('Hydraulic Diameter (m)');
      ylabel('Reynolds Number (unitless)');
      ylim([0,soln.Re*1.5]);
      box off;

    %yyaxis right;
    subplot(2,2,4)
      plot(dSeries,solnSeries.f);
      hold on;
      plot(dBest,soln.f,'o');
      hold on;
      title('Numerically Solved Duct Friction Factor vs Hydraulic Diameter');
      xlabel('Hydraulic Diameter (m)');
      ylabel('Friction Factor (unitless)');
      ylim([0,soln.f*1.5]);
      box off;

    if(usingOctave==1)
      print (fig, fullfile(plotFolder,...
                           sprintf('fig_%i_%s.pdf',ductNumber,ductName)),...
             "-dpdflatexstandalone");
    else
        fileName = sprintf('fig_%i_%s.pdf',ductNumber,ductName);
        print('-dpdf', fullfile(plotFolder,fileName));
    end
    close(fig);
end



here=1;

