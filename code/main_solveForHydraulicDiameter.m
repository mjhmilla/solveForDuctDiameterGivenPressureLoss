clc;
clear all;
close all;

%%
% Check that we are in the correct directory
%%
rootDir = pwd();
i0 = strfind(rootDir, filesep);
parentDir = rootDir((i0(1,(end-1))+1):end);
assert(strcmp(parentDir,['solveForDuctDiameterGivenPressureLoss',filesep,'code']),...
      'Error: start in the solveForDuctDiameterGivenPressureLoss/code folder');
rootDir  = rootDir(1,1:(i0(1,(end))-1));
%%
%Inputs---------------------------------------------
%%

%Beware: Octave seems to crash windows when plotting.
usingOctave       = 0;
flag_generatePlot = 1;

maxIterations =100;
numericalTolerance = 1e-12;

%%
%Set up the output folder and files
%%
outputFolder  = fullfile(rootDir,'output');
fileName      = fullfile(outputFolder, 'ductSolutionLog.txt');


%%
% Enter duct data
%%

ductStruct(1) = struct('rho',0,'nu',0,'mdot',0,'L',0,'k',0,'deltaPTarget',0,...
                       'dmin',0,'dmax',0,'ductName','');



rho     = 1.2; %kg/m^3
nu      = 0.00001524; %m^2/s


idx=1;

volumePerHour = 60; %m³/h

ductStruct(idx).ductName = 'RWk2Du';
ductStruct(idx).rho   = rho;
ductStruct(idx).nu    = nu;
ductStruct(idx).mdot  = volumePerHour*rho/3600; %kg/s
ductStruct(idx).L     = 27.2;
ductStruct(idx).k     = 0.00007;
ductStruct(idx).deltaPTarget = 3.3;
ductStruct(idx).dmin         = 0.1;
ductStruct(idx).dmax         = 0.3;





%%
% Process all ducts
%%

appendToFile = 0;
for i=1:1:length(ductStruct)
  ductParams = solveForDuctParameters(...
                ductStruct(i).ductName, ...
                i,...
                ductStruct(i).mdot,...
                ductStruct(i).L,...
                ductStruct(i).k,...
                ductStruct(i).rho,...
                ductStruct(i).nu,...
                ductStruct(i).deltaPTarget,...
                ductStruct(i).dmin,...
                ductStruct(i).dmax,...
                maxIterations,...
                numericalTolerance,...
                fileName, ...
                appendToFile);
  appendToFile=1;

  disp(sprintf('%1.3e\t%s',ductParams.d,'d'));

  disp('Checking alternative mdot formula');

  disp(sprintf('%1.3e\t%s',ductParams.mdot,'mdot'));
  disp(sprintf('%1.3e\t%s',ductParams.mdot_guess,'mdot_guess'));
  disp(sprintf('%1.3e\t%s',ductParams.mdot_error,'mdot error'));


   if(flag_generatePlot==1)
     n=101;
     dmin = ductStruct(i).dmin;
     dmax = ductStruct(i).dmax;
     mdot = ductStruct(i).mdot;
     deltaPTarget= ductStruct(i).deltaPTarget;
     L = ductStruct(i).L;
     k = ductStruct(i).k;
     dSeries = dmin + ([0:(1/(n-1)):1]' .* (dmax-dmin) ); %Initial guess for the hydralic diameter
   
     solnSeries = struct('A',zeros(n,1),'v',zeros(n,1),...
                   'Re',zeros(n,1),'f',zeros(n,1),'deltaP',zeros(n,1));
   
     for j=1:1:length(dSeries)
        %disp(j);
        solnTemp = evaluatePressureLoss(dSeries(j,1),mdot,rho,L,nu,k);
   
        solnSeries.A(j,1)   = solnTemp.A;
        solnSeries.v(j,1)   = solnTemp.v;
        solnSeries.Re(j,1)  = solnTemp.Re;
        solnSeries.f(j,1)   = solnTemp.f;
        solnSeries.deltaP(j,1) = solnTemp.deltaP;
   
   
        if(j==1)
           assert(solnSeries.deltaP(j,1) > deltaPTarget, ...
           'Pick a smaller starting diameter');
        end
        if(j==length(dSeries))
           assert(solnSeries.deltaP(j,1) < deltaPTarget, ...
           'Pick a bigger ending diameter');
        end
   
     end
   
   
   
     fig=figure;
       %yyaxis left;
       subplot(2,2,1)
        plot(dSeries,solnSeries.deltaP);
        hold on;
        plot(ductParams.d,deltaPTarget,'o');
        hold on;
        xlabel('Hydraulic Diameter (m)');
        ylabel('Pressure Loss (Pa)');
        title('Numerically Solved Duct Pressure Loss vs Hydraulic Diameter');
        ylim([0,deltaPTarget*1.5]);
        box off;
   
       subplot(2,2,2)
        plot(dSeries,solnSeries.v);
        hold on;
        plot(ductParams.d,ductParams.v,'o');
        hold on;
        xlabel('Hydraulic Diameter (m)');
        ylabel('Velocity (m/s)');
        title('Fluid velocity');
        ylim([0,ductParams.v*1.5]);
        box off;
   
       subplot(2,2,3)
        plot(dSeries,solnSeries.Re);
        hold on;
        plot(ductParams.d,ductParams.Re,'o');
        hold on;
        title('Reynolds number of fluid');
        xlabel('Hydraulic Diameter (m)');
        ylabel('Reynolds Number (unitless)');
        ylim([0,ductParams.Re*1.5]);
        box off;
   
       %yyaxis right;
       subplot(2,2,4)
        plot(dSeries,solnSeries.f);
        hold on;
        plot(ductParams.d,ductParams.f,'o');
        hold on;
        title('Numerically Solved Duct Friction Factor vs Hydraulic Diameter');
        xlabel('Hydraulic Diameter (m)');
        ylabel('Friction Factor (unitless)');
        ylim([0,ductParams.f*1.5]);
        box off;
   
       if(usingOctave==1)
        print (fig, fullfile(plotFolder,...
                         sprintf('fig_%i_%s.pdf',i,ductStruct(i).ductName)),...
              "-dpdflatexstandalone");
       else
          fileName = sprintf('fig_%i_%s.pdf',i,ductStruct(i).ductName);
          print('-dpdf', fullfile(outputFolder,fileName));
       end
       close(fig);
   end
  
end
