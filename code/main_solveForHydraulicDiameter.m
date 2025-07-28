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

useFrictionApproximation=1;
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
ductStruct(idx).useFrictionApproximation=0;

idx=idx+1;
ductStruct(idx).ductName = 'RWk2Du_fapprox';
ductStruct(idx).rho   = rho;
ductStruct(idx).nu    = nu;
ductStruct(idx).mdot  = volumePerHour*rho/3600; %kg/s
ductStruct(idx).L     = 27.2;
ductStruct(idx).k     = 0.00007;
ductStruct(idx).deltaPTarget = 3.3;
ductStruct(idx).dmin         = 0.1;
ductStruct(idx).dmax         = 0.3;
ductStruct(idx).useFrictionApproximation=1;



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
                ductStruct(i).useFrictionApproximation,...
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
       figH = plotDuctSolution(i,ductStruct(i), ductParams,...
                ductStruct(i).useFrictionApproximation,outputFolder,usingOctave);
       close(figH);
   end
  
end
