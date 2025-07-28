function figH = plotDuctSolution(ductIndex,ductStruct, ductParams,...
                useFrictionApproximation,outputFolder,usingOctave)



     n=101;

     rho = ductStruct.rho;
     nu = ductStruct.nu;
     dmin = ductStruct.dmin;
     dmax = ductStruct.dmax;
     mdot = ductStruct.mdot;
     deltaPTarget= ductStruct.deltaPTarget;
     L = ductStruct.L;
     k = ductStruct.k;
     dSeries = dmin + ([0:(1/(n-1)):1]' .* (dmax-dmin) ); %Initial guess for the hydralic diameter
   
     solnSeries = struct('A',zeros(n,1),'v',zeros(n,1),...
                   'Re',zeros(n,1),'f',zeros(n,1),'deltaP',zeros(n,1));
   
     for j=1:1:length(dSeries)
        %disp(j);
        solnTemp = evaluatePressureLoss(dSeries(j,1),mdot,rho,L,nu,k,...
                                             useFrictionApproximation);
   
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
   
   
   
     figH=figure;
       %yyaxis left;
       subplot(2,2,1)
        plot(dSeries,solnSeries.deltaP);
        hold on;
        plot(ductParams.d,deltaPTarget,'o');
        hold on;
        xlabel('Hydraulic Diameter (m)');
        ylabel('Pressure Loss (Pa)');
        title('Pressure Loss vs Hydraulic Diameter');
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
        title('Friction Factor vs Hydraulic Diameter');
        xlabel('Hydraulic Diameter (m)');
        ylabel('Friction Factor (unitless)');
        ylim([0,ductParams.f*1.5]);
        box off;
   
       if(usingOctave==1)
        print (fig, fullfile(plotFolder,...
                         sprintf('fig_%i_%s.pdf',ductIndex,ductStruct.ductName)),...
              "-dpdflatexstandalone");
       else
          fileName = sprintf('fig_%i_%s.pdf',ductIndex,ductStruct.ductName);
          print('-dpdf', fullfile(outputFolder,fileName));
       end
