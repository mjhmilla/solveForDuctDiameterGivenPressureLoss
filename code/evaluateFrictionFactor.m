function f = evaluateFrictionFactor(d,k,Re, useApproximation)

f = 0.0055*(1+ ((20000*k/d) + (1e6/Re))^(1/3));

if(useApproximation==0)
    
    f = 0.1;
    fp=f;
    g = 0.1;
    i=0;
    
    while abs(g)>1e-12 && i < 100
    
        %Eqn 79 in Duct_equations.pdf from COMIS Funadmentals
        g = -0.5/sqrt(f) - log( (k/(3.75*d)) + (2.51/(Re*sqrt(f))));
        dg_df = 0.25/(f*sqrt(f)) ...
            + 1.255/( log(10)*Re*(f*sqrt(f))*( (k/(3.75*d)) + (2.51/(Re*sqrt(f))))  );
    
        df = -g/dg_df;
        fp = f;
        f=f+df;
        
        fmin  = 0;%(-(((3.7*d)/k)*2.51)/Re)^2;
    
        if(f < fmin)
          f=0.5*(fmin+fp);
        end
    %     if(f > fmax)
    %         f = fmax;
    %     end
    
        i=i+1;
    
    end
    
    assert(abs(g)<=1e-12,'Error: failed to solve the root for the friction factor');

end
