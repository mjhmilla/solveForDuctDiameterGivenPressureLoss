function f = evaluateFrictionFactor(d,k,Re)

f = 0.01;
fp=f;
g = 0.1;
i=0;

while abs(g)>1e-12 && i < 100

    g = -0.5/sqrt(f) - log( (k/(3.7*d)) + (2.51/(Re*sqrt(f))));
    dg_df = 0.25/(f*sqrt(f)) ...
        + 1.255/( log(10)*Re*(f*sqrt(f))*( (k/(3.7*d)) + (2.51/(Re*sqrt(f))))  );

    df = -g/dg_df;
    fp = f;
    f=f+df;

    %(k/(3.7*d)) + (2.51/(Re*sqrt(f))) = 2
    %   fmin = ((2.51/(2-(k/(3.7*d))))/Re)^2
    fmin = ((2.51/(2-(k/(3.7*d))))/Re)^2;

    if(f < fmin)
      f=0.5*(fmin+fp);
    end

    i=i+1;
end

assert(abs(g)<=1e-12,'Error: failed to solve the root for the friction factor');
