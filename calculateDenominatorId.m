function [denominator] = calculateDenominatorId (Vds, u, Qav, e, Npuddle, L, w)

    denominator = L + u*abs(Vds).*((pi*(abs(Qav)/e+Npuddle)).^0.5)/w;
    
end
 