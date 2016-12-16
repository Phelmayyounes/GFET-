function [Qav] = calculateQav (beta, Ctop, Vgs, Vds, e, Nf,s)
  temp1 = Ctop*(Vgs-Vds/2) + e*Nf;
  temp2 = Ctop^2 + 4*beta*abs(temp1);
  temp3 = (-Ctop + sqrt(temp2));
  Qav = beta*s.*(temp3/(2*beta)).^2;
 
  
end
