function [Eav] = calculateEav(Qav, e, Npuddle, u, W, Id, Vsat)

  temp1 = ((abs(Qav) + e*Npuddle)*u*W).*(Id.^-1);
  temp2 = u*(Vsat.^-1);
  temp3 = temp1 + temp2;
  Eav = temp3.^-1;  %Equation 7 of the paper
  
end