function [Ctop] = calculateCtop(permitivity, thickness)
  e0 = 8.854187817620389e-12; %F/m
 
   %Ctop = e0*permitivity*area/thickness;
   Ctop = e0*permitivity/thickness; %F/m^2
   
  
end