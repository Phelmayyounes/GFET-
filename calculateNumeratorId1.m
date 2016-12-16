function [numerator] = calculateNumeratorId1 (e, u, W, Ctop, b, Npuddle, Nf, Vgs, Vds)
  z2 = Ctop*(Vgs - Vds) + e*Nf;
  z1 = Ctop*Vgs + e*Nf;
  temporal1 = 6*(Ctop^2)*b*z1 - 6*(b^2)*(z1.^2);
  temporal2 = 6*(Ctop^2)*b*z2 - 6*(b^2)*(z2.^2);
  
  if (z1 < 0)
     temporal1 = temporal1 + Ctop*((Ctop^2) - 4*b*z1).^(3/2);
  else
     temporal1 = temporal1 - Ctop*((Ctop^2) + 4*b*z1).^(3/2) + 2*(Ctop^4);
  end
      
  if (z2 < 0)
    temporal2 = temporal2 + Ctop*((Ctop^2) - 4*b*z2).^(3/2);
  else
    temporal2 = temporal2 - Ctop*((Ctop^2) + 4*b*z2).^(3/2) + 2*(Ctop^4);
  end
  
  Fnum = (-1/(12*e*(b^2)*Ctop))*(temporal2-temporal1);
  
  numerator = e*u*W*(Fnum + Npuddle*Vds);
  
end