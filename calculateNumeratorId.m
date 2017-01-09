function [numerator] = calculateNumeratorId (e, u, W, Ctop, b, Npuddle, Nf, Vgs, Vds)
  
  %Step 1 
  z2 = Ctop*(Vgs - Vds) + e*Nf;
  z1 = Ctop*Vgs + e*Nf;
  
  %Step 2
  %We take the signs of the boundaries z1, z2, to filter according to the
  %sign and apply the formula in the appendix of Fregonese paper 
  signZ1 = sign(z1);
  signZ2 = sign(z2);
  
  %Step 3
  %We make to vectors, one with the negative values of z1, and cero
  %elsewhere, and other with the positive values of z1, and cero elsewhere.
  negativeZ1 = z1.*((1-signZ1)/2); %The second coefficient drops to cero for positive values of z1, and is 1 if z<0
  positiveZ1 = z1.*((1+signZ1)/2); %The second coefficient drops to cero for negative values of z1, and is 1 if z>0
  
  %Step 4
  %Negative values of Z1
  temporal11 = 6*(Ctop^2)*b*negativeZ1 - 6*(b^2)*(negativeZ1.^2) + Ctop*((Ctop^2) - 4*b*negativeZ1).^(3/2);
  temporal11 = temporal11.*((1-signZ1)/2); %We make cero the calculations for non-negative values of z1
  
  %Step 5
  %Positive values of Z1
  temporal12 = 6*(Ctop^2)*b*positiveZ1 + 6*(b^2)*(positiveZ1.^2) - Ctop*((Ctop^2) + 4*b*positiveZ1).^(3/2) + 2*(Ctop^4);
  temporal12 = temporal12.*((1+signZ1)/2); %We make cero the calculations for non-positive values of z1
  
  %Step 6
  %We sume the contributions of negative and positive values of z1
  temporal1 = temporal11 + temporal12;
  
  %Step 3 again, for z2
  negativeZ2 = z2.*((1-signZ2)/2); %The second coefficient drops to cero for positive values of z1, and is 1 if z<0
  positiveZ2 = z2.*((1+signZ2)/2); %The second coefficient drops to cero for negative values of z1, and is 1 if z>0
  
  %Step 4 for z2
  %Negative values of z2
  temporal21 = 6*(Ctop^2)*b*negativeZ2 - 6*(b^2)*(negativeZ2.^2) + Ctop*((Ctop^2) - 4*b*negativeZ2).^(3/2);
  temporal21 = temporal21.*((1-signZ2)/2); %We make cero the calculations for non-negative values of z2
  
  %Step 5 for z2
  %Positive values of z2
  temporal22 = 6*(Ctop^2)*b*positiveZ2 + 6*(b^2)*(positiveZ2.^2) - Ctop*((Ctop^2) + 4*b*positiveZ2).^(3/2) + 2*(Ctop^4);
  temporal22 = temporal22.*((1+signZ2)/2); %We make cero the calculations for non-positive values of z2
  
  %we end the if condition in the vectorized form
  
    temporal2 = temporal21 + temporal22;
  
    Fnum = (-1/(12*e*(b^2)*Ctop))*(temporal2-temporal1);
  
    numerator = e*u*W*(Fnum + Npuddle*Vds);
  
  end