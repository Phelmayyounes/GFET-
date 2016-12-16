function [Qch] = calculateQch(e, W, Eav, Npuddle, Ctop, b, Nf, Vds, Vgs)
  

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
  negativeZ1 = z1.*floor((1-signZ1)/2); %The second coefficient drops to cero for positive values of z1, and is 1 if z<0
  positiveZ1 = z1.*ceil((1+signZ1)/2); %The second coefficient drops to cero for negative values of z1, and is 1 if z>0
  
  %Step 4
  %Negative values of Z1
  temporal11 = 6*b*(Ctop^2)*negativeZ1 - 6*(b^2)*(negativeZ1.^2) + Ctop*((Ctop^2) - 4*b*negativeZ1).^(3/2);
  temporal11 = temporal11.*floor((1-signZ1)/2); %We make cero the calculations for non-negative values of z1
  
  %Step 5
  %Positive values of Z1
  temporal12 = -6*b*(Ctop^2)*positiveZ1 - 6*(b^2)*(positiveZ1.^2) + Ctop*((Ctop^2) + 4*b*positiveZ1).^(3/2);
  temporal12 = temporal12.*ceil((1+signZ1)/2); %We make cero the calculations for non-positive values of z1
  
  %Step 6
  %We sume the contributions of negative and positive values of z1
  temporal1 = temporal11 + temporal12;
  
  %Step 3 again, for z2
  negativeZ2 = z2.*floor((1-signZ2)/2); %The second coefficient drops to cero for positive values of z1, and is 1 if z<0
  positiveZ2 = z2.*ceil((1+signZ2)/2); %The second coefficient drops to cero for negative values of z1, and is 1 if z>0
  
  %Step 4 for z2
  %Negative values of z2
  temporal21 = 6*b*(Ctop^2)*negativeZ2 - 6*(b^2)*(negativeZ2.^2) + Ctop*((Ctop^2) - 4*b*negativeZ2).^(3/2);
  temporal21 = temporal21.*floor((1-signZ2)/2); %We make cero the calculations for non-negative values of z2
  
  %Step 5 for z2
  %Positive values of z2
  temporal22 = -6*b*(Ctop^2)*positiveZ2 - 6*(b^2)*(positiveZ2.^2) + Ctop*((Ctop^2) + 4*b*positiveZ2).^(3/2);
  temporal22 = temporal22.*ceil((1+signZ2)/2); %We make cero the calculations for non-positive values of z2
    
  temporal2 = temporal21 + temporal22;
  %we end the if condition in the vectorized form
  

    
  Xnum = (-1/(12*e*(b^2)*Ctop))*(temporal2-temporal1);
  
  Qch = e*W*((Eav).^-1).*(Xnum + Npuddle*Vds);
 
  %%deprecated
  %if (z1 < 0)
  %  temporal1 =  6*b*(Ctop^2)*z1 - 6*(b^2)*(z1.^2) + Ctop*((Ctop^2) - 4*b*z1).^(3/2);
  %else
  %  temporal1 = -6*b*(Ctop^2)*z1 - 6*(b^2)*(z1.^2) + Ctop*((Ctop^2) + 4*b*z1).^(3/2);
  %end
      
  %if (z2 < 0)
  %  temporal2 = 6*b*(Ctop^2)*z2 - 6*(b^2)*(z2.^2) + Ctop*((Ctop^2) - 4*b*z2).^(3/2);
  %else
  %  temporal2 = -6*b*(Ctop^2)*z2 - 6*(b^2)*(z2.^2) + Ctop*((Ctop^2) + 4*b*z2).^(3/2);
  %end
  
end


 