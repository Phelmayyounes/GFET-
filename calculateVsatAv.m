function [VsatAV] = calculateVsatAv (w, Qav, e, Npuddle)
  temp1 = pi*abs(Qav)*e^-1 + Npuddle;
  VsatAV = w*(temp1).^-0.5;
end
 