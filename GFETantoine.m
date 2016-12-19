%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSTANTS

% beta
e = 1.6e-19 % [C]
pi = 3.1416
hBarre = 1.054e-34 % [J.s]
fermiVelocity = 1e6 % [m.s-1]
beta = (e^3)/(pi*(hBarre*fermiVelocity)^2)

% Ctop
vacuumPermitivity = 8.854e-12 % [F.m-1]
relativePermitivity = 3.4 % no unit
oxideThickness = 8.5e-9 % [m]
Ctop = vacuumPermitivity*relativePermitivity/oxideThickness % [F.m-2]

% npuddle
delta = 66.8e-3*1.6e-19 % [J]
npuddle = (delta^2)/(pi*hBarre^2*fermiVelocity^2) % [C.m-2]

% omega 
omega = 56e-3*1.6e-19/hBarre % [s-1]

% mobility
mobility = 7000e4 % [V.m-2.s-1]

% length of the channel
channelLength = 1e-6 % [m]

% width of the channel
channelWidth = 1e-6 %[m]

% Gate volage 
gateVoltage = 0 % [V]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAIN VOLTAGE ARRAY 
drainVoltageMin = 0.1 % [V]
drainVoltageMax = 1 % [V]
drainVoltageArray = linspace(drainVoltageMin,drainVoltageMax,100)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAIN CURRENT ARRAY COMPUTATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DENOMINATOR COMPUTATION

% absQnetArray
absQnetArray = abs((beta*((-Ctop+(Ctop^2.+4*beta*abs(Ctop*(gateVoltage-drainVoltageArray/2))).^0.5)))/(2*beta))

% denominatorArray
denominatorArray = channelLength + (mobility/omega)*((pi*abs(absQnetArray)/e+npuddle).^0.5).*drainVoltageArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERATOR COMPUTATION

% TnumArray
z2Array = Ctop*(gateVoltage-drainVoltageArray/2) % <0
z1 = Ctop*gateVoltage % < 0
TnumArray = (-1/(12*e*beta^2*Ctop))*((6*beta*Ctop^2*z2Array-6*beta^2*z2Array+Ctop*(Ctop^2-4*beta*z2Array).^1.5)-(6*beta*Ctop^2*z1-6*beta^2*z1+Ctop*(Ctop^2-4*beta*z1).^1.5))

% numeratorArray
numeratorArray = (e*mobility*channelWidth)*(TnumArray+npuddle*drainVoltageArray)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAIN CURRENT COMPUTATION

% drainCurrentArray
drainCurrentArray = numeratorArray./denominatorArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
plot(drainVoltageArray,drainCurrentArray)
