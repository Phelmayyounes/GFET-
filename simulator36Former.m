hold off;
clear;
%Define useful constants
   
    electronCharge = 1.602e-19; %C
    planckConstant = 6.62607e-34; %J*s
    fermiVelocity = 10^6; %m/s [30] Thiele
    beta = ((electronCharge)^3)/(pi*(fermiVelocity*planckConstant/(2*pi))^2); %
    
%Parameters

    Lox = 10e-6; %m  
    Wox = 5e-6 ; %m
    Tox = 38.2e-9; %m
    Er = 16; %No unit
    u = 0.657; %[m^2/Vs]
    Nf = 2.3e16; %m^-2
    w = (280e-3)*1.602e-19/(planckConstant/(2*pi)); %frecuency 1/
    spatialHom = (65e-3)*1.6022e-19; %J 

% Other parameters
    Npuddle = ((spatialHom)^2)/(((planckConstant/(2*pi))*fermiVelocity)^2*pi);  %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2

%Begin simulation%

    simSize = 100;
    for Vgate = [-1.25,-0.75,0.25,0.25,0.75]
        Vgs = Vgate;
        %Just an identation to localize the heart of the simulation of each
        %Vgs
          Vds = -linspace(0,1,simSize);
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
       %Now we plot 
        figure (1);
        hold on;
        plot(-Vds, -Id);
        xlabel('-Vdsi [V]');
        ylabel('-Ids [mA/um]');
    end  
    