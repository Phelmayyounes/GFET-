hold off;
clear;
%Define useful constants
   
    electronCharge = 1.602e-19; %C
    planckConstant = 6.62607e-34; %J*s
    fermiVelocity = 10^6; %m/s [30] Thiele
    beta = ((electronCharge)^3)/(pi*(fermiVelocity*planckConstant/(2*pi))^2); %
    
%Parameters

    Lox = 2.94e-6; %m  
    Wox = 4e-6 ; %m
    Tox = 28e-9; %m
    Er = 3.9; %No unit
    u = 0.68; %[m^2/Vs]
    Nf = 2.5e8; %m^-2
    w = (370e-3)*1.602e-19/(planckConstant/(2*pi)); %frecuency 1/
    spatialHom = (66.8e-3)*1.6022e-19; %J 

% Other parameters
    Npuddle = ((spatialHom)^2)/(((planckConstant/(2*pi))*fermiVelocity)^2*pi);  %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2

%Begin simulation%

    simSize = 100;
    for Vgate = [3,4,5,7]
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
        Id = Id*1000/1000000;
        figure (1);
        hold on;
        plot(-Vds, -Id);
        xlabel('-Vdsi [V]');
        ylabel('-Ids [mA/um]');
    end  
    