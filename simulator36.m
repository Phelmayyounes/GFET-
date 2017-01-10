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
    Nf = 0; %m^-2
    w = (280e-3)*1.602e-19/(planckConstant/(2*pi)); %frecuency 1/
    spatialHom = (65e-3)*1.6022e-19; %J 
    Vgate = [-1.25,-0.75,-0.25,0.25,0.75];
    VdsMin = 0;
    VdsMax = 1.5;
    

% Other parameters
    Npuddle = ((spatialHom)^2)/(((planckConstant/(2*pi))*fermiVelocity)^2*pi);  %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2

%Begin simulation%

    simSize = 100;
    for Vgs = Vgate;
        %Just an identation to localize the heart of the simulation of each
        %Vgs
          Vds = linspace(VdsMin,VdsMax,simSize);
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
       %Now we plot 
        figure (1);
        hold on;
        plot(Vds,Id,'-','DisplayName', strcat('Vgs = ',num2str(Vgs)));
        grid on;
        xlabel('-Vdsi [V]');
        ylabel('-Ids [A]');
    end
    
   legend('show')
    