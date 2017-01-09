hold off;
clear;
%Define useful constants
   
    electronCharge = 1.602e-19; %C
    planckConstant = 6.62607e-34; %J*s
    fermiVelocity = 10^6; %m/s [30] Thiele
    beta = ((electronCharge)^3)/(pi*(fermiVelocity*planckConstant/(2*pi))^2); %
    
%Parameters

    Lox = 0.44e-6; %m  
    Nf = 0; %m^-2
    Wox = 1 ; %m
    Tox = 8.5e-9; %m
    Er = 3.4; %No unit
    u = 0.7; %[m^2/Vs]
    w = (65e-3)*1.602e-19/(planckConstant/(2*pi)); %frecuency 1/
    spatialHom = (66.8e-3)*1.6022e-19; %J 

% Other parameters
    Npuddle = ((spatialHom)^2)/(((planckConstant/(2*pi))*fermiVelocity)^2*pi);  %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2

%Begin simulation%

    simSize = 100;
    for Vgate = [-0.5,-1,-1.5,-2]
        Vgs = Vgate;
        %Just an identation to localize the heart of the simulation of each
        %Vgs
          Vds = -linspace(0,0.8,simSize);
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
       %Now we plot 
        Id = Id*1000/1000000;
        figure (1);
        hold on;
        plot(-Vds, -Id,'-','DisplayName', strcat('Vgs = ',num2str(Vgs)));
        grid on;
        xlabel('-Vdsi [V]');
        ylabel('-Ids [mA/um]');
    end
    
   legend('show')
    