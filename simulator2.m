hold off;
clear;
%Define useful constants
   
    electronCharge = 1.602e-19; %C
    planckConstantEVs = 6.62607e-34; %J*s
    fermiVelocity = 10e6; %m/s [30] Thiele
    beta = ((electronCharge)^3)/(pi*(fermiVelocity*planckConstantEVs/(2*pi))^2); %
    Nf = 0; %m^-2 
    Wox = 1; %Just to ignore it since its always a factor
    u = 0.7; %[m^2/Vs] 
    Er = 3.4; %No unit
    Tox = 8.5e-9; %m 
    Lox = 1e-6; %m  
    w = (56e-3)*1.602e-19/(planckConstantEVs/(2*pi)); %frecuency 1/s
    spatialHom = (66.8e-3)*1.6022e-19; %J 
    Npuddle = ((spatialHom)^2)/((((planckConstantEVs/(2*pi))*fermiVelocity))^2*pi);  %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2

%Begin simulation%
    
    %Fig. 3  Frégonèse et al. 2013%
    simSize = 100;
    for Vgate = -[0.2, 0.5, 0.8, 1, 1.2, 1.4]
        Vgs = Vgate;
        %Just an identation to localize the heart of the simulation of each
        %Vgs
          Vds = -linspace(0,0.99,simSize);
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
       %Now we plot 
        Id = Id*1000/100000;
        figure (1);
        hold on;
        plot(-Vds, -Id);
        xlabel('-Vdsi [V]');
        ylabel('-Ids [mA/um]');
    end %Simulation for the figure 3 in paper of Fregonese 
    
    
    %Fig. 4, 5  Frégonèse et al. 2013%
    simSize = 100;
    Vgs = linspace(-2,2,simSize);
    
    for Vdrain_source = [0.1, 0.55, 1]
        Vds = Vdrain_source;
        %Now we ident the simulation of Qch as a function of Vgs
            %Fig 4.
            s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
            Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
            denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
            numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
            Id = (calculateId(numerator, denominator));
            VsatAv = calculateVsatAv(w, Qav, electronCharge, Npuddle);
            Eav = calculateEav(Qav, electronCharge, Npuddle, u, Wox, Id, VsatAv);  
            Qch = calculateQch(electronCharge, Wox, Eav, Npuddle, Ctop, beta, Nf, Vds, Vgs);
            %end Fig. 4.
            
          
            
        %Now we plot the results
        figure(2);
        hold on;
        Qch = (Qch/(Lox*100*100*electronCharge))*10^-12; %Wox is not necessary since it cancells during computations 
        plot(Vgs, Qch, '.','markersize',5); 
        xlabel('Vgsi [V]');
        ylabel('Qch/L*W*q*10e12 [cm^-2]');
        %end Fig. 4
        
        %Fig. 5
        derQch = diff(Qch)*100/4*electronCharge*10^12*10^6;  %the factor 4/100 is due to the step size of Vgs
        derQch(100) = derQch(99); %to make the vector of the correct size to be plot
        
        figure(3);
        hold on;
        %derQch = derQch; %uF/cm^-2
        plot(Vgs, -derQch, '.','markersize',5); 
        xlabel('Vgsi [V]');
        ylabel('dQch/dVgs [uF/cm^-2]');
       %end Fig 5.
        
    end %end Fig. 4, 5  Frégonèse et al. 2013%
   
    
    %Fig. 7  Frégonèse et al. 2013%
    %Redefine useful parameters
    Nf = 2.3e16; %m^-2 
    Wox = 5e-6; %m
    u = 0.6570; %[m^2/Vs] 
    Er = 16; %No unit
    Tox = 38.2e-9; %m 
    Lox = 10e-6; %m  
    w = 280e-3*1.602e-19/(planckConstantEVs/(2*pi)); %frecuency 1/s
    spatialHom = 65e-3*1.602e-19; %J 
    Npuddle = ((spatialHom)^2)/((((planckConstantEVs/(2*pi))*fermiVelocity)^2)*pi); %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2
    
    
    simSize = 100;
    for Vgate = -[-0.5, 0.1, 0.5,0.75,1]
        Vgs = Vgate;
        Vds = linspace(0,1,simSize); %We checked for both signs already and
        %We indent to show that here starts the simulation 
       
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
      
        %Now we start the plot
        figure (4);
        hold on;
        Id = Id*10^6;
        plot(Vds, Id,'.');
        xlabel('Vdsi [V]');
        ylabel('Ids [uA]');
       
    end