hold off;
clear;
%Define useful constants
   
    electronCharge = 1.602e-19; %C
    planckConstant = 6.62607e-34; %J*s
    fermiVelocity = 10^6; %m/s [30] Thiele
    beta = ((electronCharge)^3)/(pi*(fermiVelocity*planckConstant/(2*pi))^2); %
    Nf = 0; %m^-2 
    Wox = 1; %Just to ignore it since its always a factor
    u = 0.7; %[m^2/Vs] 
    Er = 3.4; %No unit
    Tox = 8.5e-9; %m 
    Lox = 0.5e-6; %m  
    w = (56e-3)*1.602e-19/(planckConstant/(2*pi)); %frecuency 1/s
    
    spatialHom = (65e-3)*1.6022e-19; %J 
    Npuddle = ((spatialHom)^2)/(((planckConstant/(2*pi))*fermiVelocity)^2*pi);  %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2

%Begin simulation%
    
    %Fig. 3  Frégonèse et al. 2013%
    simSize = 100;
    for Vgate = -[0.5, 0.8, 1, 1.2, 1.4, 2]
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
        Id = Id*1000/1000000;
        figure (1);
        hold on;
        plot(-Vds, -Id,'-','DisplayName', strcat('Vgs = ',num2str(Vgs)));
        xlabel('-Vdsi [V]');
        ylabel('-Ids [mA/um]');
    end %Simulation for the figure 3 in paper of Fregonese 
    legend('show')
    
    
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
        plot(Vgs, Qch ,'-','DisplayName', strcat('Vds = ',num2str(Vds)));
        xlabel('Vgsi [V]');
        ylabel('Qch/L*W*q*10e12 [cm^-2]');
        %end Fig. 4
        
        %Fig. 5
        derQch = diff(Qch)*100/4*electronCharge*10^12*10^6;  %the factor 4/100 is due to the step size of Vgs
        derQch(100) = derQch(99); %to make the vector of the correct size to be plot
        
        figure(3);
        hold on;
        %derQch = derQch; %uF/cm^-2
        plot(Vgs, -derQch ,'-','DisplayName', strcat('Vds = ',num2str(Vds)));
        xlabel('Vgsi [V]');
        ylabel('dQch/dVgs [uF/cm^-2]');
        
         %end Fig 5.
        
   end %end Fig. 4, 5  Frégonèse et al. 2013%
   legend('show')
    
   %Fig. 7  Frégonèse et al. 2013%
    %Redefine useful parameters, Table I
    Nf = 2.3e16; %m^-2 
    Wox = 5e-6; %m
    u = 0.657; %[m^2/Vs] 
    Er = 16; %No unit
    Tox = 38.2e-9; %m
    Lox = 10e-6; %m
    w = (280e-3)*1.602e-19/(planckConstant/(2*pi)); %frecuency 1/s
    spatialHom = (65e-3)*(1.602e-19); %J 
    Npuddle = ((spatialHom)^2)/((((planckConstant/(2*pi))*fermiVelocity)^2)*pi); %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2
    
    
    simSize = 100;
    for Vgate = -[1.25, 0.75, 0.25, -0.25, -0.75]
        Vgs = Vgate; 
        Vds = linspace(0,1.5,simSize); %We checked for both signs already and
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
        plot(Vds, Id,'-','DisplayName', strcat('Vgs = ',num2str(Vgs)));  
        xlabel('Vdsi [V]');
        ylabel('Ids [uA]');
       
    end
    legend('show')
    %end Fig. 7  Frégonèse et al. 2013%
   
    %Fig. 8  Frégonèse et al. 2013%
    %Redefine useful parameters Table II
    Nf = 1.26e16; %m^-2 
    Wox = 4e-6; %m
    u = 0.6800; %[m^2/Vs] 
    Er = 3.9; %No unit
    Tox = 28e-9; %m 
    Lox = 2.94e-6; %m  
    w = 370e-3*1.602e-19/(planckConstant/(2*pi)); %frecuency 1/s
    spatialHom = 65e-3*1.602e-19; %J 
    Npuddle = ((spatialHom)^2)/((((planckConstant/(2*pi))*fermiVelocity)^2)*pi); %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2
    
    
    simSize = 100;
    for Vgate = 0:0.5:3
        Vgs = -Vgate;
        Vds = linspace(0,3,simSize); %We checked for both signs already and
        %We indent to show that here starts the simulation 
       
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
      
        %Now we start the plot
        figure (5);
        hold on;
        Id = Id*10^6/(Wox*10^6); %We can also divide by Wox instead of Lox, to get the correct units and correct values
        plot(Vds, Id,'-','DisplayName', strcat('Vgs = ',num2str(Vgs)));
        xlabel('Vdsi [V]');
        ylabel('Ids [uA/um]');
       
    end
    legend('show')
    %end Fig. 8  Frégonèse et al. 2013%
    
   
    %Fig. 9, 9.1  Frégonèse et al. 2013%
    %Redefine useful parameters Table II
    Nf = 2.5e15; %m^-2 

    
    simSize = 100;
    Vgs = linspace(0,3,simSize);
    
    for Vdrain_source = [1.05, 1.55, 2.05, 2.55, 3.05]
        Vds = Vdrain_source;
        %Now we ident the simulation of Ids as a function of Vgs
          
            s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
            Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
            denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
            numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
            Id = (calculateId(numerator, denominator));
              
        %Now we plot the results
        figure(6);
        hold on;
        Ids = Id*10^6/(Wox*10^6); %We could divide per Wox if that fits the values, since the author did not specify it 
        plot(Vgs, Ids);
        xlabel('Vgs [V]');
        ylabel('Ids [uA/um]');
        
        %gm = Ids./Vgs; %Transconductance, we are not sure if this is the formula
        gm = diff(Ids)*100/4;  %the factor 4/100 is due to the step size of Vgs
        gm(100) = gm(99); %to make the vector of the correct size to be plot
        
        
        figure(7);
        hold on;
        plot(Vgs, gm,'-','DisplayName', strcat('Vds = ',num2str(Vds)));
        xlabel('Vgs [V]');
        ylabel('gm [mS/mm]');
        
    end %end Fig. 9, 9.1  Frégonèse et al. 2013%
    legend('show')
    %Fig. 10.1, 10.2, 10.3  Frégonèse et al. 2013%
    %Redefine useful parameters Table II
    Nf = 0; %m^-2 
    Wox = 1; %The value cancels out and we report it per unit of length
    u = 0.7000; %[m^2/Vs] 
    Er = 3.5; %No unit
    Tox = 8.5e-9; %m 
    Lox = 440e-9; %m  
    w = 56e-3*1.602e-19/(planckConstant/(2*pi)); %frecuency 1/s
    spatialHom = 66.8e-3*1.602e-19; %J 
    Npuddle = ((spatialHom)^2)/((((planckConstant/(2*pi))*fermiVelocity)^2)*pi); %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2
    
    %Fig. 10.1 
    simSize = 100;
    for Vgate = 0:0.25:2
        Vgs = Vgate;
        Vds = linspace(0,1.5,simSize); %We checked for both signs already and
        %We indent to show that here starts the simulation 
       
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
      
        %Now we start the plot
        figure (8);
        hold on;
        Id = Id*10^3/(10^6); 
        plot(Vds, Id,'-','DisplayName', strcat('Vgs = ',num2str(Vgs)));
        xlabel('Vdsi [V]');
        ylabel('Ids [mA/um]');
        
    end
    legend('show')
    %end Fig. 10.1%
    
    %Fig. 10.2
    Lox = 1e-6; %m
    simSize = 100;
    for Vgate = 0:0.25:2
        Vgs = Vgate;
        Vds = linspace(0,1.5,simSize); %We checked for both signs already and
        %We indent to show that here starts the simulation 
       
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
      
        %Now we start the plot
        figure (9);
        hold on;
        Id = Id*10^3/(10^6); 
        plot(Vds, Id,'-','DisplayName', strcat('Vgs = ',num2str(Vgs)));
        xlabel('Vdsi [V]');
        ylabel('Ids [mA/um]');
        
    end
    legend('show')
    %end Fig. 10.2%
    
    %Fig. 10.3
    Lox = 3e-6; %m
    simSize = 100;
    for Vgate = 0:0.25:2
        Vgs = Vgate;
        Vds = linspace(0,1.5,simSize); %We checked for both signs already and
        %We indent to show that here starts the simulation 
       
          s = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
          denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id = (calculateId(numerator, denominator));
      
        %Now we start the plot
        figure (10);
        hold on;
        Id = Id*10^3/(10^6); 
        plot(Vds, Id,'-','DisplayName', strcat('Vgs = ',num2str(Vgs)));
        xlabel('Vdsi [V]');
        ylabel('Ids [mA/um]');
        
    end
    legend('show')
    %}
%end Fig. 10.3%
