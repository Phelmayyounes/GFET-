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
    Lox = 3e-6; %m  
    w = (56e-3)*1.602e-19/(planckConstantEVs/(2*pi)); %frecuency 1/s
    spatialHom = (66.8e-3)*1.6022e-19; %J 
    Npuddle = ((spatialHom)^2)/((((planckConstantEVs/(2*pi))*fermiVelocity))^2*pi);  %1/m^2
    Ctop = calculateCtop(Er, Tox); %F/m^2

     
%Flow calculation example

    %s = sign(Vgs, Vds, Ctop, Nf, electronCharge);
    %Qav = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s);
    %denominator = calculateDenominatorId(Vds, u, Qav, electronCharge, Npuddle, Lox, Wox); ##this is just the denominator in ec. (5)
    %numerator = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds ); #This is the intregral in the numerator of ec 5 
    %Id = calculateId(numerator, denominator);
    %VsatAv = calculateVsatAv(w, Qav, electronCharge, Npuddle);
    %Eav = calculateEav(Qav, e, Npuddle, u, Wox, Id, VsatAv);  
    %Qch = calculateQch(e, Wox, Eav, Npuddle, Ctop, beta, Nf, Vds, Vgs);

%%now we perform the simulation for some values of Vgate 
%%We obtain a plot of Ids vs. Vds. regarding Fig. 3 Frégonèse et al.
%%2013
    
    
    simSize = 100;

    for Vgate = 0:0.2:1.6
        Vgs = Vgate;
        s = zeros(simSize,1);
        Qav = zeros(simSize, 1);
        denominator = zeros(simSize, 1);
        numerator = zeros(simSize, 1);
        Id = zeros(simSize, 1);
        Fnum = zeros(simSize, 1);
        Vds = linspace(1/100,1,100);
        
       
          s(i) = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav(i) = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s(i));
          denominator(i) = calculateDenominatorId(Vds, u, Qav(i), electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          [numerator(i), Fnum(i)] = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id(i) = (calculateId(numerator(i), denominator(i)));
        

        Id = Id*1000/100000;
        figure (1);
        hold on;
        plot(Vds, Id);
       
    end

 %%now we calculate some values of Qch as functionn of Vds, Vgs as in the
 %%figure 4
    simSize = 100;
    s = zeros(simSize,1);
    Qav = zeros(simSize, 1);
    denominator = zeros(simSize, 1);
    numerator = zeros(simSize, 1);
    Id = zeros(simSize, 1);
    Fnum = zeros(simSize, 1);
    VsatAv = zeros(simSize, 1);
    Eav = zeros(simSize, 1);
    Qch = zeros(simSize, 1);
   
    
    
    for Vdrain_source = [0.1, 0.55, 1]
        Vds = Vdrain_source;
        Vgs = linspace(-2,2,simSize); 
        for i = 1:simSize
            
            s(i) = sign(Ctop*(Vgs(i)-Vds/2) + electronCharge*Nf);
            Qav(i) = calculateQav(beta, Ctop, Vgs(i), Vds, electronCharge, Nf, s(i));
            denominator(i) = calculateDenominatorId(Vds, u, Qav(i), electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
            [numerator(i), Fnum(i)] = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs(i), Vds); %This is the intregral in the numerator of ec 5 
            Id(i) = (calculateId(numerator(i), denominator(i)));
            VsatAv(i) = calculateVsatAv(w, Qav(i), electronCharge, Npuddle);
            Eav(i) = calculateEav(Qav(i), electronCharge, Npuddle, u, Wox, Id(i), VsatAv(i));  
            Qch(i) = calculateQch(electronCharge, Wox, Eav(i), Npuddle, Ctop, beta, Nf, Vds, Vgs(i));

        end

        figure(2);
        hold on;
        Qch = (Qch/(Lox*100*100*electronCharge))*10^-12; %Wox is not necessary since it cancells during computations 
        plot(Vgs, Qch, '.','markersize',10); %scaling of Qch??
    
    end
    
    
%%now we perform the simulation for some values of Vgate 
%%We obtain a plot of Vgs vs. Vds. regarding Fig. 7

%Re-Define useful parameters
    Nf = 2.3e14; %m^-2 
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
    
    for Vgate = [-0.75,-0.5,-0.25,0,0.25,0.5,0.75]
        Vgs = Vgate;
        s = zeros(simSize,1);
        Qav = zeros(simSize, 1);
        denominator = zeros(simSize, 1);
        numerator = zeros(simSize, 1);
        Id = zeros(simSize, 1);
        Fnum = zeros(simSize, 1);

        for i = 1:simSize
          Vds = -i*1/100;
          s(i) = sign(Ctop*(Vgs-Vds/2) + electronCharge*Nf);
          Qav(i) = calculateQav(beta, Ctop, Vgs, Vds, electronCharge, Nf, s(i));
          denominator(i) = calculateDenominatorId(Vds, u, Qav(i), electronCharge, Npuddle, Lox, w); %this is just the denominator in ec. (5)
          [numerator(i), Fnum(i)] = calculateNumeratorId(electronCharge, u, Wox, Ctop, beta, Npuddle, Nf, Vgs, Vds); %This is the intregral in the numerator of ec 5 
          Id(i) = (calculateId(numerator(i), denominator(i)));
        end
        
        Id = Id*10^6;
        Vds = linspace(1/100,1,100);
        figure (3);
        hold on;
        plot(-Vds, -Id,'.');
       
    end