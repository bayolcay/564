clear all;
delete all;

%% Permeability vs Field Strength/Flux Density vs Field Strength Curves for Kool M?® MAX (?26)

Hx=0:1:700;
Ux=0;
Bx=0;
% fit formulas for B and ?
for i=1:length(Hx)
    Ux(i)=1/(0.01 + 5.700e-08*Hx(i)^2.205 );
    Bx(i,1)= Hx(i)*100;                                                     % creation of B-H curve for FEMM unit conversion to AT/m
    Bx(i,2)=(  (8.741e-02+1.634e-02*Hx(i)+7.844e-04*Hx(i)*Hx(i))  / (1+1.044e-01*Hx(i)+6.576e-04*Hx(i)*Hx(i)) )^1.814;
end
% plot(Ux)
 
%figure; 
% plot(Bx(:,2),'LineWidth',2)
xlabel('H (AT/cm)');
ylabel('B (T) ')
title('B-H Curve of Kool Mu^® MAX (26u)') 


annotation('line',[.125 .42],[.1 .69])
annotation('line',[.43 .8],[.69 .92])
annotation('line',[.42 .42],[.125 .75],'LineStyle',':','LineWidth',2)

 
TxtLoc1=[0.15 0.35 0.01 0.3]             % x   y   w    h
TxtLoc2=[0.5 0.4 0.01 0.3]               % x   y   w    h

annotation('textbox',TxtLoc1,'String','Linear Region','FitBoxToText','on','FontSize',10);
annotation('textbox',TxtLoc2,'String', 'Saturation Region','FitBoxToText','on','FontSize',10);

 grid on;
 grid minor;
 
 %% Intitialization data for largest product in the Kool M?® MAX family 
 Partno=79337;
 U0=4*pi*1e-7
 Uinit=26;
 
 CoreHeight=25.4*1e-3;                                                      %m
 InnerDia=78.60*1e-3;                                                       %m
 OuterDia=132.60*1e-3;                                                      %m
 WindowArea=4710*1e-6;                                                      % m2
 CrossSection= 678*1e-6;                                                    % m2
 MeanPathLength=324*1e-3;                                                   % m
 
 Temperature=100;                                                           %C According to Datasheet, KoolMu-26u permeability does not change considerably wrt T
 IDC=40;
 
 %% Selection of operating point
 % since powder core materials show "soft saturation"  characteristics
 % it is not an easy process to select a saturation region
 % bacause of that an operating point where B-H curve slope considerably
 % reduced and winding factor stays below 0.5 has been selected

 
 FieldStrength_AtIdc=300*1e2;                                               %A*T/m
 TurnsCount=FieldStrength_AtIdc*MeanPathLength/IDC;                         % too many turns!!! 

%% PartA.1     Inductance Calculation with homogeneous flux distribution and linear BH characteristics

 Ur=Uinit*U0;
 HomoLinearFluxDensity=Ur*TurnsCount*IDC/MeanPathLength;
 Total_Flux_Homo=HomoLinearFluxDensity*CrossSection
 L_Homo_Linear= TurnsCount*Total_Flux_Homo/IDC;
 
%% PartA.2     Inductance Calculation with non-homogeneous flux distribution and linear BH characteristics

FluxDensity_Linear=  @(r) Ur*TurnsCount*IDC./(pi*r)

Total_Flux_NonHomo = CoreHeight*integral(@(r)FluxDensity_Linear(r),InnerDia*0.5,OuterDia*0.5)

L_NonHomo_Linear = TurnsCount*Total_Flux_NonHomo/IDC

 
%% PartA.3     Inductance Calculation with homogeneous flux distribution, non-linear BH characteristics and DC current is increased by 50%
IDC1_5= IDC*1.5;
FieldStrength_At1_5_Idc=IDC1_5*TurnsCount/MeanPathLength;
Unonlinear=0.01*1/(0.01 + 5.700e-08*(FieldStrength_At1_5_Idc*1e-2)^2.205 ); %calculation of permeability at IDC with given permeability fitter
                                                                            %function given in terms of cm                                                                             %output given as %
Ur=Unonlinear*Uinit*U0;
HomononLinearFluxDensity=Ur*FieldStrength_At1_5_Idc;
Total_Flux_Homo_Nonlinear=HomononLinearFluxDensity*CrossSection
L_HomoNonlinear= TurnsCount*Total_Flux_Homo_Nonlinear/IDC1_5;
 
%% PartA.4     Inductance Calculation with non-homogeneous flux distribution, non-linear BH characteristics and DC current is increased by 50%
 
FieldStrength_atR=@(r) IDC1_5*TurnsCount./(pi*r*1e2);

Permeability_atR=@(r)  Uinit*U0*0.01 ./ ( 0.01 + 5.700e-08 * FieldStrength_atR(r).^2.205); 

FluxDensity_NonHomo_NonLinear= @(r) Permeability_atR(r)*TurnsCount*IDC1_5./(pi*r) 
Total_Flux_NonHomo_nonLinear = CoreHeight*integral(@(r)FluxDensity_NonHomo_NonLinear(r),InnerDia*0.5,OuterDia*0.5);
 
 L_NonHomo_nonLinear = TurnsCount*Total_Flux_NonHomo_nonLinear/IDC1_5
 
 
 %% PartA.5     Inductance Calculation of Gapped Toroid with homogeneous flux distribution, linear BH characteristics 
 
 gap= 2e-3;
 Ur=Uinit*U0;
 Total_Flux_nonFringing=CrossSection*TurnsCount*IDC/(gap/U0+(MeanPathLength-gap)/Ur)
 L_NonFringing = TurnsCount*Total_Flux_nonFringing/IDC
 
%% PartA.6     Inductance Calculation of Gapped Toroid with homogeneous flux distribution, linear BH characteristics 
%air gap crossection should be taken larger 
%fringing area would be proportional to the gap length
%therefor gap dimensions can be enlargen by the gap length
%as the total flux leaving in the core would be equal to the flux in the air gap

%F=TotalFlux*totalReuctance

%L_Fringing = TurnsCount*TotalFlux/IDC;
 

 
 
 