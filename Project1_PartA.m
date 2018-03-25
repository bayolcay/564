clear all;
delete all;

%% Permeability vs Field Strength/Flux Density vs Field Strength Curves for Kool M?® MAX (?26)

Hx=0:1:700;
Ux=0;
Bx=0;
% fit formulas for B and ?
for i=1:length(Hx)
    Ux(i)=1/(0.01 + 5.700e-08*Hx(i)^2.205 );
    Bx(i)=(  (8.741e-02+1.634e-02*Hx(i)+7.844e-04*Hx(i)*Hx(i))  / (1+1.044e-01*Hx(i)+6.576e-04*Hx(i)*Hx(i)) )^1.814;
end
% plot(Ux)
 
%figure; 
% plot(Bx,'LineWidth',2)
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
 
 Temperature=100;                                                           %C According to Datasheet curves KoolMu-26u permeability does not change considerably wrt T
 IDC=40;
 
 %% Selection of operating point
 % since powder core materials show "soft saturation"  characteristics
 % it is not an easy process to select a saturation region
 % bacause of that an operating point where B-H curve slope considerably
 % reduced and winding factor stays below 0.5 has been selected
 
 FieldStrength_AtIdc=300*1e2;                                               %A*T/m
 TurnsCount=FieldStrength_AtIdc*MeanPathLength/IDC;                          %mm to cm conversion of MeanPathLength

%% PartA.1     Inductance Calculation with homogeneous flux distribution and linear BH characteristics

 Ur=Uinit*U0;
 HomoFluxDensity=Ur*TurnsCount*IDC/MeanPathLength;
 Total_Flux_Homo=HomoFluxDensity*CrossSection
 L_Homo_Linear= TurnsCount*Total_Flux_Homo/IDC;
 
%% PartA.2     Inductance Calculation with non-homogeneous flux distribution and linear BH characteristics

FluxDensity=  @(r) Ur*TurnsCount*IDC./(pi*r)

Total_Flux_NonHomo = CoreHeight*integral(@(r)FluxDensity(r),InnerDia*0.5,OuterDia*0.5)

L_NonHomo_Linear = TurnsCount*Total_Flux_NonHomo/IDC

 
%% PartA.3     Inductance Calculation with non-linear core and the DC current is increased by 50%

 Ur=Uinit*U0;
 HomoFluxDensity=Ur*TurnsCount*IDC/MeanPathLength;
 Total_Flux_Homo=HomoFluxDensity*CrossSection
 L_Homo_Linear= TurnsCount*Total_Flux_Homo/IDC;
 
 
 
 
 
 
 
 
 
 
 
 