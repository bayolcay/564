%% EE564 - Design of Electrical Machines
%% Project-3: Induction Generator Design
%% Name: Olcay BAY
%% ID: 1673672
%% INTRODUCTION
% In this project, design of a squirrel cage asynchronous
% generator for wind turbines will be performed

%% CONTENTS
%
% # Project specifications and selected main design inputs
% # Calculation of main dimensions
% # Selection of main dimensions and validation of machine loading
% # Selection of stator slot number and turn numbers
% # Validating the results and iterations
% # Calculation of MMF, flux density, winding factors and the resultant
% induced voltage
% # Selection of rotor slot number
% # Selection of stator and rotor conductors
% # Calculation of stator slot and rotor bar dimensions
% # Calculation of equivalent core length and effective air gap distance
% # Calculation of winding and bar resistances, leakage inductances and
% magnetizing inductance
% # Calculation of copper and aluminium losses
% # Calculation of copper and aluminium masses
% # Calculation of iron mass
% # Calculation of core losses and the core loss resistance of equivalent
% circuit
% # Calculation of stator to rotor turns ratio
% # Calculation of efficiency
% # Torque-speed characteristics
% # Determination of basic parameters like starting torque, maximum torque
% etc.
% # Motoranalysis
% # Conclusions
% # References


%%  Specifications of asynchronous squirrel cage induction motor
%%
% Rated Power Output: 250 kW
% Line to line voltage: 400 V
% Power factor: 0.87
% Rated Wind Speed: 14 m/s
% Rated Turbine Speed: 24.3 rpm
% Gear Ratio: 31.2
% Frequency: 50 Hz
% Rated Speed: 758 rpm
% Gearbox: (Coupled from wind turbine blade)
% Insulation Class: F (Max. avg:145, Hotspot: 155 )
% Duty: Continuous running duty (S1)
% Efficiency: IE2,  95-96%
% Ingress protection: IP54
% Connection of motor windings:Y
%% Design Procedure
% Design procedure is given by the following flowchart[1]: 
clear all;
clc;
close all;
I = imread('Design Procedure.png');
figure;
imshow(I);
title('Shear Stress For Different Machines','FontSize',16,'FontWeight','Bold');

%% Main Design Inputs
%% 
Pout_n = 250e3;% Nominal Output Power, watts
pole = 8;%Pole Number
p = pole/2;%Pole Pair 
m = 3;%Number of phases
Vn = 400; % Nominal line voltage (l-l)
Vphase = Vn/sqrt(3); % Nominal phase voltage
fn = 50; %Nominal Frequency, Hz
Nsync = 120*fn/pole; % rpm
Nrated = 758; %Nominal rotor speed, rpm
gear_ratio = 31.2; %
pf_expected = 0.86; %
n_expected = 0.95; %
%
u0 = 4*pi*1e-7;
rho_cu = 1.7e-8; % Ohm*m @25C
Kfe=0.92; % Lamination Factor, using 3mm laminations
%
%% Machine Size
%%
Sout_expected=Pout_n/pf_expected;%Expected apparent output power
Pmec = Pout_n/n_expected; % watts 
Irated = Sout_expected/(sqrt(3)*Vn); % amps
omega_rot = Nrated*2*pi/60; % rated angular speed of rotor, rad/sec
Trated = Pmec/omega_rot; % Rated torque,
Power_pp=Pout_n/p; % Power per pole-pair, watts
% Power_pp:62500, p:4  from graph given below, Cmech between 200-250.
I = imread('Cmec.png');
figure;
imshow(I);
title('Specific Machine Constant (Cmec) vs power/pole-pair','FontSize',16,'FontWeight','Bold');
%taken as
Cmech = 220; % kWs/m^3
% Calculation of the aspect ratio 
X = (pi/pole)*(p)^(1/3); % aspect ratio
fprintf('Aspect ratio is %g\n',X);
% diameter^2*length can be calculated by using the Cmec relation:
D2L = Pout_n*1e-3/(Cmech*(Nsync/60)); %
% From these two information, the inner diameter and length can be calculated:
Di_s = (D2L/X)^(1/3); % Stator inner diameter, m
Di_s=round(Di_s,3); % round up to 3 decimal places 
L = Di_s*X; % Rotor-Stator Length, m
L=round(L,3); % round up to 3 decimal places 
% For a 8 pole machine, outer diameter of the stator can be found using:  
% [T.Miller - Electric Machine Design Course, Lecture-5, Slide4]
% The air gap calculation a scale factor of 1.6 is added for heavy duty operation.
g = 1e-3*1.6*(0.18+0.006*Pout_n^0.4); % m
g=  round(g, 4); % round up to 4 decimal places
fprintf('Inner diameter of the stator is %g m\n',Di_s);
fprintf('Length of the machine is %g m\n',L);
fprintf('\nAir gap distance is %g mm\n\n',g);


%% Stress Factors
%%
% Calculation of the tangential force using rated torque:
Ftan = Trated/(Di_s*0.5); % Newton
surface_area = pi*Di_s*L; % m^2
% calculation of tangantial stress:
StressTangent = 1e-3*Ftan/surface_area; % kPa
fprintf('Tangential force is %g Newtons\n',Ftan);
fprintf('Tangential stress is %g kPascals\n',StressTangent);
% Resulting factor (Ftan_s: 23 kPa) is around given average value in figure:
I = imread('shear stress.png');
figure;
imshow(I);
title('Shear Stress For Different Machines','FontSize',16,'FontWeight','Bold');
%---------------------------------------------------------------------------------------------------------------------------------------
% The magnetic loading of the machine is selected from the table 6.2 of textbook

I = imread('Permitted Flux Densities.png');
figure;
imshow(I);
title('Reasonable Flux Densities','FontSize',16,'FontWeight','Bold');
Bgap = 0.75; % peak magnetic loading of air gap , T 
fprintf('\nAir gap magnetic loading is: %g Tesla\n',Bgap);
% Then the electric loading becomes 
electric_loading = StressTangent/Bgap; % kA/m
fprintf('Resultant electric loading is %g kA/m\n',electric_loading);
% Comparing with the reference values given in Table 6.3 in texbook
I = imread('Electrical Loading.png');
figure;
imshow(I);
title('current densities(J) and linear current densities(A) for various electrical machines','FontSize',16,'FontWeight','Bold');
% As seen, resulting value is acceptable  

%% Selection of stator slot number 
% Stator winding is decided to be a double layer winding, 
%stator slot pitch for asynchronous machines is typically in the range of 7-45 mm
% Then minimum and maximum number of slots can be calculated by:
StatorInnerCircumference = pi*Di_s; % m
Qs_max = floor(StatorInnerCircumference/0.007);
Qs_min = ceil(StatorInnerCircumference/0.045);
fprintf('Possible stator slot numbers: %g -%g \n',Qs_min,Qs_max);
% Then possible stator slot numbers resulting an integer slot machine 
%can be found by
k=1;
Qs=0;
qs=1;
fprintf('Possible stator slot numbers(Qs) and resulting qs values:\n');
while (Qs<Qs_max)
    Qs=m*pole*qs;
    if (Qs>Qs_min&&Qs<Qs_max)
        fprintf(' Qs=%d, qs = %d\n',Qs,qs);
    end
    qs=qs+1;
end 

% In order to construct a smooth MMF waveform, a large Qs is desirable,
% Considering manufacturing difficulties, 96 or 120 would be optimal choice
% As a starting point Qs taken as 96
Qs=96;
qs=Qs/m/pole;
%for pole>4 recomended rotor slot number,Qr given by Qr=1.2Qs (from lecture notes) 
Qr=round(1.2*Qs);
Qr=Qr-mod(Qr,2); % Odd values may cause vibration

fprintf(' Selected slot numbera Qs=%d, qs=%d, Qr=%d\n',Qs,qs,Qr);
%However Qr must be checked in order to prevent unwanted Qs-Qr pairs causing slot harmonics:
if ( Qr==(Qs+2*p) ||  Qr==(Qs-2*p) || Qr==(2*Qs+2*p) || Qr==(2*Qs-2*p) ||  Qr==(Qs+p) || Qr==(Qs-p) ||   Qr==(0.5*Qs+p) || Qr==(0.5*Qs-p)    )
        fprintf(' -----------------------------------improper Qs-Qr pair!!! \n');
end
fprintf('---------------------------------------------------------\n');


 %Then for two layer winding,
PitchRatio=1;
fprintf('\n\nWinding Factors for a coil span of 12/12:  \n');

for HarmonicOrder=1:2:13
    SlotAngle=(2*pi/(Qs/p));%  radians
    kd=sin(qs*HarmonicOrder*(SlotAngle*0.5)) / (  qs*sin(HarmonicOrder*(SlotAngle*0.5)) );
    kp=sin(HarmonicOrder*pi/2*PitchRatio);
    kw=kp*kd;
    if mod(HarmonicOrder,3)>0
            fprintf('kd%d=%4.3f, kp%d=%4.3f,  kw%d=%4.3f\n',HarmonicOrder,kd, HarmonicOrder,kp, HarmonicOrder, kw);
    end
end
 
PitchRatio=11/12;
fprintf('\n\nWinding Factors for a coil span of 11/12: \n');
for HarmonicOrder=1:2:13
    SlotAngle=(2*pi/(Qs/p));%  radians
    kd=sin(qs*HarmonicOrder*(SlotAngle*0.5)) / (  qs*sin(HarmonicOrder*(SlotAngle*0.5)) );
    kp=sin(HarmonicOrder*pi/2*PitchRatio);
    kw=kp*kd;
    if mod(HarmonicOrder,3)>0
            fprintf('kd%d=%4.3f, kp%d=%4.3f,  kw%d=%4.3f\n',HarmonicOrder,kd, HarmonicOrder,kp, HarmonicOrder, kw);
    end
end

PitchRatio=10/12;
fprintf('\n\nWinding Factors for a coil span of 10/12:  \n');

for HarmonicOrder=1:2:13
    SlotAngle=(2*pi/(Qs/p));%  radians
    kd=sin(qs*HarmonicOrder*(SlotAngle*0.5)) / (  qs*sin(HarmonicOrder*(SlotAngle*0.5)) );
    kp=sin(HarmonicOrder*pi/2*PitchRatio);
    kw=kp*kd;
    if mod(HarmonicOrder,3)>0
            fprintf('kd%d=%4.3f, kp%d=%4.3f,  kw%d=%4.3f\n',HarmonicOrder,kd, HarmonicOrder,kp, HarmonicOrder, kw);
    end
end

%%
% When compared to full-pitch coil span, 10/12  seems to provide optimum attenuation values for harmonics 
% the 5th harmonic and 7th harmonics are almost eliminated.
% 3rd harmonic and multiples will be eliminated on the line-to-line voltages.
%However, elimination of harmonics with distribution and  fractional pitch
%factor resulted in approximately 8% loss on the fundamental component
% This can be easily compansated increasing the number of stator turns or by increasing pole area
fprintf('---------------------------------------------------------\n');
fprintf('Selected coil span: 10/12  \n');
kd1=sin(qs*(SlotAngle*0.5)) / (  qs*sin(SlotAngle*0.5) );
kp1=sin(pi/2*PitchRatio);
kw1=kp1*kd1;
fprintf('---------------------------------------------------------\n');

%% Calculation of stator number of turns (Ns)
% Under the rated phase voltage, required number of series turns in a phase can be
% found using expected air gap flux and determined machine dimensions and winding factor:
% pole flux 2/pi*Bpeak*Ts*L
Flux_PerPole = 2*Di_s*L*Bgap/pole; % weber
Ns = Vphase/(4.44*fn*Flux_PerPole*kw1);
% number of conductorsin a phase:
Ns_cu=Ns*2;
fprintf('required number of series turns in a phase,Ns : %g\n', Ns);
% For a double layer winding
SemiSlotPerPhase=2*Qs/m;
% Then the number of turns per semi slot (SzQ) can be found by:
SzQ = round(Ns_cu/SemiSlotPerPhase);
zQs= SzQ*2;
fprintf('number of series conductors  in a slot,zQ: %g\n', zQs);
Ns_old=Ns;
Ns=zQs*p*qs;
fprintf('----updated Ns:%g\n', Ns);

% Resulting winding diagram is given in figure:
I = imread('Winding Diagram.png');
figure;
imshow(I);
title('Stator Winding Diagram','FontSize',16,'FontWeight','Bold');

% Although updated Ns is very close to calculated value, in order to keep air gap flux density at desired level,
% Machine dimensions must be scaled accordingly:
Di_s_new=sqrt(Di_s^2*X*Ns_old/(X*Ns));
L_new=X*Di_s_new;
Di_s=round(Di_s_new,3); % round up to 3 decimal places
L_s=round(L_new,3);

fprintf('Updated inner diameter of the stator is %g m\n',Di_s);
fprintf('Updated length of the machine is %g m\n',L);
Bgap_actual = Vphase*pole/(4.44*Ns*fn*kw1*4*L*Di_s*0.5); % Tesla
Bgap=Bgap_actual;% Update Bgap
Flux_PerPole = 2*Di_s*L*Bgap_actual/pole; % weber
Ns = Vphase/(4.44*fn*Flux_PerPole*kw1);
Ts=pi*Di_s/Qs; %Slot pitch
%for simplicity initially saturation factor taken as 1 

%% STATOR SLOT DIMENSIONS
I = imread('Stator Slot Dimensions.png');
figure;
imshow(I);
title('Stator Slot Dimensions','FontSize',16,'FontWeight','Bold');
% MMF calculation is important to check whether the air gap flux density

Js=5e6;  %current density for stator (A/m2)
Acu=Irated/Js;  %required copper cross section (mm2)
%Skin depth for Copper at %50 Hz 10mm, copper diameter must be smaller than
%20mm, considering slot opening:
dcu=2.91e-3;  %wire diameter (m) 
fprintf('wire diameter %g(mm), AWG9\n',dcu*1e3);
Strands= ceil(Acu/ (pi*dcu*dcu/4));
fprintf('Number of strands=%g\n',Strands);
Kf=0.3; %copper fill factor
Ass=Acu/Kf; % Slot area, mm2
fprintf('Slot Area=%g mm2\n',Ass*1e6);

Bt=1.5; % Average tooth flux From Table 6.2 (1.4T-2.1T)
b_os=3e-3; %m, slot opening typically selected as 3mm opening must allow wire insertion during manufacturing
h_os=1e-3;%m
h_w=1e-3;%m 
b_ts=((Bgap*Ts)/(Kfe*Bt)) ;% stator tooth width, bds,m
b_s1=(pi*(Di_s+2*(h_os+h_w))/Qs)-(b_ts); % b4 (book notation)
b_s2 = 1.2*b_s1; % from FEA model                                                                                                           !!!!!
h_s=2*Ass/(b_s1+b_s2);
fprintf('The height of slot=%g mm\n',h_s*1e3);

% Outer diameter of machine typically defined as 1.66*inner diameter for an
% 8 pole machine, however a more optimum selection can be obtained as
B_sbc=1.4;% Back core flux density
% For the calculation of the height of the stator back iron or yoke (hcs),
h_cs=Flux_PerPole/(2*L*B_sbc);% m
fprintf('The height of stator back iron is (h_cs) is %g mm\n',h_cs*1e3);
Do_s = (Di_s+2*(h_os+h_w+h_s+h_cs)); % m
Do_s=round(Do_s,3); % round up to 3 decimal places
fprintf('Outer diameter of the stator %g m\n',Do_s);

%%
peak_MMF=Bgap*g/u0;
Imag=peak_MMF/(4/pi*(m/2)*(Ns/pole)*kw1*sqrt(2)); % amps
fprintf('The peak MMF is %g Amps\n',Imag);
%% ROTOR SLOT DIMENSIONS
% Rotor Bar Current Calculations pg320 in [1]
I = imread('Rotor Slot Dimensions.png');
figure;
imshow(I);
title('Rotor Slot Dimensions','FontSize',16,'FontWeight','Bold');

Ir=zQs*Qs/Qr*Irated*pf_expected;
Iring=Ir/(2*sin(pi*p/Qr));
Jr=4*1e6;  %current density for rotor bar (A/m2)
Jring=4*1e6;  %current density for rotor bar (A/m2)
Ar=Ir/Jr; % Rotor Bar Crossection, m2
Aring=Iring/Jring; % Rotor Bar Crossection, m2
Tr=(pi*(Di_s-(2*g)))/Qr; % Rotor slot pitch, m
Btr=1.6 ;%T, Rotor tooth flux density
b_tr=Bgap*Tr/(Btr*Kfe) ;%Rotor tooth width to prevent rotor flux density from saturization

% Typical values assigned for slot opening and tooth hat
h_or = 2*1e-3; % m
b_or = 4*1e-3; % m
% upper diameter and lower diameter of the slot can be approximated as
d_r1=Tr-b_tr;
d_r2=0.5*d_r1;
%approximating tear-shaped bar as a triangle with a safety factor of 1.2, Ar= (h_r*(2*d_r1)/2*1.2
h_r=Ar/1.2/d_r1;
% Flux path in the rotor follows a similar behaviour as in the stator
% Then, calculation of rotor back core height (h_cr):
B_rbc=1.2;% Rotor back core flux density
h_cr = Flux_PerPole/(2*L*B_rbc); % m
fprintf('The height of the rotor back iron (hcr) is %g mm.\n',1e3*h_cr);
Di_r=Di_s-2*(g+h_r+h_cr+d_r1+d_r2);
%% todo Shaft Diameter








close all;
