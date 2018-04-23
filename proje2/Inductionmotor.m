% Input Parameters

clear all;
CosPhi=0.85;                                                                %
u0=4*pi*1e-7;
HarmonicOrder=1;
RatedVoltage=400;                                                           %l-l
RatedPhaseVoltage=RatedVoltage/(sqrt(3));                                   % star connected
NumberofPhases=3;
Frequency=50;
ExpectedEfficiency=0.8;
PolePairs=2;
Poles=2*PolePairs;
SyncSpeed=50/PolePairs;                                                     %rps
RatedSlip=0.01;
RatedSpeed=(1-RatedSlip)*SyncSpeed;
NumberofSlots=36;
CopperResistance=0.021;                                                     % Ohm at 80C
RotorDiameter=55e-3;                                                        % meter
StatorOuterDiameter=90e-3;
StatorTeethWidth=2.25e-3;
StatorTeethHeight=12e-3;
StatorCoreDepth=0.5*(StatorOuterDiameter-RotorDiameter-2*StatorTeethHeight);

ToothPerPole=NumberofSlots/(PolePairs*2)
IronDensity=7.8e3;

SlotPerPolePerPhase=(NumberofSlots/(2*PolePairs*NumberofPhases));
PitchRatio=8/9;
SlotAngle=(2*pi/(NumberofSlots/PolePairs))
kd=sin(SlotPerPolePerPhase*HarmonicOrder*(SlotAngle*0.5)) / (    SlotPerPolePerPhase*sin(HarmonicOrder*(SlotAngle*0.5)) );
kp=sin(HarmonicOrder*pi/2*PitchRatio);
kw=kp*kd;

StatorSlotArea=41e-6;
CurrentDensity=4e6;                                                         %A/m2
FillFactor=0.5;
CarterFactor=1.05;
%% 

MachineConstant=80;                                                         %kWs/m3
AspectRatio= PolePairs^(1/3)*pi/(2*PolePairs);%
EffectiveLength= AspectRatio*RotorDiameter;                                 % meter
ExpectedMecPower=MachineConstant*RotorDiameter^2*EffectiveLength*SyncSpeed;

RatedRmsCurrent=(ExpectedMecPower*1e3)/(RatedVoltage*1.73*ExpectedEfficiency);% Amper
WireSize= RatedRmsCurrent/CurrentDensity;                                   % m2
NumberofConductorsinaSlot= floor(StatorSlotArea*FillFactor/(WireSize));
NumberofConductorsinaSlot= NumberofConductorsinaSlot-mod(NumberofConductorsinaSlot,2);
NumberofTurnsinaPhase=NumberofConductorsinaSlot*PolePairs*SlotPerPolePerPhase;
Airgap= (0.18+0.006*(ExpectedMecPower*1000)^0.4)/1000;                             % meter
EffectiveAirgap=Airgap*CarterFactor;

ActualFlux=RatedPhaseVoltage/(4.44*Frequency*NumberofTurnsinaPhase*kw);
Bgav=ActualFlux/(pi*RotorDiameter*EffectiveLength/Poles);
Bgpeak=Bgav*pi/2;
ElectricalLoading=NumberofConductorsinaSlot*RatedRmsCurrent*NumberofSlots/(pi*RotorDiameter);
RatedTorque=(ExpectedMecPower*1e3)/(2*pi*RatedSpeed)
MeanTurnLength= 2*EffectiveLength+2.3*pi*RotorDiameter*PitchRatio/(2*PolePairs)+0.08;
StatorResistance=MeanTurnLength*NumberofTurnsinaPhase*CopperResistance;

PhaseInductance=(2/pi)*(pi*RotorDiameter/(2*PolePairs))*EffectiveLength*u0/(Airgap*CarterFactor)*(4*SlotPerPolePerPhase/pi)*(NumberofPhases/NumberofSlots)*(kw*NumberofTurnsinaPhase)^2; 
MagnetizingInductance=PhaseInductance*(NumberofPhases/2);
LeakageInductance=MagnetizingInductance*0.011;

StatorCopperLoss=3*StatorResistance*RatedRmsCurrent^2;

MeanTeethFluxDensity= ActualFlux/(StatorTeethWidth*EffectiveLength*ToothPerPole);
PeakTeethFluxDensity=pi/2*MeanTeethFluxDensity;
TeethLossFactor=35;                                                         %W/kg
StatorTeethVolume=EffectiveLength*StatorTeethWidth*NumberofSlots*StatorTeethHeight;
TeethMass=StatorTeethVolume*IronDensity;
TeethLoss=TeethMass*TeethLossFactor;

StatorCrossSection=EffectiveLength*StatorCoreDepth;
StatorCoreFluxDensity=ActualFlux/(2*StatorCrossSection);
StatorVolume=(0.25*StatorOuterDiameter^2-0.25*(StatorOuterDiameter-StatorCoreDepth)^2)*EffectiveLength;
StatorLossFactor=35;
StatorMass=StatorVolume*IronDensity;


