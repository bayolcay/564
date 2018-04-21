
CosPhi=0.85;                                                                %
RatedVoltage=400;                                                           %l-l
NumberofPhases=3;
NumberofSlots=36;
Efficiency=0.8;
PolePairs=2;
StatorSlotArea=41e-6;
CurrentDensity=3e6;                                                        %A/m2
FillFactor=0.6;
SlotPerPolePerPhase=(NumberofSlots/(2*PolePairs*NumberofPhases))

MachineConstant=80;                                                        %kWs/m3
SyncSpeed=50/PolePairs;                                                    %rps
AspectRatio= PolePairs^(1/3)*pi/(2*PolePairs);%
RotorDiameter=55e-3;                                                        %
EffectiveLength= AspectRatio*RotorDiameter;
ExpectedMecPower=MachineConstant*RotorDiameter^2*EffectiveLength*SyncSpeed;
RatedRmsCurrent=ExpectedMecPower*1e3/(RatedVoltage*1.73*Efficiency);        % Amper
WireSize= RatedRmsCurrent/CurrentDensity;                                   %m2
ConductorCountinaSlot=StatorSlotArea*FillFactor/(WireSize);




SlotAngle=(2*pi/(NumberofSlots/PolePairs))
kd=sin(SlotPerPolePerPhase*HarmonicOrder*(SlotAngle*0.5))  / (    SlotPerPolePerPhase*sin(HarmonicOrder*(SlotAngle*0.5)) )