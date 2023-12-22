clear all;
clc;

[DSSStartOK, DSSObj, DSSText] = DSSStartup;

if DSSStartOK
    

   

DSSObj.AllowForms = false;    

DSSText.Command = 'Set DataPath = C:\Users\rjha\Dropbox\2019\IEEE123SinglephaseValidation\IEEE123_1P_VR_C_DG10\withall4VR4C\IEEE123_1P_VR_C_OpenDSS_DG10';

DSSText.Command = 'compile IEEE123SinglePh.dss';    
    
%% modelling of capacitor for constant Q 
DSSText.Command = 'New generator.c83 bus1 = 83 phases=1 kv=2.4018 kw =0.0001 kvar =200';
DSSText.Command = 'New generator.c88 bus1 = 88 phases=1 kv=2.4018 kw =0.0001 kvar =16.667';
DSSText.Command = 'New generator.c90 bus1 = 90 phases=1 kv=2.4018 kw =0.0001 kvar =16.667';
DSSText.Command = 'New generator.c92 bus1 = 92 phases=1 kv=2.4018 kw =0.0001 kvar =16.667';

%% Voltage regulators with fix taps

ta1 = 1;   %% tap position of VR 1
ta2 = -5;     %% tap position of VR 5
ta3 = 4;     %% tap position of VR 4
ta4 = -4;    %% tap position of VR -4

% ta1 = 0;   %% tap position of VR 0
% ta2 = 0;     %% tap position of VR 0
% ta3 = 0;     %% tap position of VR 0
% ta4 = 0;    %% tap position of VR 0

DSSText.Command = 'New "Transformer.reg1a" windings=2 Xhl=0.001 %loadloss=1E-005 ppm_antifloat=0 phases=1 conns=[wye, wye, ] buses=[150, 150r, ] kVs=[2.4018, 2.4018, ] kVAs=[1666.7, 1666.7, ] normhkVA=1833.3 emerghkVA=2500';
DSSText.Command = 'New "RegControl.creg1a" transformer=reg1a winding=2 tapwinding=2 vreg=120 band=2 ptratio=20 CTprim=700 R=3 X=7.5';
fileID1=fopen('MyFile1.txt','w');
fprintf(fileID1,'Transformer.reg1a.wdg=2 Tap=(0.00625 %d * 1 +)',ta1);
formatSpec = '%c';
fileID1_1 = fopen('MyFile1.txt','r');
r_tap1 = fscanf(fileID1_1,formatSpec);
DSSText.Command = r_tap1 ; 


DSSText.Command = 'New "Transformer.reg2a" windings=2 bank=reg2 Xhl=0.01 %loadloss=1E-005 ppm_antifloat=0 phases=1 conns=[wye, wye, ] buses=[9, 9r, ] kVs=[2.402, 2.402, ] kVAs=[2000, 2000, ] normhkVA=2200 emerghkVA=3000';
DSSText.Command = 'New "RegControl.creg2a" transformer=reg2a winding=2 tapwinding=2 vreg=120 band=2 ptratio=20 CTprim=50 R=0.4 X=0.4';
fileID2=fopen('MyFile2.txt','w');
fprintf(fileID2,'Transformer.reg2a.wdg=2 Tap=(0.00625 %d * 1 +)',ta2);
formatSpec = '%c';
fileID2 = fopen('MyFile2.txt','r');
r_tap2 = fscanf(fileID2,formatSpec);
DSSText.Command = r_tap2 ; 

DSSText.Command = 'New "Transformer.reg3a" windings=2 bank=reg3 Xhl=0.01 %loadloss=1E-005 ppm_antifloat=0 phases=1 conns=[wye, wye, ] buses=[25, 25r, ] kVs=[2.402, 2.402, ] kVAs=[2000, 2000, ] normhkVA=2200 emerghkVA=3000';
DSSText.Command = 'New "RegControl.creg3a" transformer=reg3a winding=2 tapwinding=2 vreg=120 band=1 ptratio=20 CTprim=50 R=0.4 X=0.4';
fileID3=fopen('MyFile3.txt','w');
fprintf(fileID3,'Transformer.reg3a.wdg=2 Tap=(0.00625 %d * 1 +)',ta3);
formatSpec = '%c';
fileID3 = fopen('MyFile3.txt','r');
r_tap3 = fscanf(fileID3,formatSpec);
DSSText.Command = r_tap3 ;

DSSText.Command = 'New "Transformer.reg4a" windings=2 bank=reg4 Xhl=0.01 %loadloss=1E-005 ppm_antifloat=0 phases=1 conns=[wye, wye, ] buses=[160, 160r, ] kVs=[2.402, 2.402, ] kVAs=[2000, 2000, ] normhkVA=2200 emerghkVA=3000';
DSSText.Command = 'New "RegControl.creg4a" transformer=reg4a winding=2 tapwinding=2 vreg=124 band=2 ptratio=20 CTprim=300 R=0.6 X=1.3';
fileID4=fopen('MyFile4.txt','w');
fprintf(fileID4,'Transformer.reg4a.wdg=2 Tap=(0.00625 %d * 1 +)',ta4);
formatSpec = '%c';
fileID4 = fopen('MyFile4.txt','r');
r_tap4 = fscanf(fileID4,formatSpec);
DSSText.Command = r_tap4 ;

fclose all;
DSSText.Command='Set controlmode=OFF';

DSSText.Command =  'RegControl.creg1a.maxtapchange=0' ;  %Fixed at present tap
DSSText.Command =  'RegControl.creg2a.maxtapchange=0' ;  %Fixed at present tap
DSSText.Command =  'RegControl.creg3a.maxtapchange=0' ;  %Fixed at present tap
DSSText.Command =  'RegControl.creg4a.maxtapchange=0' ;  %Fixed at present tap


%% defining PV in opendss 
DSSText.Command = 'New PVSystem.9 phases=1 bus1=9 kV=2.4 kVA=16 irradiance=1 Pmpp=13.33 kVAr =0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.28 phases=1 bus1=28 kV=2.4 kVA=16 irradiance=1 Pmpp=13.33 kVAr =0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.48 phases=1 bus1=48 kV=2.4 kVA=84 irradiance=1 Pmpp=70 kVAr =0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.58 phases=1 bus1=58 kV=2.4 kVA=8 irradiance=1 Pmpp=6.667 kVAr =0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.65 phases=1 bus1=65 kV=2.4 kVA=56 irradiance=1 Pmpp=46.667 kVAr =0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.76 phases=1 bus1=76 kV=2.4 kVA=98 irradiance=1 Pmpp=81.666 kVAr =0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.87 phases=1 bus1=87 kV=2.4 kVA=16 irradiance=1 Pmpp=13.333 kVAr =0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.100 phases=1 bus1=100 kV=2.4 kVA=16 irradiance=1 Pmpp=13.333 kVAr =0 %cutin=0.1 %cutout=0.1';


DSSCircuit      = DSSObj.ActiveCircuit;
DSSSolution     = DSSCircuit.Solution;
DSSSolution.Solve;                          %% solving the circuit without PV
Buses=DSSCircuit.Allbusnames;               %% knowing all the buses in the system
DSSCircuit.SetActiveElement('Vsource.SOURCE');                  %% power generated by source
Power_subsopend = -1*(DSSCircuit.ActivecktElement.Powers);
Vallnodes = DSSCircuit.AllBusVmagPU;

% DSSText.Command = 'Show Voltage LN Nodes';
% DSSText.Command = 'Show Powers kVA Elem';
% DSSText.Command = 'Show taps';

%% Voltage as per ascending node number

load busnamematlab.txt;
busname = busnamematlab;

bus_voltage = [busname Vallnodes(2:129)'];         %% busname and corresponding volatge
Vasced = sort(busname);                     %% sorting the bus name in ascending order
VoltOpenDSS = [];
for i = 1:size(busname,1)
indx = find(i ==bus_voltage(:,1));
Vo = bus_voltage(indx,2);
VoltOpenDSS = [VoltOpenDSS Vo];
end

Voltagewithsubstaion = [Vallnodes(1) VoltOpenDSS];

          a = 'DSS start properly';
else
    
    a = 'DSS did not start properly';
    disp(a)
end

