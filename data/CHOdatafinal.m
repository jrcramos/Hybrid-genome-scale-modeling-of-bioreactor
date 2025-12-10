%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHETIC DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nbatch, batch]=CHOdatafinal
rng default;
    
  nbatch=9;
  % Central Composite design of Experiences DoE
  % The controllers are u1(Desired miu pre induction) between 0-0.25 and
  % u2(Feeding rate post induction) between 0-0.025
  %Induction starts at 96h
  mat=ccdesign(2); %9 experiments
  miupreDoE=rescale(mat(1:9,1),0,0.01); %Rescale matrix elements between 0 and 0.01 (desired miu) Feed=miupre*X0*V0*exp(miupre*t)
  FeedpostDoE=rescale(mat(1:9,2),0,0.03); %Rescale matrix elements between 0 and 0.03 (Feed rate post induction)
  for i=1:nbatch
    tind=96;
    miupre=miupreDoE(i);
    Feedpost=FeedpostDoE(i);
    X0=0.2;
    V0=2;
    f=[];
    for t=1:8 %pre induction feed rate
        f(t,1)=(miupre*X0*V0*exp(miupre*t));
    end
    for t=9:21
        f(t,1)=Feedpost;
    end
    f=[f];
    
    batch(i).id=sprintf('batch%u',i);
    
    batch(i).u=f; %controller parameters

    doplot=1;
    [t,species,vol,raterule,ualongtime,sdspecies,species_true,rvol_true,rann_true]=...
        protrue(batch(i).u,doplot);
    
    batch(i).np=           length(t);
    batch(i).t=            t;
    
    %species
    batch(i).cnoise=       species;
    %batch(i).cnoise=       species_true; %TESTING
    batch(i).cnoise(1,:)=    species_true(1,:); %true initial values
 
    %compartment size
    batch(i).vol(:,1)=          vol;
    
    %raterules
    batch(i).raterule=     raterule;
    
    %state
    %batch(i).state =  [species, vol, raterule];
    batch(i).state =  [species, vol];
    %batch(i).state =  [species_true, vol]; %TEST ONLY
    batch(i).state(1,:)=    [species_true(1,:),vol(1)]; %true initial values
   
    batch(i).sc=      [sdspecies, ones(batch(i).np,1)];
    
    %control actions
    batch(i).ualongtime=   ualongtime;
    
    %true process data only important for plot purpose
    % if not available, assign null
    batch(i).c_true=       species_true;
    batch(i).rvol_true=    rvol_true;
    batch(i).rann_true=    rann_true; 
    
    
    tspan=batch(i).t;
    x=batch(i).c_true;
    x(:,47)=batch(i).u;
%     
%    Name=sprintf('Batch %i-1',i);
%    f=figure('Name',Name);
%     
% movegui(f,'north');
% subplot(5,5,1), plot(tspan,x(:,1)), ylabel('X')
% subplot(5,5,2), plot(tspan,x(:,2)), ylabel('mAB')
% subplot(5,5,3), plot(tspan,x(:,3)), ylabel('ALA')
% subplot(5,5,4), plot(tspan,x(:,4)), ylabel('ARG')
% subplot(5,5,5), plot(tspan,x(:,5)), ylabel('ASN')
% subplot(5,5,6), plot(tspan,x(:,6)), ylabel('ASP')
% subplot(5,5,7), plot(tspan,x(:,7)), ylabel('CYS')
% subplot(5,5,8), plot(tspan,x(:,8)), ylabel('EGLC')
% subplot(5,5,9), plot(tspan,x(:,9)), ylabel('EGLN')
% subplot(5,5,10), plot(tspan,x(:,10)), ylabel('EGLU')
% subplot(5,5,11), plot(tspan,x(:,11)), ylabel('GLY')
% subplot(5,5,12), plot(tspan,x(:,12)), ylabel('HIS')
% subplot(5,5,13), plot(tspan,x(:,13)), ylabel('ILE')
% subplot(5,5,14), plot(tspan,x(:,14)), ylabel('LAC')
% subplot(5,5,15), plot(tspan,x(:,15)), ylabel('LEU')
% subplot(5,5,16), plot(tspan,x(:,16)), ylabel('LYS')
% subplot(5,5,17), plot(tspan,x(:,17)), ylabel('MET')
% subplot(5,5,18), plot(tspan,x(:,18)), ylabel('NH4')
% subplot(5,5,19), plot(tspan,x(:,19)), ylabel('PHE')
% subplot(5,5,20), plot(tspan,x(:,20)), ylabel('PRO')
% subplot(5,5,21), plot(tspan,x(:,21)), ylabel('SER')
% subplot(5,5,22), plot(tspan,x(:,22)), ylabel('THR')
% subplot(5,5,23), plot(tspan,x(:,23)), ylabel('TYR')
% subplot(5,5,24), plot(tspan,x(:,24)), ylabel('VAL')
% 
%    Name=sprintf('Batch %i-2',i);
%    f=figure('Name',Name);
% movegui(f,'east');
% subplot(5,5,1), plot(tspan,x(:,25)), ylabel('ACCOA')
% subplot(5,5,2), plot(tspan,x(:,26)), ylabel('ADP')
% subplot(5,5,3), plot(tspan,x(:,27)), ylabel('AKG')
% subplot(5,5,4), plot(tspan,x(:,28)), ylabel('AMP')
% subplot(5,5,5), plot(tspan,x(:,29)), ylabel('ATP')
% subplot(5,5,6), plot(tspan,x(:,30)), ylabel('CIT')
% subplot(5,5,7), plot(tspan,x(:,31)), ylabel('F6P')
% subplot(5,5,8), plot(tspan,x(:,32)), ylabel('G6P')
% subplot(5,5,9), plot(tspan,x(:,33)), ylabel('GAP')
% subplot(5,5,10), plot(tspan,x(:,34)), ylabel('GLU')
% subplot(5,5,11), plot(tspan,x(:,35)), ylabel('MAL')
% subplot(5,5,12), plot(tspan,x(:,36)), ylabel('NAD')
% subplot(5,5,13), plot(tspan,x(:,37)), ylabel('NADH')
% subplot(5,5,14), plot(tspan,x(:,38)), ylabel('NADP')
% subplot(5,5,15), plot(tspan,x(:,39)), ylabel('NADPH')
% subplot(5,5,16), plot(tspan,x(:,40)), ylabel('OXA')
% subplot(5,5,17), plot(tspan,x(:,41)), ylabel('PEP')
% subplot(5,5,18), plot(tspan,x(:,42)), ylabel('PYR')
% subplot(5,5,19), plot(tspan,x(:,43)), ylabel('R5P')
% subplot(5,5,20), plot(tspan,x(:,44)), ylabel('SUC')
% subplot(5,5,21), plot(tspan,x(:,45)), ylabel('X5P')
% subplot(5,5,22), plot(tspan,x(:,46)), ylabel('Feed')
    

end

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [t,cnoise,vol,raterule,ucontrol, sc, c, rvol_true,rann_true]=protrue(upars,doplot)
tspan=0:12:240;

%    X    mAB   ALA    ARG    ASN     ASP    CYS    EGLC    EGLN   EGLU  EPYR   GLY    HIS    ILE   LAC    LEU    LYS   MET   NH4   PHE   PRO   SER    THR   TYR   VAL    ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P    GAP     GLU      MAL       NAD      NADH     NADP       NADPH    OXA      PEP      PYR    R5P      SUC  X5P   ||| V
x0= [4    0    1.11    6.06   11.17   4.47   2.5    72.55   6.21   0.97  0.23   0.65   1.50   8.17  5.22   10.84  5.75  2.23  2.12  3.42  8.98  10.70  5.20  2.32  7.62   4.4e-3    6.7e-5    1.9e-7     3.1e-7     7e-6    6e-4   4.2e-8   5.6e-8  4e-7    2.2e-4   1.5e-6    8.9e-6   7.5e-7   9.5e-7    9.3e-7   2e-7    8.8e-7    5e-7  2.4e-8    8e-8  3E-9     0.270];

std=[1.2  70    1     1.5   0.6     0.25    0.2    90        1    0.12   2      1.1    0.5    2.5    1.1    1.3    0.3   0.2   2.1   1.1   2     3.3    1.2   0.12   2.7    1e-8       1e-8     1e-8       1e-8    1e-8    1e-8     1e-8    1e-8    1e-8     1e-8       1e-8    1e-8      1e-8     1e-8     1e-8    1e-8       1e-8   1e-8     1e-8  1e-8   0.1];
odeopts=odeset('AbsTol',1e-8,'RelTol',1e-7,'NonNegative',1:47);
ofun=@(t,s)odefun(t,s,upars);
[t,x]=ode15s(ofun,tspan,x0,odeopts);
np=length(t);
% concentrations
c=x(:,1:25);
%volume
vol=x(:,47);
%standard deviation of concentrations
sc=repmat(std(1:25),np,1); 
% add noise to concentrations
cnoise=zeros(np,25);
for i=1:25
   cnoise(:,i)=c(:,i)+randn(np,1)*std(i);
end
% raterules
raterule = [];
% kinetic rates
rann_true =zeros(np,3); %unknown kinetics
rvol_true=zeros(np,4);
ucontrol=zeros(np,2);
if doplot
    f=figure(1);
    movegui(f,'north');
    subplot(5,5,1), plot(tspan,x(:,1)); ylabel('X'); hold on
    subplot(5,5,2), plot(tspan,x(:,2)); ylabel('mAB'); hold on
    subplot(5,5,3), plot(tspan,x(:,3)); ylabel('ALA'); hold on
    subplot(5,5,4), plot(tspan,x(:,4)); ylabel('ARG'); hold on
    subplot(5,5,5), plot(tspan,x(:,5)); ylabel('ASN'); hold on
    subplot(5,5,6), plot(tspan,x(:,6)); ylabel('ASP'); hold on
    subplot(5,5,7), plot(tspan,x(:,7)); ylabel('CYS'); hold on
    subplot(5,5,8), plot(tspan,x(:,8)); ylabel('EGLC'); hold on
    subplot(5,5,9), plot(tspan,x(:,9)); ylabel('EGLN'); hold on
    subplot(5,5,10), plot(tspan,x(:,10)); ylabel('EGLU'); hold on
    subplot(5,5,11), plot(tspan,x(:,11)); ylabel('EPYR'); hold on
    subplot(5,5,12), plot(tspan,x(:,12)); ylabel('GLY'); hold on
    subplot(5,5,13), plot(tspan,x(:,13)); ylabel('HIS'); hold on
    subplot(5,5,14), plot(tspan,x(:,14)); ylabel('ILE'); hold on
    subplot(5,5,15), plot(tspan,x(:,15)); ylabel('LAC'); hold on
    subplot(5,5,16), plot(tspan,x(:,16)); ylabel('LEU'); hold on
    subplot(5,5,17), plot(tspan,x(:,17)); ylabel('LYS'); hold on
    subplot(5,5,18), plot(tspan,x(:,18)); ylabel('MET'); hold on
    subplot(5,5,19), plot(tspan,x(:,19)); ylabel('NH4'); hold on
    subplot(5,5,20), plot(tspan,x(:,20)); ylabel('PHE'); hold on
    subplot(5,5,21), plot(tspan,x(:,21)); ylabel('PRO'); hold on
    subplot(5,5,22), plot(tspan,x(:,22)); ylabel('SER'); hold on
    subplot(5,5,23), plot(tspan,x(:,23)); ylabel('THR'); hold on
    subplot(5,5,24), plot(tspan,x(:,24)); ylabel('TYR'); hold on
    subplot(5,5,25), plot(tspan,x(:,25)); ylabel('VAL'); hold on
    
    f=figure(2);
    movegui(f,'east');
    subplot(5,5,1), plot(tspan,x(:,26)); ylabel('ACCOA'); hold on
    subplot(5,5,2), plot(tspan,x(:,27)); ylabel('ADP'); hold on
    subplot(5,5,3), plot(tspan,x(:,28)); ylabel('AKG'); hold on
    subplot(5,5,4), plot(tspan,x(:,29)); ylabel('AMP'); hold on
    subplot(5,5,5), plot(tspan,x(:,30)); ylabel('ATP'); hold on
    subplot(5,5,6), plot(tspan,x(:,31)); ylabel('CIT'); hold on
    subplot(5,5,7), plot(tspan,x(:,32)); ylabel('F6P'); hold on
    subplot(5,5,8), plot(tspan,x(:,33)); ylabel('G6P'); hold on
    subplot(5,5,9), plot(tspan,x(:,34)); ylabel('GAP'); hold on
    subplot(5,5,10), plot(tspan,x(:,35)); ylabel('GLU'); hold on
    subplot(5,5,11), plot(tspan,x(:,36)); ylabel('MAL'); hold on
    subplot(5,5,12), plot(tspan,x(:,37)); ylabel('NAD'); hold on
    subplot(5,5,13), plot(tspan,x(:,38)); ylabel('NADH'); hold on
    subplot(5,5,14), plot(tspan,x(:,39)); ylabel('NADP'); hold on
    subplot(5,5,15), plot(tspan,x(:,40)); ylabel('NADPH'); hold on
    subplot(5,5,16), plot(tspan,x(:,41)); ylabel('OXA'); hold on
    subplot(5,5,17), plot(tspan,x(:,42)); ylabel('PEP'); hold on
    subplot(5,5,18), plot(tspan,x(:,43)); ylabel('PYR'); hold on
    subplot(5,5,19), plot(tspan,x(:,44)); ylabel('R5P'); hold on
    subplot(5,5,20), plot(tspan,x(:,45)); ylabel('SUC'); hold on
    subplot(5,5,21), plot(tspan,x(:,46)); ylabel('X5P'); hold on
    subplot(5,5,22), plot(tspan,x(:,47)); ylabel('Feed'); hold on
end
end
%-------------------------------------------------------------------------
% PARAMETERIZATION OF CONTROL VARIABLES
function [u]=control_function_CHO(t,upars)%------------------------------------
% define feeding strategy
%           X    mAB   ALA     ARG     ASN     ASP     CYS    EGLC      EGLN    EGLU  EPYR  GLY    HIS    ILE    LAC    LEU      LYS     MET   NH4   PHE     PRO    SER      THR    TYR    VAL     ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P   GAP     GLU      MAL       NAD      NADH     NADP       NADPH  OXA      PEP       PYR     R5P      SUC  X5P   ||| V
feed_conc= [0    0       0     29.71   0       0       3.32   1892.8    23.53   0.00  39.67  0.00   11.73   46.96  0      27.71    1.02     3.57  0    23.89   41.31  72.97    32.12  0.0    56.21    0     0     0     0     0    0      0    0    0     0     0    0   0    0    0      0    0    0     0    0    0       0];
%           X    mAB   ALA     ARG   ASN    ASP    CYS    EGLC    EGLN   EGLU EPYR   GLY    HIS    ILE   LAC    LEU   LYS   MET   NH4  PHE   PRO   SER    THR  TYR  VAL||  ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P   GAP     GLU      MAL       NAD      NADH     NADP       NADPH  OXA      PEP       PYR     R5P      SUC  X5P   ||| V
feed_flag= [0     0     0      1      1      1       1     1       1     1     1      1      1      1     0      1     1     1     0     1     1      1    1    1    1      0        0         0         0          0       0        0        0     0       0        0          0        0      0           0      0         0       0      0         0     0        0 ];

nupars=size(upars,1);
tf=240;
dt=tf/nupars;
i = min(floor(t/dt)+1,nupars);
Feed = upars(i); %feedrate
Induc=0;
if t>=96
    Induc=1;
end
u =[Feed;Induc];
end%-----------------------------------------------------------------------



function f=odefun(t,state,upars)%------------------------------------------
%Concentrations of Extracellular species (X in 10^6 cells/mL; mAB in mg/L;
%the rest in mM (1 mmol/mL)
state =max(0,state);
X       =state(1) ; %Conversion factor is 2.93x10-3 mmol.10-6 cells for biomass
mAB       =state(2); %Conversion factor is 9.17x10-3 mmol.mg-1 for mAB
ALA      =state(3);
ARG      =state(4);
ASN      =state(5);
ASP      =state(6);
CYS      =state(7);
EGLC      =state(8);
EGLN      =state(9);
EGLU      =state(10);
EPYR      =state(11);
GLY      =state(12);
HIS      =state(13);
ILE      =state(14);
LAC      =state(15);
LEU      =state(16);
LYS      =state(17);
MET      =state(18);
NH4      =state(19);
PHE      =state(20);
PRO      =state(21);
SER      =state(22);
THR      =state(23);
TYR      =state(24);
VAL      =state(25);
%Concentrations of Extracellular species (all in mmol.10^-6 cells)
%state(25:45)=state(25:45)/1e6;
ACCOA      =state(26);
ADP      =state(27);
AKG      =state(28);
AMP      =state(29);
ATP      =state(30);
CIT      =state(31);
F6P      =state(32);
G6P      =state(33);
GAP      =state(34);
GLU      =state(35);
MAL      =state(36);
NAD      =state(37);
NADH      =state(38);
NADP      =state(39);
NADPH      =state(40);
OXA      =state(41);
PEP      =state(42);
PYR      =state(43);
R5P      =state(44);
SUC      =state(45);
X5P      =state(46);

%Compartment - Reactor Volume (L)
V        =state(47);

% Controllers NONE AT THE TIME
% [u]=control_function_parkramirez(t,upars);
% Feed =u(1);
% Sin = u(2);

%Reaction rates
[VHK, VPGI, VPFK, VPGK, VPK, VLDH, VG6PDH, VEP, VTK, VPDH, VCS, VCITS, VAKGDH, VSDH, VMLD, VPC, VME, VGLNT, VGLDH, VALATA, VASN, VASTA, VAATOSUC, VHISARGTA, VGLUT, VSDHH, VATPASE, VNADPHOX, VRESP, VLEAK, VAK, VPPRIBP, MU, VMAB, KdX,NADbasal,NADPbasal,AMPbasal,THDR,VPYRT, Feed]=kinetics(X,mAB,ALA,ARG, ASN, ASP, CYS, EGLC, EGLN, EGLU,EPYR, GLY, HIS, ILE, LAC, LEU, LYS, MET, NH4, PHE, PRO, SER, THR, TYR, VAL, ACCOA, ADP, AKG, AMP, ATP, CIT, F6P, G6P, GAP, GLU, MAL, NAD, NADH, NADP, NADPH, OXA, PEP, PYR, R5P, SUC, X5P, V,t,upars);

[u]=control_function_CHO(t,upars);
Feed=u(1);


%Material balances
% define feeding strategy
%           X    mAB   ALA     ARG     ASN     ASP     CYS    EGLC      EGLN    EGLU  EPYR  GLY    HIS    ILE    LAC    LEU      LYS     MET   NH4   PHE     PRO    SER      THR    TYR    VAL     ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P   GAP     GLU      MAL       NAD      NADH     NADP       NADPH  OXA      PEP       PYR     R5P      SUC  X5P   ||| V
feed_conc= [0    0       0     29.71   0       0       3.32   42.8    23.53   0.00  39.67  0.00   11.73   46.96  0      27.71    1.02     3.57  0    23.89   41.31  72.97    32.12  0.0    56.21    0     0     0     0     0    0      0    0    0     0     0    0   0    0    0      0    0    0     0    0    0       0];
%           X    mAB   ALA       ARG   ASN    ASP    CYS    EGLC    EGLN   EGLU EPYR   GLY    HIS    ILE   LAC    LEU   LYS   MET   NH4  PHE   PRO   SER    THR  TYR  VAL||  ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P   GAP     GLU      MAL       NAD      NADH     NADP       NADPH  OXA      PEP       PYR     R5P      SUC  X5P   ||| V
feed_flag= [0     0     0      1      1      1       1     1       1     1     1      1      1      1     0      1     1     1     0     1     1      1    1    1    1      0        0         0         0          0       0        0        0     0       0        0          0        0      0           0      0         0       0      0         0     0        0 ];
Feed_mass=Feed*feed_conc.*feed_flag; % mmol
D=Feed/V;

%Conversion factor is 9.17x10-3 mmol.mg-1 for mAB
%Conversion factor is 2.93x10-3 mmol.10-6 cells for biomass
%Extracellular
f(1,1)= (MU-KdX)*X                                                                    -D*X    ;%X 
f(2,1)= (VMAB)*X*1e3                                                                       -D*mAB  ;%mAB 
f(3,1)= (VALATA   -5.4e-4*VMAB      -9.70e-5*MU)*X*1e3 + Feed_mass(3)/V               -D*ALA  ;%ALA
f(4,1)= (-VHISARGTA  -1.94e-4*VMAB   -5.34e-5*MU)*X*1e3  + Feed_mass(4)/V             -D*ARG  ;%ARG
f(5,1)= (-VASN  -3.6e-4*VMAB        -4.51e-5*MU)*X*1e3    + Feed_mass(5)/V            -D*ASN  ;%ASN   
f(6,1)= (VASN - VASTA -2*VPPRIBP  -3.18e-4*VMAB -5.57e-5*MU)*X*1e3 + Feed_mass(6)/V   -D*ASP  ;%ASP
f(7,1)= (-2.21e-4*VMAB -1.8e-5*MU)*X*1e3                             + Feed_mass(7)/V   -D*CYS  ;%CYS 
f(8,1)= ( -VHK)*X*1e3                                                   + Feed_mass(8)/V   -D*EGLC ;%EGLC 
f(9,1)= (-VGLNT -2*VPPRIBP  -4.15e-4*VMAB -7.63e-5*MU)*X*1e3       + Feed_mass(9)/V   -D*EGLN ;%EGLN 
f(10,1)= (VGLUT   -4.29e-4*VMAB -6.55e-5*MU)*X*1e3                 + Feed_mass(10)/V  -D*EGLU ;%EGLU
f(11,1)= (-VPYRT)*X*1e3                                            + Feed_mass(11)/V  -D*EPYR ;%EPYR
f(12,1)= (+THDR-VPPRIBP + 2*VPPRIBP  -6.09e-4*VMAB -1.03e-4*MU)*X*1e3   + Feed_mass(12)/V  -D*GLY  ;%GLY
f(13,1)= (-VHISARGTA   -1.80e-4*VMAB   -2.08e-5*MU)*X*1e3          + Feed_mass(13)/V  -D*HIS  ;%HIS
f(14,1)= (-VAATOSUC  -1.94e-4*VMAB -4.76e-5*MU)*X*1e3              + Feed_mass(14)/V  -D*ILE  ;%ILE
f(15,1)= (VLDH)*X*1e3                                                                      -D*LAC  ;%LAC  
f(16,1)= (-VAATOSUC -6.23e-4*VMAB -9.23e-5*MU)*X*1e3               + Feed_mass(16)/V  -D*LEU  ;%LEU
f(17,1)= (-VAATOSUC  -6.78e-4*VMAB -8.35e-5*MU)*X*1e3              + Feed_mass(17)/V  -D*LYS  ;%LYS
f(18,1)=  ( -8.30e-5*VMAB  -2.28e-5*MU)*X*1e3                      + Feed_mass(18)/V  -D*MET  ;%MET 
f(19,1)= (VGLNT + VGLDH + VASN + VHISARGTA + VSDHH)*X*1e3                                  -D*NH4  ;%NH4
f(20,1)=   (-2.77e-4*VMAB -4.04e-5*MU)*X*1e3                       + Feed_mass(20)/V  -D*PHE  ;%PHE 
f(21,1)=  (-6.64e-4*VMAB -4.80e-5*MU)*X*1e3                        + Feed_mass(21)/V  -D*PRO  ;%PRO
f(22,1)=  (-VSDHH   -1.22e-3*VMAB -6.98e-5*MU)*X*1e3               + Feed_mass(22)/V  -D*SER  ;%SER  
f(23,1)=   ( -THDR-7.61e-4*VMAB  -5.52e-5*MU)*X*1e3                     + Feed_mass(23)/V  -D*THR  ;%THR 
f(24,1)= (-VAATOSUC  -4.29e-4*VMAB -2.84e-5*MU)*X*1e3              + Feed_mass(24)/V  -D*TYR  ;%TYR 
f(25,1)= (-VAATOSUC  -9.16e-4*VMAB -6.97e-5*MU)*X*1e3              + Feed_mass(25)/V  -D*VAL  ;%VAL  

%Intracellular
PO=2.5;
f(26,1)= VPDH - VCS + 8*VAATOSUC+ THDR   - MU*ACCOA                                                                                                         ;%ACCOA
f(27,1)= VHK + VPFK - VPGK - VPK - VSDH + VPC - VGLNT + VAATOSUC + VATPASE - 2*PO*VRESP + 2*VAK + 2*VPPRIBP + 6.53e-3*MU + 2.75e-2*VMAB - MU*ADP  ;%ADP 
f(28,1)= VCITS - VAKGDH + VGLDH + VALATA - VASTA -7*VAATOSUC - VHISARGTA - MU*AKG                                                                      ;%AKG
f(29,1)= - VAK + VPPRIBP+AMPbasal - MU*AMP                                                                                                                      ;%AMP  
f(30,1)= -VHK - VPFK + VPGK + VPK + VSDH - VPC + VGLNT - VAATOSUC - VATPASE + 2*PO*VRESP - VAK -2*VPPRIBP -6.53e-3*MU -2.75e-2*VMAB    -MU*ATP    ;%ATP
f(31,1)= VCS - VCITS - 8.64e-5*MU -MU*CIT                                                                                                         ;%CIT 
f(32,1)= VPGI - VPFK + 2*VTK -MU*F6P                                                                                                                   ;%F6P
f(33,1)= VHK - VPGI - VG6PDH -5.5e-6*MU -MU*G6P                                                                                                   ;%G6P
f(34,1)= 2*VPFK - VPGK + VTK  -MU*GAP                                                                                                                  ;%GAP
f(35,1)= VGLNT - VGLDH - VALATA + VASTA + 4*VAATOSUC + 4*VHISARGTA - VGLUT -MU*GLU                                                                     ;%GLU
f(36,1)= VSDH - VMLD -VME + VAATOSUC + 2*VPPRIBP -MU*MAL                                                                                               ;%MAL
f(37,1)= - VPGK + VLDH - VPDH - VCITS - VAKGDH - 2*VSDH/3 - VMLD - VME - VGLDH -9*VAATOSUC + 2*VRESP + 2*VLEAK+NADbasal-THDR -MU*NAD                                 ;%NAD
f(38,1)= VPGK - VLDH + VPDH + VCITS + VAKGDH + 2*VSDH/3 + VMLD + VME + VGLDH + 9*VAATOSUC - 2*VRESP - 2*VLEAK+THDR -MU*NADH                                 ;%NADH
f(39,1)= - 2*VG6PDH - 2*VAATOSUC + VNADPHOX+NADPbasal -MU*NADP                                                                                                   ;%NADP
f(40,1)= 2*VG6PDH + 2*VAATOSUC - VNADPHOX -MU*NADPH                                                                                                    ;%NADPH
f(41,1)= - VCS + VMLD + VPC + VASTA -MU*OXA                                                                                                            ;%OXA
f(42,1)= VPGK - VPK -MU*PEP                                                                                                                            ;%PEP
f(43,1)= VPK - VLDH - VPDH - VPC + VME - VALATA + VSDHH+VPYRT -MU*PYR                                                                                        ;%PYR
f(44,1)= VG6PDH - VEP - VTK - 0.6*VPPRIBP -2.63e-5*MU -MU*R5P                                                                                     ;%R5P
f(45,1)= VAKGDH - VSDH + 3*VAATOSUC -MU*SUC                                                                                                            ;%SUC
f(46,1)= VEP - 2*VTK -MU*X5P                                                                                                                           ;%X5P
%f(25:45,1)=f(25:45,1)*1e6;

% f(25:45,1)
% pause
%Volume
f(47,1)=Feed; %V
end%-----------------------------------------------------------------------
%------------------------------------------------------------

% kinetic model CHO Mab----------------------------------------
function [VHK, VPGI, VPFK, VPGK, VPK, VLDH, VG6PDH, VEP, VTK, VPDH, VCS, VCITS, VAKGDH, VSDH, VMLD, VPC, VME, VGLNT, VGLDH, VALATA, VASN, VASTA, VAATOSUC, VHISARGTA, VGLUT, VSDHH, VATPASE, VNADPHOX, VRESP, VLEAK, VAK, VPPRIBP, MU, VMAB, KdX, NADbasal,NADPbasal,AMPbasal,THDR,VPYRT, Feed]=kinetics(X,mAB,ALA,ARG, ASN, ASP, CYS, EGLC, EGLN, EGLU, EPYR, GLY, HIS, ILE, LAC, LEU, LYS, MET, NH4, PHE, PRO, SER, THR, TYR, VAL, ACCOA, ADP, AKG, AMP, ATP, CIT, F6P, G6P, GAP, GLU, MAL, NAD, NADH, NADP, NADPH, OXA, PEP, PYR, R5P, SUC, X5P, V, t,upars);

aAMP = 0.47;
aF6P = 4.1;
bAMP = 10.46;
bF6P = 1.7;
KAAMP = 0.09;
KAF6P = 4.9e-6; %mmol.10e-6cells
KdCIT = 7.4e-4; %mmol.10e-6cells
KdEGLN = 3.5; %mM
KdG6P = 1.710841e-8; %mmol.10e-6cells
KdLAC = 6.3; %mM 
KdPEP = 3.7e-7;%mmol.10e-6cells
KdPYR = 4.5e-7;%mmol.10e-6cells
KdLACg = 170;%mM
KdNH4g = 5.40;%mM
KmACCOA = 8.6e-7;%mmol.10e-6cells
KmADP = 7.2e-7;%mmol.10e-6cells
KmAKG = 1.1e-6;%mmol.10e-6cells
KmALA = 18; %mM
KmAMP = 2.4e-9;%mmol.10e-6cells
KmARG = 7.3e-2;%mM
KmASP = 4e-2;%mM
KmASN = 1.2;%mM
KmATP = 1.56e-5;%mmol.10e-6cells
KmCIT = 1.6e-7;%mmol.10e-6cells
KmEGLC = 20;%mM
KmEGLN = 2.3;%mM
KmEGLU = 6.6;%mM

KmF6P = 2.6e-6;%mmol.10e-6cells
KmG6P = 2.866e-8;%mmol.10e-6cells
KmGAP = 7.8e-7;%mmol.10e-6cells
KmGLU = 7.5e-5;%mmol.10e-6cells
KmGLY = 5e-2;%mmol.10e-6cells
KmHIS = 0.58;%mM
KmILE = 0.038;%mM
KmLAC = 2.9;%mM
KmLEU = 1e-6; %mM
KmLYS = 9.8e-2;%mM
KmMAL = 2.7e-6;%mmol.10e-6cells
KmNADH = 2.6e-6;%mmol.10e-6cells
KmNADPH = 8.4e-8;%mmol.10e-6cells
KmNH4 = 1.1;%mM
KmOXA = 2.2e-7;%mmol.10e-6cells
KmPEP = 1.7e-7;%mmol.10e-6cells
KmPYR = 3.9e-6;%mmol.10e-6cells
KmR5P = 7e-8;%mmol.10e-6cells
KmSER = 0.22; %mM
KmSUC = 2.3e-7; %mmol.10e-6cells
KmTYR = 4.8e-2;%mM
KmVAL = 1;%mM
KmX5P = 1.1e-8; %mmol.10e-6cells


KdALA = 0.1; %mM????

Pratio = 2.1;
VfmaxAK = 1.8e-8;%mmol.10-6cells.h-1
VrmaxAK = 1.2e-8;%mmol.10-6cells.h-1
VmaxAKGDH = 2.6e-4;%mmol.10-6cells.h-1
VfmaxALATA = 1.7e-4;%mmol.10-6cells.h-1
VrmaxALATA = 2.4e-4;%mmol.10-6cells.h-1
VmaxASN = 6.2e-6;%mmol.10-6cells.h-1
VmaxASTA = 7.8e-5;%mmol.10-6cells.h-1
VmaxATPASE = 9.14e-4;%mmol.10-6cells.h-1
VmaxCITS = 3.9e-5;%mmol.10-6cells.h-1
VmaxCS = 8.8e-5;%mmol.10-6cells.h-1
VmaxEP = 1.3e-5;%mmol.10-6cells.h-1
VmaxG6PDH = 1.5e-5;%mmol.10-6cells.h-1
VfmaxGLDH = 1.2e-6;%mmol.10-6cells.h-1
VrmaxGLDH = 6.3e-5;%mmol.10-6cells.h-1
VfmaxGLNT = 1.27e-4;%mmol.10-6cells.h-1
VrmaxGLNT = 1.9e-5;%mmol.10-6cells.h-1
VfmaxGLUT = 1.8e-6;%mmol.10-6cells.h-1
VrmaxGLUT = 2.5e-6;%mmol.10-6cells.h-1

Vmaxgrowth = 6.6e-2; %h-1

VmaxHISARGTA = 1.9e-5; %mmol.10-6cells.h-1
VmaxHK = 6.6e-4;%mmol.10-6cells.h-1
VfmaxLDH = 1.7e-3;%mmol.10-6cells.h-1
VrmaxLDH = 5.4e-4;%%mmol.10-6cells.h-1
Vmaxleak = 2.9e-5;%mmol.10-6cells.h-1
VmaxAATOSUC = 1.3e-4;%mmol.10-6cells.h-1
VmaxMAB = 6.2e-2;%mmol.10-6cells.h-1
VmaxME = 1.3e-5;%mmol.10-6cells.h-1
VmaxMLD = 9.5e-5;%mmol.10-6cells.h-1
VmaxNADPHOX = 1.4e-5;%mmol.10-6cells.h-1
VmaxPC = 2.1e-4;%mmol.10-6cells.h-1
VmaxPDH = 7.1e-4;%mmol.10-6cells.h-1
VmaxPFK = 1.5e-3;%mmol.10-6cells.h-1
VfmaxPGI = 7.6e-4;%mmol.10-6cells.h-1
VrmaxPGI = 3.5e-4;%mmol.10-6cells.h-1
VmaxPGK = 1.065e-3;%mmol.10-6cells.h-1
VmaxPK = 1.376e-3;%mmol.10-6cells.h-1
VmaxPPRIBP = 2.6e-9;%mmol.10-6cells.h-1
Vmaxresp = 0.0032;%mmol.10-6cells.h-1
VmaxSDH = 4e-4;%mmol.10-6cells.h-1
VmaxSDHH = 8.5e-6;%mmol.10-6cells.h-1
VmaxTK = 2.3e-5;%mmol.10-6cells.h-1
KmALAg = 7.6e-2;%mM
KmARGg = 1.5e-2;%mM
KmASNg = 4.2e-3; %mM
KmASPg = 1e-2; %mM
KmATPg = 1.4e-7;%mM
KmCITg = 3.5e-9; %mmol.10-6cells
KmCYSg = 1e-3;%mM
KmEGLNg = 1.4e-3;%mM
KmG6Pg = 9.9e-13; %mmol.10-6cells
KmEGLUg = 1.7e-6; %mM
KmGLYg = 1e-2;%mM
KmHISg = 1.7e-2;%mM
KmILEg = 1e-2;%mM
KmLYSg = 2.7e-2;%mM
KmMETg = 3.5e-2;%mM
KmPHEg = 1.6e-2;%mM
KmPROg = 1e-2;%mM
KmR5Pg = 1.1e-10;%mmol.10-6cells
KmSERg = 1.1e-3;%mM
KmTHRg = 9.6e-3;%mM
KmTYRg = 1.2e-2;%mM
KmVALg = 1.3e-3;%mM
KmALAMAB = 0.13;%mM
KmARGMAB = 6.2e-2;%mM
KmASNMAB = 1.4e-6;%mM
KmASPMAB = 5.6e-2;%mM
KmATPMAB = 1.1e-7;%mM
KmCYSMAB = 1e-3;%mM
KmEGLNMAB = 1.3e-3;%mM
KmEGLUMAB = 0.22;%mmol.10-6cells
KmGLYMAB = 5e-2;%mM
KmHISMAB = 5.7e-2;%mM
KmILEMAB = 8.3e-2;%mM
KmLEUMAB = 1e-6; %mM
KmLYSMAB = 5e-2;%mM
KmMETMAB = 8e-2;%mM
KmPHEMAB = 7.1e-2;%mM
KmPROMAB = 1.1e-2;%mM
KmSERMAB = 1e-3;%mM
KmTHRMAB = 6.6e-2;%mM
KmTYRMAB = 9.8e-2;%mM
KmVALMAB = 7.4e-2;%mM
KmADPATP = 0.569423;
KmATPADP = 0.091;
KmNADNADH = 4.8e-2;
KmNADHNAD = 0.3;
KmNADPNADPH = 2;

kdMin=2.6e-6; %h-1
kdMax=3.3e-3; %h-1
beta = 1e-3; %h-1
Tind= 96; % h
F0=6e-06; %Lh-1xcell
FPI= 0.0002; % Lh-1
alfa=0.066; % -
K_bNAD=1e-5;
keq_NAD=1e-5;%-
K_bAMP=1e-5;
keq_AMP=1e-4;
K_bNADP=1e-5;
keq_NADP=1e-5;%-
keq_AK=10;%-
km_AK=1e-6;
VrmaxTHDR=2e-7;%mmol.10-6cells.h-1
km_THDR=1;
Km_TYR_VAA=2;
vrmaxPYRT=3e-6;%mmol.10-6cells.h-1
k_EPYR=5e-3;
k_PYR=4e-7;


e1=EGLC * (1+(bAMP *AMP/ATP)/(aAMP*KAAMP));
e2=KmEGLC * (1+(AMP/ATP)/(KAAMP));
e3=EGLC * (1+(AMP/ATP)/(aAMP*KAAMP));
e4=(ATP/ADP)/(KmATPADP+(ATP/ADP));
e5=EGLC/(KmEGLC+EGLC);
e6=KdG6P/(KdG6P+G6P);
VHK = VmaxHK * (e1/(e2+e3))*e4*e5*e6; %closed. %mmol.10-6cells.h-1

VPGI = VfmaxPGI * (G6P/(KmG6P+G6P)) * (KdPEP/(KdPEP+PEP)) - VrmaxPGI * (F6P/(KmF6P+F6P)); %closed. %mmol.10-6cells.h-1

e1=F6P * (1+(bAMP *AMP/ATP)/(aAMP*KAAMP));
e2=KmF6P * (1+(AMP/ATP)/(KAAMP));
e3=F6P * (1+(AMP/ATP)/(aAMP*KAAMP));
e4=(ATP/ADP)/(KmATPADP+(ATP/ADP));
e5=KdLAC/(KdLAC+LAC);
VPFK = VmaxPFK * (e1/(e2+e3))*e4*e5; %closed.% mmol.10-6cells.h-1

VPGK = VmaxPGK * (GAP/(KmGAP+GAP)) * ((ADP/ATP)/(KmADPATP+(ADP/ATP))) * ((NAD/NADH)/(KmNADNADH+(NAD/NADH))); %closed. %mmol.10-6cells.h-1

e1=PEP * (1+(bF6P *F6P)/(aF6P*KAF6P));
e2=KmPEP * (1+F6P/(KAF6P));
e3=F6P * (1+F6P/KAF6P);
e4=(ADP/ATP)/(KmADPATP+(ADP/ATP));
e5=KdALA/(KdALA+ALA);
VPK = VmaxPK * (e1/(e2+e3))*e4*e5; %closed.% mmol.10-6cells.h-1
% VPK = VmaxPK * PEP/(1e-8+PEP);%

e1=PYR * (1+(bAMP *AMP/ATP)/(aAMP*KAAMP));
e2=KmPYR * (1+(AMP/ATP)/(KAAMP));
e3=PYR * (1+(AMP/ATP)/(aAMP*KAAMP));
e4=(NADH/NAD)/(KmNADHNAD+(NADH/NAD));
e5=LAC/(KmLAC+LAC);
e6=(NAD/NADH)/(KmNADNADH+(NAD/NADH));
e7=KdPYR/(KdPYR+PYR);
VLDH = VfmaxLDH * (e1/(e2+e3))*e4- VrmaxLDH*e5*e6*e7; %closed. %mmol.10-6cells.h-1

VG6PDH = VmaxG6PDH * G6P / (KmG6P + G6P) * (NADP/NADPH) / (KmNADPNADPH + (NADP/NADPH)); %closed.  %mmol.10-6cells.h-1

VEP = VmaxEP * R5P / (KmR5P + R5P); %closed.  %mmol.10-6cells.h-1

VTK = VmaxTK * R5P / (KmR5P + R5P) * X5P / (KmX5P + X5P); %closed.  %mmol.10-6cells.h-1

VPDH = VmaxPDH * PYR / (KmPYR + PYR) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %mmol.10-6cells.h-1

VCS = VmaxCS * OXA / (KmOXA + OXA) * ACCOA / (KmACCOA + ACCOA);%closed.  %mmol.10-6cells.h-1

VCITS = VmaxCITS * CIT / (KmCIT + CIT) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %mmol.10-6cells.h-1

VAKGDH = VmaxAKGDH * AKG / (KmAKG + AKG) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %mmol.10-6cells.h-1

VSDH = VmaxSDH * SUC / (KmSUC + SUC) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) * (ADP/ATP) / (KmADPATP + (ADP/ATP));%closed.  %mmol.10-6cells.h-1

VMLD = VmaxMLD * MAL / (KmMAL + MAL) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %mmol.10-6cells.h-1

VPC = VmaxPC * PYR / (KmPYR +PYR) * (ADP/ATP) / (KmADPATP + (ADP/ATP));%closed.  %mmol.10-6cells.h-1

VME = VmaxME * MAL / (KmMAL + MAL) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %mmol.10-6cells.h-1

VGLNT = VfmaxGLNT*0.1 * EGLN / (KmEGLN + EGLN)* (ATP/ADP) / (KmATPADP + (ATP/ADP)) - VrmaxGLNT * GLU / (KmGLU + GLU) * (ADP/ATP) / (KmADPATP + (ADP/ATP)) * NH4/(KmNH4 + NH4); %closed.  %mmol.10-6cells.h-1
% VGLNT      = 0.01*VfmaxGLNT * EGLN / (KmEGLN + EGLN);

% VGLDH = VfmaxGLDH * GLU / (KmGLU + GLU) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) - VrmaxGLDH * AKG / (KmAKG + AKG) * (NADH/NAD) / (KmNADHNAD + (NADH/NAD)) * NH4/(KmNH4 + NH4); %closed.  %mmol.10-6cells.h-1
VGLDH = VfmaxGLDH * GLU / (KmGLU + GLU) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) - VrmaxGLDH * AKG / (KmAKG + AKG) * (NADH/NAD) / (KmNADHNAD + (NADH/NAD)) * KmNH4/(KmNH4 + NH4); %closed.  %mmol.10-6cells.h-1

VALATA = VfmaxALATA * GLU / (KmGLU + GLU) * PYR / (KmPYR +PYR) - VrmaxALATA * ALA / (KmALA + ALA) * AKG / (KmAKG + AKG) * KdEGLN / (KdEGLN + EGLN);%closed.  %mmol.10-6cells.h-1

VASN = VmaxASN * ASN / (KmASN + ASN); %closed.  %mmol.10-6cells.h-1

VASTA = VmaxASTA * AKG / (KmAKG + AKG) * ASP / (KmASP + ASP); %closed.  %mmol.10-6cells.h-1

% VAATOSUC = VmaxAATOSUC * LYS /(KmLYS + LYS) * ILE / (KmILE + ILE) * AKG / (KmAKG + AKG) * VAL / (KmVAL + VAL) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) * (NADP/NADPH) / (KmNADPNADPH + (NADP/NADPH)) * (ATP/ADP) / (KmATPADP + (ATP/ADP));%closed.  %mmol.10-6cells.h-1
VAATOSUC = 0.05*VmaxAATOSUC * LYS /(KmLYS + LYS) * ILE / (KmILE + ILE) * AKG / (KmAKG + AKG) * VAL / (KmVAL + VAL) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) * (NADP/NADPH) / (KmNADPNADPH + (NADP/NADPH)) * (ATP/ADP) / (KmATPADP + (ATP/ADP))*TYR / (Km_TYR_VAA + TYR);%closed.  %mmol.10-6cells.h-1

VHISARGTA = VmaxHISARGTA * HIS / (KmHIS + HIS) * ARG /(KmARG + ARG) * LEU / (KmLEU + LEU) * AKG / (KmAKG + AKG); %closed.  %mmol.10-6cells.h-1
% VHISARGTA = VmaxHISARGTA * HIS / (KmHIS + HIS) * ARG /(KmARG + ARG) * LEU / (KmLEU + LEU) * AKG / (KmAKG + AKG)* KmNH4/(KmNH4 + NH4); %closed.  %mmol.10-6cells.h-1

VGLUT = VfmaxGLUT * GLU / (KmGLU + GLU) - VrmaxGLUT * EGLU / (KmEGLU + EGLU);%closed.  %mmol.10-6cells.h-1

VSDHH = VmaxSDHH * SER / (KmSER + SER); %closed.  %mmol.10-6cells.h-1
% VSDHH = VmaxSDHH * SER / (KmSER + SER)*kiNH4/(kiNH4+NH4); %closed.  %mmol.10-6cells.h-1

VATPASE = VmaxATPASE * ATP / (KmATP + ATP);%closed.  %mmol.10-6cells.h-1

VNADPHOX = VmaxNADPHOX * NADPH / (KmNADPH + NADPH); %closed.  %mmol.10-6cells.h-1

VRESP = Vmaxresp * (ADP/ATP) / (KmADPATP + (ADP/ATP)) * NADH / (KmNADH + NADH);%closed.  %mmol.10-6cells.h-1

VLEAK = Vmaxleak * NADH / (KmNADH + NADH); %closed.  %mmol.10-6cells.h-1

% VAK = VfmaxAK * AMP / (KmAMP + AMP) * ATP / (KmATP + ATP) - VrmaxAK * ADP / (KmADP + ADP); %closed.  %mmol.10-6cells.h-1
VAK=1e3*VfmaxAK*(AMP-ADP/keq_AK)/(km_AK+AMP+ADP/keq_AK);

VPPRIBP = VmaxPPRIBP * R5P / (KmR5P + R5P) * EGLN / (KmEGLN + EGLN) * ASP / (KmASP + ASP) * GLY / (KmGLY + GLY) * ATP / (KmATP + ATP); %closed.  %mmol.10-6cells.h-1

MU = Vmaxgrowth * (G6P / (KmG6Pg + G6P)) * (R5P / (KmR5Pg + R5P)) * (LYS / (KmLYSg + LYS)) * (CIT / (KmCITg + CIT)) * (ILE / (KmILEg + ILE)) * (VAL / (KmVALg + VAL)) * (TYR / (KmTYRg + TYR)) * (ASN / (KmASNg + ASN)) * (ASP / (KmASPg + ASP)) * (ALA / (KmALAg + ALA)) * (ARG / (KmARGg + ARG)) * (EGLN / (KmEGLNg + EGLN)) * (EGLU / (KmEGLUg + EGLU)) * (GLY / (KmGLYg + GLY)) * (HIS / (KmHISg + HIS)) * (PHE / (KmPHEg + PHE)) * (PRO / (KmPROg + PRO)) * (THR / (KmTHRg + THR)) * (KdLACg / (KdLACg + LAC)) * (KdNH4g / (KdNH4g + NH4))* (ATP / (KmATPg + ATP)); %closed.  %mmol.10-6cells.h-1 
% MU_ =  Vmaxgrowth * (LYS / (KmLYSg + LYS)) * (CIT / (KmCITg + CIT))* (ILE / (KmILEg + ILE)) * (VAL / (KmVALg + VAL)) * (TYR / (KmTYRg + TYR)) * (ASN / (KmASNg + ASN)) * (ASP / (KmASPg + ASP)) * (ALA / (KmALAg + ALA))* (ARG / (KmARGg + ARG)) * (EGLN / (KmEGLNg + EGLN)) * (EGLU / (KmEGLUg + EGLU)) * (GLY / (KmGLYg + GLY)) * (HIS / (KmHISg + HIS)) * (PHE / (KmPHEg + PHE)) * (PRO / (KmPROg + PRO)) * (THR / (KmTHRg + THR))* (ATP / (KmATPg + ATP))    *EGLC/(EGLC+0.98); 
if t>Tind
    VMAB = VmaxMAB * (LYS / (KmLYSMAB + LYS)) * (ILE / (KmILEMAB + ILE)) * (VAL / (KmVALMAB + VAL)) * (TYR / (KmTYRMAB + TYR)) * (ASN / (KmASNMAB + ASN)) * (ASP / (KmASPMAB + ASP)) * (ALA / (KmALAMAB + ALA)) * (ARG / (KmARGMAB + ARG)) * (EGLN / (KmEGLNMAB + EGLN)) * (EGLU / (KmEGLUMAB + EGLU)) * (GLY / (KmGLYMAB + GLY)) * (HIS / (KmHISMAB + HIS)) * (PHE / (KmPHEMAB + PHE)) * (PRO / (KmPROMAB + PRO)) * (THR / (KmTHRMAB + THR)) * (ATP / (KmATPMAB + ATP)); %closed.  %mmol.10-6cells.h-1
else
    VMAB=1e-12;
end
% f = (1-exp(-4*MU/((Vmaxgrowth))));
KdX=kdMin + kdMax*(beta/(beta+MU))^2;
NADbasal=K_bNAD*MU/Vmaxgrowth*(1-NAD/(keq_NAD));
NADPbasal=K_bNADP*MU/Vmaxgrowth*(1-NADP/(keq_NADP));
AMPbasal=K_bAMP*MU/Vmaxgrowth*(1-AMP/(keq_AMP));
THDR=VrmaxTHDR*THR/(km_THDR+km_THDR)*(NAD/NADH);
% VPYRT=vrmaxPYRT*(EPYR/k_EPYR)^2*(1-(PYR/EPYR/keq_PYRT))/(1+(EPYR/k_EPYR)^2+(PYR/k_PYR)^2);
VPYRT=vrmaxPYRT*EPYR/(1+EPYR+k_EPYR)*k_PYR/(k_PYR+PYR);
VHK=VHK*MU/Vmaxgrowth+(-0.071*0.5-0.076*0.5)/-0.364*6.8333e-05;

% Feed rate, 40-90h before induction exponential feed, >90h constant feed
% if t>30 %&& t<32 
%     Feed=F0*X;%6e-06*X;%1e-4;%;%* exp (alfa*(t-40));
% % elseif t>32
% %     Feed=FPI;
% else
%     Feed=0;
% end
Feed=[];
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [t,x]=ode_simulation_with_feed(ofun,tspan,x0,odeopts)
% % It is assumed that 6 feeds occur 70, 100, 130, 160, 190, 220h
% % Concentration on the feed (feed_conc) is assumed as the same as the initial media composition, x0
% % Feeding volumes (feed_vol) are 5 or 10 mL
% % feed_flag contains the info about which compounds were feeded
% % feed_conc, feed_flag, feed_vol are the same size as the number of model
% % state and are in the same order (ex: X, mAB, ALA, ...), making it easier to add the feeds to the corresponding ODE equation
% 
% %           X    mAB   ALA     ARG   ASN    ASP    CYS    EGLC    EGLN   EGLU   GLY    HIS    ILE  LAC    LEU   LYS   MET   NH4  PHE   PRO   SER    THR  TYR  VAL||  ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P   GAP     GLU      MAL       NAD      NADH     NADP       NADPH  OXA      PEP       PYR     R5P      SUC  X5P   ||| V
% feed_conc= [0.2  0     0.96    2.4   0.48    2     1.1    26      3.7    0.8    0.78   0.85   2.2  0.75   2.2   2.1   2.43  0.2  1.1   0.65  2.83   1.1  1.3  0.7    1.4e-7    6.7e-5    1e-7     3.1e-7     7e-6    1e-6    2.2e-7   1.6e-7  4e-7    2.2e-4   1.5e-6    6.9e-7   8.5e-6   3.5e-7    2.3e-7   1e-6    8.8e-7    2e-6  5.4e-8    8e-8  3E-7     2];
% %           X    mAB   ALA     ARG   ASN    ASP    CYS    EGLC    EGLN   EGLU   GLY    HIS    ILE  LAC    LEU   LYS   MET   NH4  PHE   PRO   SER    THR  TYR  VAL||  ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P   GAP     GLU      MAL       NAD      NADH     NADP       NADPH  OXA      PEP       PYR     R5P      SUC  X5P   ||| V
% feed_flag= [ 0     0     0     1      1      1       1     1       1     1      1      1      1     0      1     1     1     0     1     1      1    1    1    1      0        0         0         0          0       0        0        0     0       0        0          0        0      0           0      0         0        1      0         0     0        0 ];
% feed_vol=  [ 0     0     0     5      5      5       5     10      10    5      5      5      5     0      5     5     5     0     5     5      5    5    5    5      0        0         0         0          0       0        0        0     0       0        0          0        0      0           0      0         0        5      0         0     0        0 ];
% 
% fedd_mass=feed_conc.*feed_flag.*feed_vol*1e-3; % calculate the total mmol feeded for each model state
% 
% feed_time= [ 70 100 130 160 190 220]; % time that the feed occurs
% feed_time=feed_time(feed_time<max(tspan)); % if feed_times must be within tspan wanted
% 
% % tspan=sort([tspan feed_time
% %
% 
% 
% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    