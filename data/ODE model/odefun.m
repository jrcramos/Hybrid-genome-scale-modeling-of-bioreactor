function f=odefun(t,state,pars)%------------------------------------------

state =max(0,state);
X       =state(1) ; 
mAB       =state(2); 
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
V        =state(47);


%Reaction rates
[VHK, VPGI, VPFK, VPGK, VPK, VLDH, VG6PDH, VEP, VTK, VPDH, VCS, VCITS, VAKGDH, VSDH, VMLD, VPC, VME, VGLNT, VGLDH, VALATA, VASN, VASTA, VAATOSUC, VHISARGTA, VGLUT, VSDHH, VATPASE, VNADPHOX, VRESP, VLEAK, VAK, VPPRIBP, MU, VMAB, KdX,NADbasal,NADPbasal,AMPbasal,THDR,VPYRT, Feed]=kinetics(X,mAB,ALA,ARG, ASN, ASP, CYS, EGLC, EGLN, EGLU,EPYR, GLY, HIS, ILE, LAC, LEU, LYS, MET, NH4, PHE, PRO, SER, THR, TYR, VAL, ACCOA, ADP, AKG, AMP, ATP, CIT, F6P, G6P, GAP, GLU, MAL, NAD, NADH, NADP, NADPH, OXA, PEP, PYR, R5P, SUC, X5P, V,t,pars);


%Material balances
% define feeding strategy
%           X    mAB   ALA     ARG     ASN     ASP     CYS    EGLC      EGLN    EGLU  EPYR  GLY    HIS    ILE    LAC     LEU      LYS      MET   NH4  PHE     PRO    SER      THR    TYR    VAL      ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P   GAP     GLU      MAL       NAD      NADH     NADP       NADPH  OXA      PEP       PYR     R5P      SUC   X5P      V
feed_conc= [0    0       0     29.71   0       0       3.32   42.8     23.53    0.00  39.67  0.00   11.73  46.96  0      27.71    1.02     3.57  0    23.89   41.31  72.97    32.12  0.0    56.21    0          0        0         0         0       0         0        0     0       0        0         0        0         0          0     0        0         0       0        0    0       0];
%           X    mAB   ALA       ARG   ASN    ASP    CYS    EGLC    EGLN   EGLU EPYR   GLY    HIS    ILE   LAC    LEU   LYS   MET   NH4  PHE   PRO   SER    THR  TYR  VAL||  ACCOA     ADP       AKG       AMP       ATP     CIT      F6P      G6P   GAP     GLU      MAL       NAD      NADH     NADP       NADPH  OXA      PEP       PYR     R5P      SUC  X5P   ||| V
feed_flag= [0     0     0      1      1      1       1     1       1     1     1      1      1      1     0      1     1     1     0     1     1      1    1    1    1      0        0         0         0          0       0        0        0     0       0        0          0        0        0           0      0        0        0       0         0    0        0 ];
Feed_mass=Feed*feed_conc.*feed_flag; % mmol
D=Feed/V;


%Extracellular
f=zeros(47,1);
f(1,1)= (MU-KdX)*X                                                                         -D*X    ;%X 
f(2,1)= (VMAB)*X*1e3                                                                       -D*mAB  ;%mAB 
f(3,1)= (VALATA   -5.4e-4*VMAB      -9.70e-5*MU)*X*1e3 + Feed_mass(3)/V                    -D*ALA  ;%ALA
f(4,1)= (-VHISARGTA  -1.94e-4*VMAB   -5.34e-5*MU)*X*1e3  + Feed_mass(4)/V                  -D*ARG  ;%ARG
f(5,1)= (-VASN  -3.6e-4*VMAB        -4.51e-5*MU)*X*1e3    + Feed_mass(5)/V                 -D*ASN  ;%ASN   
f(6,1)= (VASN - VASTA -2*VPPRIBP  -3.18e-4*VMAB -5.57e-5*MU)*X*1e3 + Feed_mass(6)/V        -D*ASP  ;%ASP
f(7,1)= (-2.21e-4*VMAB -1.8e-5*MU)*X*1e3                             + Feed_mass(7)/V      -D*CYS  ;%CYS 
f(8,1)= ( -VHK)*X*1e3                                                   + Feed_mass(8)/V   -D*EGLC ;%EGLC 
f(9,1)= (-VGLNT -2*VPPRIBP  -4.15e-4*VMAB -7.63e-5*MU)*X*1e3       + Feed_mass(9)/V        -D*EGLN ;%EGLN 
f(10,1)= (VGLUT   -4.29e-4*VMAB -6.55e-5*MU)*X*1e3                 + Feed_mass(10)/V       -D*EGLU ;%EGLU
f(11,1)= (-VPYRT)*X*1e3                                            + Feed_mass(11)/V       -D*EPYR ;%EPYR
f(12,1)= (+THDR-VPPRIBP + 2*VPPRIBP  -6.09e-4*VMAB -1.03e-4*MU)*X*1e3   + Feed_mass(12)/V  -D*GLY  ;%GLY
f(13,1)= (-VHISARGTA   -1.80e-4*VMAB   -2.08e-5*MU)*X*1e3          + Feed_mass(13)/V       -D*HIS  ;%HIS
f(14,1)= (-VAATOSUC  -1.94e-4*VMAB -4.76e-5*MU)*X*1e3              + Feed_mass(14)/V       -D*ILE  ;%ILE
f(15,1)= (VLDH)*X*1e3                                                                      -D*LAC  ;%LAC  
f(16,1)= (-VAATOSUC -6.23e-4*VMAB -9.23e-5*MU)*X*1e3               + Feed_mass(16)/V       -D*LEU  ;%LEU
f(17,1)= (-VAATOSUC  -6.78e-4*VMAB -8.35e-5*MU)*X*1e3              + Feed_mass(17)/V       -D*LYS  ;%LYS
f(18,1)=  ( -8.30e-5*VMAB  -2.28e-5*MU)*X*1e3                      + Feed_mass(18)/V       -D*MET  ;%MET 
f(19,1)= (VGLNT + VGLDH + VASN + VHISARGTA + VSDHH)*X*1e3                                  -D*NH4  ;%NH4
f(20,1)=   (-2.77e-4*VMAB -4.04e-5*MU)*X*1e3                       + Feed_mass(20)/V       -D*PHE  ;%PHE 
f(21,1)=  (-6.64e-4*VMAB -4.80e-5*MU)*X*1e3                        + Feed_mass(21)/V       -D*PRO  ;%PRO
f(22,1)=  (-VSDHH   -1.22e-3*VMAB -6.98e-5*MU)*X*1e3               + Feed_mass(22)/V       -D*SER  ;%SER  
f(23,1)=   ( -THDR-7.61e-4*VMAB  -5.52e-5*MU)*X*1e3                + Feed_mass(23)/V       -D*THR  ;%THR 
f(24,1)= (-VAATOSUC  -4.29e-4*VMAB -2.84e-5*MU)*X*1e3              + Feed_mass(24)/V       -D*TYR  ;%TYR 
f(25,1)= (-VAATOSUC  -9.16e-4*VMAB -6.97e-5*MU)*X*1e3              + Feed_mass(25)/V       -D*VAL  ;%VAL  

%Intracellular
PO=2.5;
f(26,1)= VPDH - VCS + 8*VAATOSUC+ THDR   - MU*ACCOA                                                                                               ;%ACCOA
f(27,1)= VHK + VPFK - VPGK - VPK - VSDH + VPC - VGLNT + VAATOSUC + VATPASE - 2*PO*VRESP + 2*VAK + 2*VPPRIBP + 6.53e-3*MU + 2.75e-2*VMAB - MU*ADP  ;%ADP 
f(28,1)= VCITS - VAKGDH + VGLDH + VALATA - VASTA -7*VAATOSUC - VHISARGTA - MU*AKG                                                                 ;%AKG
f(29,1)= - VAK + VPPRIBP+AMPbasal - MU*AMP                                                                                                        ;%AMP  
f(30,1)= -VHK - VPFK + VPGK + VPK + VSDH - VPC + VGLNT - VAATOSUC - VATPASE + 2*PO*VRESP - VAK -2*VPPRIBP -6.53e-3*MU -2.75e-2*VMAB    -MU*ATP    ;%ATP
f(31,1)= VCS - VCITS - 8.64e-5*MU -MU*CIT                                                                                                         ;%CIT 
f(32,1)= VPGI - VPFK + 2*VTK -MU*F6P                                                                                                              ;%F6P
f(33,1)= VHK - VPGI - VG6PDH -5.5e-6*MU -MU*G6P                                                                                                   ;%G6P
f(34,1)= 2*VPFK - VPGK + VTK  -MU*GAP                                                                                                             ;%GAP
f(35,1)= VGLNT - VGLDH - VALATA + VASTA + 4*VAATOSUC + 4*VHISARGTA - VGLUT -MU*GLU                                                                ;%GLU
f(36,1)= VSDH - VMLD -VME + VAATOSUC + 2*VPPRIBP -MU*MAL                                                                                          ;%MAL
f(37,1)= - VPGK + VLDH - VPDH - VCITS - VAKGDH - 2*VSDH/3 - VMLD - VME - VGLDH -9*VAATOSUC + 2*VRESP + 2*VLEAK+NADbasal-THDR -MU*NAD              ;%NAD
f(38,1)= VPGK - VLDH + VPDH + VCITS + VAKGDH + 2*VSDH/3 + VMLD + VME + VGLDH + 9*VAATOSUC - 2*VRESP - 2*VLEAK+THDR -MU*NADH                       ;%NADH
f(39,1)= - 2*VG6PDH - 2*VAATOSUC + VNADPHOX+NADPbasal -MU*NADP                                                                                    ;%NADP
f(40,1)= 2*VG6PDH + 2*VAATOSUC - VNADPHOX -MU*NADPH                                                                                               ;%NADPH
f(41,1)= - VCS + VMLD + VPC + VASTA -MU*OXA                                                                                                       ;%OXA
f(42,1)= VPGK - VPK -MU*PEP                                                                                                                       ;%PEP
f(43,1)= VPK - VLDH - VPDH - VPC + VME - VALATA + VSDHH+VPYRT -MU*PYR                                                                             ;%PYR
f(44,1)= VG6PDH - VEP - VTK - 0.6*VPPRIBP -2.63e-5*MU -MU*R5P                                                                                     ;%R5P
f(45,1)= VAKGDH - VSDH + 3*VAATOSUC -MU*SUC                                                                                                       ;%SUC
f(46,1)= VEP - 2*VTK -MU*X5P                                                                                                                      ;%X5P
f(47,1)=Feed; %V
end%-----------------------------------------------------------------------
%------------------------------------------------------------

% kinetic model CHO Mab----------------------------------------
function [VHK, VPGI, VPFK, VPGK, VPK, VLDH, VG6PDH, VEP, VTK, VPDH, VCS, VCITS, VAKGDH, VSDH, VMLD, VPC, VME, VGLNT, VGLDH, VALATA, VASN, VASTA, VAATOSUC, VHISARGTA, VGLUT, VSDHH, VATPASE, VNADPHOX, VRESP, VLEAK, VAK, VPPRIBP, MU, VMAB, KdX, NADbasal,NADPbasal,AMPbasal,THDR,VPYRT, Feed]=kinetics(X,mAB,ALA,ARG, ASN, ASP, CYS, EGLC, EGLN, EGLU, EPYR, GLY, HIS, ILE, LAC, LEU, LYS, MET, NH4, PHE, PRO, SER, THR, TYR, VAL, ACCOA, ADP, AKG, AMP, ATP, CIT, F6P, G6P, GAP, GLU, MAL, NAD, NADH, NADP, NADPH, OXA, PEP, PYR, R5P, SUC, X5P, V, t,pars)

aAMP = pars(1);
aF6P = pars(2);
bAMP = pars(3);
bF6P = pars(4);
KAAMP = pars(5);
KAF6P = pars(6); 
KdCIT = pars(7); 
KdEGLN = pars(8); %mM
KdG6P = pars(9); 
KdLAC = pars(10); %mM 
KdPEP = pars(11);
KdPYR = pars(12);
KdLACg = pars(13);%mM
KdNH4g = pars(14);%mM
KmACCOA = pars(15);
KmADP = pars(16);
KmAKG = pars(17);
KmALA = pars(18); %mM
KmAMP = pars(19);
KmARG = pars(20);%mM
KmASP = pars(21);%mM
KmASN = pars(22);%mM
KmATP = pars(23);
KmCIT = pars(24);
KmEGLC = pars(25);%mM
KmEGLN = pars(26);%mM
KmEGLU = pars(27);%mM

KmF6P = pars(28);
KmG6P = pars(29);
KmGAP = pars(30);
KmGLU = pars(31);
KmGLY = pars(32);
KmHIS = pars(33);%mM
KmILE = pars(34);%mM
KmLAC = pars(35);%mM
KmLEU = pars(36); %mM
KmLYS = pars(37);%mM
KmMAL = pars(38);
KmNADH = pars(39);
KmNADPH = pars(40);
KmNH4 = pars(41);%mM
KmOXA = pars(42);
KmPEP = pars(43);
KmPYR = pars(44);
KmR5P = pars(45);
KmSER = pars(46); %mM
KmSUC = pars(47); 
KmTYR = pars(48);%mM
KmVAL = pars(49);%mM
KmX5P = pars(50); 


KdALA = pars(51); %mM

Pratio = pars(52);
VfmaxAK = pars(53);%
VrmaxAK = pars(54);%
VmaxAKGDH = pars(55);%
VfmaxALATA = pars(56);%
VrmaxALATA = pars(57);%
VmaxASN = pars(58);%
VmaxASTA = pars(59);%
VmaxATPASE = pars(60);%
VmaxCITS = pars(61);%
VmaxCS = pars(62);%
VmaxEP = pars(63);%
VmaxG6PDH = pars(64);%
VfmaxGLDH = pars(65);%
VrmaxGLDH = pars(66);%
VfmaxGLNT = pars(67);%
VrmaxGLNT = pars(68);%
VfmaxGLUT = pars(69);%
VrmaxGLUT = pars(70);%

Vmaxgrowth = pars(71); %h-1

VmaxHISARGTA = pars(72); %
VmaxHK = pars(73);%
VfmaxLDH = pars(74);%
VrmaxLDH = pars(75);%%
Vmaxleak = pars(76);%
VmaxAATOSUC = pars(77);%
VmaxMAB = pars(78);%
VmaxME = pars(79);%
VmaxMLD = pars(80);%
VmaxNADPHOX = pars(81);%
VmaxPC = pars(82);%
VmaxPDH = pars(83);%
VmaxPFK = pars(84);%
VfmaxPGI = pars(85);%
VrmaxPGI = pars(86);%
VmaxPGK = pars(87);%
VmaxPK = pars(88);%
VmaxPPRIBP = pars(89);%
Vmaxresp = pars(90);%
VmaxSDH = pars(91);%
VmaxSDHH = pars(92);%
VmaxTK = pars(93);%
KmALAg = pars(94);%mM
KmARGg = pars(95);%mM
KmASNg = pars(96); %mM
KmASPg = pars(97); %mM
KmATPg = pars(98);%mM
KmCITg = pars(99); %
KmCYSg = pars(100);%mM
KmEGLNg = pars(101);%mM
KmG6Pg = pars(102); %
KmEGLUg = pars(103); %mM
KmGLYg = pars(104);%mM
KmHISg = pars(105);%mM
KmILEg = pars(106);%mM
KmLYSg = pars(107);%mM
KmMETg = pars(108);%mM
KmPHEg = pars(109);%mM
KmPROg = pars(110);%mM
KmR5Pg = pars(111);%
KmSERg = pars(112);%mM
KmTHRg = pars(113);%mM
KmTYRg = pars(114);%mM
KmVALg = pars(115);%mM
KmALAMAB = pars(116);%mM
KmARGMAB = pars(117);%mM
KmASNMAB = pars(118);%mM
KmASPMAB = pars(119);%mM
KmATPMAB = pars(120);%mM
KmCYSMAB = pars(121);%mM
KmEGLNMAB = pars(122);%mM
KmEGLUMAB = pars(123);%
KmGLYMAB = pars(124);%mM
KmHISMAB = pars(125);%mM
KmILEMAB = pars(126);%mM
KmLEUMAB = pars(127); %mM
KmLYSMAB = pars(128);%mM
KmMETMAB = pars(129);%mM
KmPHEMAB = pars(130);%mM
KmPROMAB = pars(131);%mM
KmSERMAB = pars(132);%mM
KmTHRMAB = pars(133);%mM
KmTYRMAB = pars(134);%mM
KmVALMAB = pars(135);%mM
KmADPATP = pars(136);
KmATPADP = pars(137);
KmNADNADH = pars(138);
KmNADHNAD = pars(139);
KmNADPNADPH = pars(140);

kdMin= pars(141); %h-1
kdMax= pars(142); %h-1
beta = pars(143); %h-1
Tind= pars(144); % h
F0= pars(145); %Lh-1xcell
FPI= pars(146); % Lh-1
alfa= pars(147); % -
K_bNAD= pars(148);
keq_NAD= pars(149);%-
K_bAMP= pars(150);
keq_AMP= pars(151);
K_bNADP= pars(152);
keq_NADP= pars(153);%-
keq_AK= pars(154);%-
km_AK= pars(155);
VrmaxTHDR= pars(156);%
km_THDR= pars(157);
Km_TYR_VAA= pars(158);
vrmaxPYRT= pars(159);%
k_EPYR= pars(160);
k_PYR= pars(161);

e1=EGLC * (1+(bAMP *AMP/ATP)/(aAMP*KAAMP));
e2=KmEGLC * (1+(AMP/ATP)/(KAAMP));
e3=EGLC * (1+(AMP/ATP)/(aAMP*KAAMP));
e4=(ATP/ADP)/(KmATPADP+(ATP/ADP));
e5=EGLC/(KmEGLC+EGLC);
e6=KdG6P/(KdG6P+G6P);
VHK = VmaxHK * (e1/(e2+e3))*e4*e5*e6; %closed. %

VPGI = VfmaxPGI * (G6P/(KmG6P+G6P)) * (KdPEP/(KdPEP+PEP)) - VrmaxPGI * (F6P/(KmF6P+F6P)); %closed. %

e1=F6P * (1+(bAMP *AMP/ATP)/(aAMP*KAAMP));
e2=KmF6P * (1+(AMP/ATP)/(KAAMP));
e3=F6P * (1+(AMP/ATP)/(aAMP*KAAMP));
e4=(ATP/ADP)/(KmATPADP+(ATP/ADP));
e5=KdLAC/(KdLAC+LAC);
VPFK = VmaxPFK * (e1/(e2+e3))*e4*e5; %closed.% 

VPGK = VmaxPGK * (GAP/(KmGAP+GAP)) * ((ADP/ATP)/(KmADPATP+(ADP/ATP))) * ((NAD/NADH)/(KmNADNADH+(NAD/NADH))); %closed. %

e1=PEP * (1+(bF6P *F6P)/(aF6P*KAF6P));
e2=KmPEP * (1+F6P/(KAF6P));
e3=F6P * (1+F6P/KAF6P);
e4=(ADP/ATP)/(KmADPATP+(ADP/ATP));
e5=KdALA/(KdALA+ALA);
VPK = VmaxPK * (e1/(e2+e3))*e4*e5; %closed.% 
% VPK = VmaxPK * PEP/(1e-8+PEP);%

e1=PYR * (1+(bAMP *AMP/ATP)/(aAMP*KAAMP));
e2=KmPYR * (1+(AMP/ATP)/(KAAMP));
e3=PYR * (1+(AMP/ATP)/(aAMP*KAAMP));
e4=(NADH/NAD)/(KmNADHNAD+(NADH/NAD));
e5=LAC/(KmLAC+LAC);
e6=(NAD/NADH)/(KmNADNADH+(NAD/NADH));
e7=KdPYR/(KdPYR+PYR);
VLDH = VfmaxLDH * (e1/(e2+e3))*e4- VrmaxLDH*e5*e6*e7; %closed. %

VG6PDH = VmaxG6PDH * G6P / (KmG6P + G6P) * (NADP/NADPH) / (KmNADPNADPH + (NADP/NADPH)); %closed.  %

VEP = VmaxEP * R5P / (KmR5P + R5P); %closed.  %

VTK = VmaxTK * R5P / (KmR5P + R5P) * X5P / (KmX5P + X5P); %closed.  %

VPDH = VmaxPDH * PYR / (KmPYR + PYR) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %

VCS = VmaxCS * OXA / (KmOXA + OXA) * ACCOA / (KmACCOA + ACCOA);%closed.  %

VCITS = VmaxCITS * CIT / (KmCIT + CIT) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %

VAKGDH = VmaxAKGDH * AKG / (KmAKG + AKG) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %

VSDH = VmaxSDH * SUC / (KmSUC + SUC) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) * (ADP/ATP) / (KmADPATP + (ADP/ATP));%closed.  %

VMLD = VmaxMLD * MAL / (KmMAL + MAL) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %

VPC = VmaxPC * PYR / (KmPYR +PYR) * (ADP/ATP) / (KmADPATP + (ADP/ATP));%closed.  %

VME = VmaxME * MAL / (KmMAL + MAL) * (NAD/NADH) / (KmNADNADH + (NAD/NADH));%closed.  %

VGLNT = VfmaxGLNT*0.1 * EGLN / (KmEGLN + EGLN)* (ATP/ADP) / (KmATPADP + (ATP/ADP)) - VrmaxGLNT * GLU / (KmGLU + GLU) * (ADP/ATP) / (KmADPATP + (ADP/ATP)) * NH4/(KmNH4 + NH4); %closed.  %
% VGLNT      = 0.01*VfmaxGLNT * EGLN / (KmEGLN + EGLN);

% VGLDH = VfmaxGLDH * GLU / (KmGLU + GLU) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) - VrmaxGLDH * AKG / (KmAKG + AKG) * (NADH/NAD) / (KmNADHNAD + (NADH/NAD)) * NH4/(KmNH4 + NH4); %closed.  %
VGLDH = VfmaxGLDH * GLU / (KmGLU + GLU) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) - VrmaxGLDH * AKG / (KmAKG + AKG) * (NADH/NAD) / (KmNADHNAD + (NADH/NAD)) * KmNH4/(KmNH4 + NH4); %closed.  %

VALATA = VfmaxALATA * GLU / (KmGLU + GLU) * PYR / (KmPYR +PYR) - VrmaxALATA * ALA / (KmALA + ALA) * AKG / (KmAKG + AKG) * KdEGLN / (KdEGLN + EGLN);%closed.  %

VASN = VmaxASN * ASN / (KmASN + ASN); %closed.  %

VASTA = VmaxASTA * AKG / (KmAKG + AKG) * ASP / (KmASP + ASP); %closed.  %

% VAATOSUC = VmaxAATOSUC * LYS /(KmLYS + LYS) * ILE / (KmILE + ILE) * AKG / (KmAKG + AKG) * VAL / (KmVAL + VAL) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) * (NADP/NADPH) / (KmNADPNADPH + (NADP/NADPH)) * (ATP/ADP) / (KmATPADP + (ATP/ADP));%closed.  %
VAATOSUC = 0.05*VmaxAATOSUC * LYS /(KmLYS + LYS) * ILE / (KmILE + ILE) * AKG / (KmAKG + AKG) * VAL / (KmVAL + VAL) * (NAD/NADH) / (KmNADNADH + (NAD/NADH)) * (NADP/NADPH) / (KmNADPNADPH + (NADP/NADPH)) * (ATP/ADP) / (KmATPADP + (ATP/ADP))*TYR / (Km_TYR_VAA + TYR);%closed.  %

VHISARGTA = VmaxHISARGTA * HIS / (KmHIS + HIS) * ARG /(KmARG + ARG) * LEU / (KmLEU + LEU) * AKG / (KmAKG + AKG); %closed.  %
% VHISARGTA = VmaxHISARGTA * HIS / (KmHIS + HIS) * ARG /(KmARG + ARG) * LEU / (KmLEU + LEU) * AKG / (KmAKG + AKG)* KmNH4/(KmNH4 + NH4); %closed.  %

VGLUT = VfmaxGLUT * GLU / (KmGLU + GLU) - VrmaxGLUT * EGLU / (KmEGLU + EGLU);%closed.  %

VSDHH = VmaxSDHH * SER / (KmSER + SER); %closed.  %
% VSDHH = VmaxSDHH * SER / (KmSER + SER)*kiNH4/(kiNH4+NH4); %closed.  %

VATPASE = VmaxATPASE * ATP / (KmATP + ATP);%closed.  %

VNADPHOX = VmaxNADPHOX * NADPH / (KmNADPH + NADPH); %closed.  %

VRESP = Vmaxresp * (ADP/ATP) / (KmADPATP + (ADP/ATP)) * NADH / (KmNADH + NADH);%closed.  %

VLEAK = Vmaxleak * NADH / (KmNADH + NADH); %closed.  %

% VAK = VfmaxAK * AMP / (KmAMP + AMP) * ATP / (KmATP + ATP) - VrmaxAK * ADP / (KmADP + ADP); %closed.  %
VAK=1e3*VfmaxAK*(AMP-ADP/keq_AK)/(km_AK+AMP+ADP/keq_AK);

VPPRIBP = VmaxPPRIBP * R5P / (KmR5P + R5P) * EGLN / (KmEGLN + EGLN) * ASP / (KmASP + ASP) * GLY / (KmGLY + GLY) * ATP / (KmATP + ATP); %closed.  %

MU = Vmaxgrowth * (G6P / (KmG6Pg + G6P)) * (R5P / (KmR5Pg + R5P)) * (LYS / (KmLYSg + LYS)) * (CIT / (KmCITg + CIT)) * (ILE / (KmILEg + ILE)) * (VAL / (KmVALg + VAL)) * (TYR / (KmTYRg + TYR)) * (ASN / (KmASNg + ASN)) * (ASP / (KmASPg + ASP)) * (ALA / (KmALAg + ALA)) * (ARG / (KmARGg + ARG)) * (EGLN / (KmEGLNg + EGLN)) * (EGLU / (KmEGLUg + EGLU)) * (GLY / (KmGLYg + GLY)) * (HIS / (KmHISg + HIS)) * (PHE / (KmPHEg + PHE)) * (PRO / (KmPROg + PRO)) * (THR / (KmTHRg + THR)) * (KdLACg / (KdLACg + LAC)) * (KdNH4g / (KdNH4g + NH4))* (ATP / (KmATPg + ATP)); %closed.  % 
% MU_ =  Vmaxgrowth * (LYS / (KmLYSg + LYS)) * (CIT / (KmCITg + CIT))* (ILE / (KmILEg + ILE)) * (VAL / (KmVALg + VAL)) * (TYR / (KmTYRg + TYR)) * (ASN / (KmASNg + ASN)) * (ASP / (KmASPg + ASP)) * (ALA / (KmALAg + ALA))* (ARG / (KmARGg + ARG)) * (EGLN / (KmEGLNg + EGLN)) * (EGLU / (KmEGLUg + EGLU)) * (GLY / (KmGLYg + GLY)) * (HIS / (KmHISg + HIS)) * (PHE / (KmPHEg + PHE)) * (PRO / (KmPROg + PRO)) * (THR / (KmTHRg + THR))* (ATP / (KmATPg + ATP))    *EGLC/(EGLC+0.98); 
if t>Tind
    VMAB = VmaxMAB * (LYS / (KmLYSMAB + LYS)) * (ILE / (KmILEMAB + ILE)) * (VAL / (KmVALMAB + VAL)) * (TYR / (KmTYRMAB + TYR)) * (ASN / (KmASNMAB + ASN)) * (ASP / (KmASPMAB + ASP)) * (ALA / (KmALAMAB + ALA)) * (ARG / (KmARGMAB + ARG)) * (EGLN / (KmEGLNMAB + EGLN)) * (EGLU / (KmEGLUMAB + EGLU)) * (GLY / (KmGLYMAB + GLY)) * (HIS / (KmHISMAB + HIS)) * (PHE / (KmPHEMAB + PHE)) * (PRO / (KmPROMAB + PRO)) * (THR / (KmTHRMAB + THR)) * (ATP / (KmATPMAB + ATP)); %closed.  %
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

upars=pars(end-20:end);
nupars=size(upars,1);
tf=240;
dt=tf/nupars;
i = min(floor(t/dt)+1,nupars);
Feed = upars(i); %feedrate

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
