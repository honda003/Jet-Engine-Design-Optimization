clc 
clear all
close all

%*******************************Design point*******************************

%Design Flight Condition
h=10000;
Mo=0.7;

%Constants
gc=1.4;
gt=1.33;
hpr=47000000;
Pi_d=0.97;
Pi_b=0.97;
Pi_n=0.97;
eita_b=0.97;
eita_m=0.99;
ec=0.9;
et=0.9;

R=287;
Cpc=(gc*R)/(gc-1);
Cpt=(R*gt)/(gt-1);

%Design Parameters
Pi_c=20;
Pi_ch=sqrt(20/0.95);
Pi_cl=0.95*Pi_ch;
Tt4=1000;
Weight=110000;
S=25;

%Intake 
[To, ao, Po, rho]=atmosisa(h);
rho_o=Po/(R*To);
vo=Mo*sqrt(gc*R*To);
taw_r=(1+(((gc-1)/2)*(Mo^2)));
Tto=To*taw_r;
Pi_r=((1+(((gc-1)/2)*(Mo^2)))^(gc/(gc-1)));
Pto=Po*Pi_r;

%Diffuser
Tt2=Tto;
Pt2=Pto*Pi_d;

%Compressor Low
taw_cl=(Pi_cl)^((gc-1)/(gc*ec));
Tt2half=Tt2*taw_cl;
Pt2half=Pi_cl*Pt2;
Work_CompressorLow=Cpc*(Tt2half-Tt2);

%Compressor High
taw_ch=(Pi_ch)^((gc-1)/(gc*ec));
Tt3=Tt2half*taw_ch;
Pt3=Pi_ch*Pt2half;
Work_CompressorHigh=Cpc*(Tt3-Tt2half);

%Burner
f=((Cpt*Tt4)-(Cpc*Tt3))/((hpr*eita_b)-(Cpt*Tt4));
Pt4=Pt3*Pi_b;
  
%Turbine High
syms Tt4req
eqn0=Cpc*(Tt3-Tt2half)-Cpt*eita_m*(1+f)*(Tt4-Tt4req);
Tt4half=vpasolve(eqn0==0,Tt4req);
taw_th=Tt4half/Tt4;
Pi_th=(taw_th)^(gt/((gt-1)*et));
Pt4half=Pt4*Pi_th;

%Turbine Low
syms Tt5req
eqn01=Cpc*(Tt2half-Tt2)-Cpt*eita_m*(1+f)*(Tt4half-Tt5req);
Tt5=vpasolve(eqn01==0,Tt5req);
taw_tl=Tt5/Tt4half;
Pi_tl=(taw_tl)^(gt/((gt-1)*et));
Pt5=Pt4half*Pi_tl;

%Nozzle
Tt9=Tt5;
Pt9_per_Pa=Pi_n*(Pi_tl*Pi_th)*Pi_b*(Pi_cl*Pi_ch)*Pi_d*Pi_r;

if Pt9_per_Pa <= 1.85
    
    P9=Po;
    Pt9=Pt9_per_Pa*P9;
    M9=sqrt((2/(gt-1))*((Pt9_per_Pa^((gt-1)/gt))-1));
    T9=Tt9/(1+(((gt-1)*(M9^2))/2));
    v9=sqrt(2*Cpt*Tt9*(1-((1/Pt9_per_Pa)^((gt-1)/gt))));
    Nozzle_Type=sprintf('To achieve Full expansion a Convergent Nozzle is required');
    
else
    
    P9=Po;
    Pt9=Pt9_per_Pa*P9;
    M9=sqrt((2/(gt-1))*((Pt9_per_Pa^((gt-1)/gt))-1));
    T9=Tt9/(1+(((gt-1)*(M9^2))/2));
    v9=sqrt(2*Cpt*Tt9*(1-((1/Pt9_per_Pa)^((gt-1)/gt))));
    Nozzle_Type=sprintf('To achieve Full expansion a Convergent-Divergent Nozzle is required');
end

%Nozzle Expansion Ratio
Mth=1;
E=(Mth/M9)*((1+((gt-1)*(M9^2)/2))/(1+((gt-1)*(Mth^2)/2)))^((gt+1)/((2*gt)-2));

%Parameters
F_per_mo=vpa(((1+f)*v9)-(vo));
Specific_FuelConsumption=vpa(f/F_per_mo);
Thermal_Efficiency=vpa((((1+f)*(v9^2))-((vo^2)))/(2*f*hpr));
Propulsive_Efficiency=vpa((2*F_per_mo*vo)/(((1+f)*(v9^2))-((vo^2))));
Overall_Efficiency=Thermal_Efficiency*Propulsive_Efficiency;

%Thrust Required
Lift=Weight;
CL=Lift/(0.5*rho_o*(vo^2)*S);
CD=0.022+(0.2*(CL^2));
Freq=CD*(0.5*rho_o*(vo^2)*S);

%Mass Flow Rate
mo=Freq/F_per_mo;

%Low Compressor Rotational Speed
syms Cceta2 

rm=0.5;                                   %mean radius
loading_coefficient=0.4;
eqn=Cceta2/loading_coefficient;

eqn1=vpa(Work_CompressorLow/Cceta2);
U1=sqrt(eqn*eqn1);
N1=vpa((U1*60)/(2*pi*rm));

%High Compressor Rotational Speed
eqn2=vpa(Work_CompressorHigh/Cceta2);
U2=sqrt(eqn*eqn2);
N2=vpa((U2*60)/(2*pi*rm));

%Engine Geometery 

Ao=mo/(rho_o*vo);
rho_9=P9/(R*T9);
A9=(mo*(1+f))/(rho_9*v9);

%%************************************2**************************************

%Off Design
%Chocked Nozzle with fixed area ratio
MFP8_rel=1;
Pi_th_rel=1;
taw_th_rel=1;
Pi_tl_rel=1;
taw_tl_rel=1;
MFP4_rel=1;

gc=1.4;
gt=1.33;
hpr=47000000;
Pi_d=0.97;
Pi_b=0.97;
Pi_n=0.97;
eita_b=0.97;
eita_m=0.99;
ec=0.9;
et=0.9;

c2=0;

syms xl xh

Tt4_off=Tt4; %Ref Parameter
Tt4_rel=1;

Tt4half_off=Tt4half;
Tt4half_rel=1;

h_off=0:1000:10000;
Mo_off=0.01:0.069:0.7;


syms xl xh M9_offs Cceta2


for i=1:length(Mo_off)
    c2=c2+1;
[To_off(c2), ao_off(c2), Po_off(c2), rho_off(c2)]=atmosisa(h_off(c2));

To_rel(c2)=To_off(c2)/To;
Po_rel(c2)=Po_off(c2)/Po;

vo_off(c2)=Mo_off(c2)*sqrt(gc*R*To_off(c2));
vo_rel(c2)=vo_off(c2)/vo;

taw_r_off(c2)=(1+(((gc-1)/2)*(Mo_off(c2)^2)));
Tto_off(c2)=To_off(c2)*taw_r_off(c2);

Pi_r_off(c2)=((1+(((gc-1)/2)*(Mo_off(c2)^2)))^(gc/(gc-1)));
Pi_r_rel(c2)=Pi_r_off(c2)/Pi_r;
Pto_off(c2)=Po_off(c2)*Pi_r_off(c2);

%Diffuser

Tt2_off(c2)=Tto_off(c2);
Pt2_off(c2)=Pto_off(c2)*Pi_d;

Pt2_rel(c2)=Pt2_off(c2)/Pt2;

Tt2_rel(c2)=Tt2_off(c2)/Tt2;
Tt4_rel=Tt4_off/Tt4;

Pi_tl_off=Pi_tl_rel*Pi_tl;
taw_tl_off=taw_tl_rel*taw_tl;

eqn5(c2)=(xl-1)-(((Cpt*Tt4half_off)*(1-taw_tl_off))/(Cpc*Tt2_off(c2)));
taw_cl_off(c2)=vpasolve(eqn5(c2)==0,xl);
taw_cl_rel(c2)=taw_cl_off(c2)/taw_cl;
Pi_cl_off(c2)=(taw_cl_off(c2))^((gc*ec)/(gc-1));
Pi_cl_rel(c2)=Pi_cl_off(c2)/Pi_cl;
Tt2half_off(c2)=taw_cl_off(c2)*Tt2_off(c2);

Pi_th_off=Pi_th_rel*Pi_th;
taw_th_off=taw_th_rel*taw_th;

eqn6(c2)=(xh-1)-(((Cpt*Tt4_off)*(1-taw_th_off))/(Cpc*Tt2half_off(c2)));
taw_ch_off(c2)=vpasolve(eqn6(c2)==0,xh);
taw_ch_rel(c2)=taw_ch_off(c2)/taw_ch;
Pi_ch_off(c2)=(taw_ch_off(c2))^((gc*ec)/(gc-1));
Pi_ch_rel(c2)=Pi_ch_off(c2)/Pi_ch;
Tt3_off(c2)=taw_ch_off(c2)*Tt2half_off(c2);

Pi_t_off=Pi_th_off*Pi_tl_off;

Pi_c_off(c2)=Pi_cl_off(c2)*Pi_ch_off(c2);
Pi_c_rel(c2)=Pi_c_off(c2)/Pi_c;

f_off(c2)=((Cpt*Tt4_off)-(Cpc*Tt3_off(c2)))/((hpr*eita_b)-(Cpt*Tt4_off));
f_rel(c2)=f_off(c2)/f;


MFP2_rel(c2)=(Pi_c_rel(c2)*sqrt((Pi_c^((gc-1)/gc))-1))/sqrt(((Pi_c_rel(c2)*Pi_c)^((gc-1)/gc))-1);

mo_rel(c2)=(MFP2_rel(c2)*Pt2_rel(c2))/sqrt(Tt2_rel(c2));
mo_off(c2)=mo_rel(c2)*mo;


Work_CompressorLow_off(c2)=Cpc*(Tt2half_off(c2)-Tt2_off(c2));
Work_CompressorHigh_off(c2)=Cpc*(Tt3_off(c2)-Tt2half_off(c2));


eqn1_off(c2)=vpa(Work_CompressorLow_off(c2)/Cceta2);
U1_off(c2)=sqrt(eqn*eqn1_off(c2));
N1_off(c2)=vpa((U1_off(c2)*60)/(2*pi*rm));
N1_rel(c2)=N1_off(c2)/N1;

eqn2_off(c2)=vpa(Work_CompressorHigh_off(c2)/Cceta2);
U2_off(c2)=sqrt(eqn*eqn2_off(c2));
N2_off(c2)=vpa((U2_off(c2)*60)/(2*pi*rm));
N2_rel(c2)=N2_off(c2)/N2;

Tt5_off(c2)=Tt4half_off*taw_tl_off;


Mth_off=1;
eqn3=(Mth_off/M9_offs)*(((1+((gt-1)*(M9_offs^2)/2))/(1+((gt-1)*(Mth_off^2)/2)))^((gt+1)/((2*gt)-2)));
M9_off=vpasolve(eqn3==E,M9_offs);

Tt9_off(c2)=Tt5_off(c2);
T9_off(c2)=Tt9_off(c2)/((1+((gt-1)*(M9_off^2)/2)));

Pt9_per_Pa_off(c2)=Pi_n*Pi_t_off*Pi_b*Pi_c_off(c2)*Pi_d*Pi_r_off(c2);

Pt9_off(c2)=Pt9_per_Pa_off(c2)*Po_off(c2);
P9_off(c2)=Pt9_off(c2)/((1+((gt-1)*(M9_off^2)/2))^(gt/(gt-1)));

rho9_off(c2)=P9_off(c2)/(R*T9_off(c2));
v9_off(c2)=sqrt(2*Cpt*Tt9_off(c2)*(1-((P9_off(c2)/Pt9_off(c2))^((gt-1)/gt))));

A9_off(c2)=A9;

F_off(c2)=(mo_off(c2)*(((1+f_off(c2))*v9_off(c2))-(vo_off(c2))))+(A9_off(c2)*(P9_off(c2)-Po_off(c2)));
F_per_mo_off(c2)=F_off(c2)/mo_off(c2);
Specific_FuelConsumption_off(c2)=vpa(f_off(c2)/F_per_mo_off(c2));
Thermal_Efficiency_off(c2)=vpa((((1+f_off(c2))*(v9_off(c2)^2))-((vo_off(c2)^2)))/(2*f_off(c2)*hpr));
Propulsive_Efficiency_off(c2)=vpa((2*F_per_mo_off(c2)*vo_off(c2))/(((1+f_off(c2))*(v9_off(c2)^2))-((vo_off(c2)^2))));
Overall_Efficiency_off(c2)=Thermal_Efficiency_off(c2)*Propulsive_Efficiency_off(c2);


F_rel(c2)=F_off(c2)/Freq;
Overall_Efficiency_rel(c2)=Overall_Efficiency_off(c2)/Overall_Efficiency;

Y(c2)=vpa(mo_rel(c2)/(F_rel(c2)*sqrt(To_rel(c2))));
X(c2)=vpa(F_rel(c2)/Po_rel(c2));
Z(c2)=vpa(N2_rel(c2)/sqrt(To_rel(c2)));
end

r=0;
figure(1)
for i=1:length(Mo_off)
r=r+1;
    hold on

plot3(MFP2_rel,Pi_c_rel,Mo_off)
xlabel('MFP2_r'), ylabel('πc_r')
caption=sprintf('H=%d Km',h_off(r)/1000);
    t=text(MFP2_rel(r),Pi_c_rel(r), caption);
    t.Color='blue';
    t.FontSize=6;
    
end


figure(2)
r=0;
for i=1:length(Mo_off)
r=r+1;
    hold on

plot3(F_per_mo_off,Specific_FuelConsumption_off,Mo_off)
xlabel('F/mo)off'), ylabel('S)off')
caption=sprintf('H=%d Km',h_off(r)/1000);
t=text(F_per_mo_off(r),Specific_FuelConsumption_off(r), caption);
t.Color='red';
t.FontSize=6;
end

figure(3)
r=0;
for i=1:length(Mo_off)
r=r+1;
    hold on

plot(Overall_Efficiency_off,Mo_off)
xlabel('OverallEfficiency)off'), ylabel('Mo)off')
caption=sprintf('H=%d Km  Mo=%0.2f',h_off(r)/1000,Mo_off(r));
t=text(Overall_Efficiency_off(r), Mo_off(r), caption);
t.Color='black';
t.FontSize=6;
end

