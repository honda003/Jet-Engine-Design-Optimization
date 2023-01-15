clc 
clear all
close all

c1=0;

%Flight Conditions
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

%Design Parameters 

Pi_c=4:1:20;
Tt4=900:50:1700;


%Intake 


[To, ao, Po, rho]=atmosisa(h);

R=287;

Cpc=(gc*R)/(gc-1);

Cpt=(R*gt)/(gt-1);

vo=Mo*sqrt(gc*R*To);

taw_r=(1+(((gc-1)/2)*(Mo^2)));
Tto=To*taw_r;

Pi_r=((1+(((gc-1)/2)*(Mo^2)))^(gc/(gc-1)));
Pto=Po*Pi_r;

%Diffuser

Tt2=Tto;
Pt2=Pto*Pi_d;

%Compressor
for k=1:length(Tt4)
    for i=1:length(Pi_c) 
        
        c1=c1+1;
    
taw_c(c1)=(Pi_c(i))^((gc-1)/(gc*ec));
   Tt3(c1)=Tt2*taw_c(c1);
   Pt3(c1)=Pi_c(i)*Pt2;
   
   Power_Compressor(c1)=Cpc*(Tt3(c1)-Tt2);


%Burner

  f(c1)=((Cpt*Tt4(k))-(Cpc*Tt3(c1)))/((hpr*eita_b)-(Cpt*Tt4(k)));
  Pt4(c1)=Pt3(c1)*Pi_b;
  
%Turbine
Power_Turbine(c1)=(Power_Compressor(c1))/eita_m;

Tt5(c1)=Tt4(k)-((Cpc*(Tt3(c1)-Tt2))/(eita_m*(1+f(c1))*Cpt));

taw_t(c1)=Tt5(c1)/Tt4(k);
Pi_t(c1)=(taw_t(c1))^(gt/((gt-1)*et));

Pt5(c1)=Pt4(c1)*Pi_t(c1);

%Nozzle

Tt9(c1)=Tt5(c1);
Pt9_per_Pa(c1)=Pi_n*(Pi_t(c1))*Pi_b*(Pi_c(i))*Pi_d*Pi_r;

if Pt9_per_Pa(c1) <= 1.85
    %%Convergent Nozzle with Full expansion (P9=Pa)
    
    P9(c1)=Po;
    Pt9(c1)=Pt9_per_Pa(c1)*P9(c1);
    
    M9(c1)=sqrt((2/(gt-1))*((Pt9_per_Pa(c1)^((gt-1)/gt))-1));
    
    T9(c1)=Tt9(c1)/(1+(((gt-1)*(M9(c1)^2))/2));
    
    v9(c1)=sqrt(2*Cpt*Tt9(c1)*(1-((1/Pt9_per_Pa(c1))^((gt-1)/gt))));
    
    Nozzle_Type(c1)=1;
else
    %%Convergent Divergent nozzle with Full expansion (P9=Pa)
    P9(c1)=Po;
    Pt9(c1)=Pt9_per_Pa(c1)*P9(c1);
    
    M9(c1)=sqrt((2/(gt-1))*((Pt9_per_Pa(c1)^((gt-1)/gt))-1));
    
    T9(c1)=Tt9(c1)/(1+(((gt-1)*(M9(c1)^2))/2));
    
    v9(c1)=sqrt(2*Cpt*Tt9(c1)*(1-((1/Pt9_per_Pa(c1))^((gt-1)/gt))));
    
    Nozzle_Type(c1)=2;
end

%Parameters

F_per_mo(c1)=vpa(((1+f(c1))*v9(c1))-(vo));
Specific_FuelConsumption(c1)=vpa(f(c1)/F_per_mo(c1));
Thermal_Efficiency(c1)=vpa((((1+f(c1))*(v9(c1)^2))-((vo^2)))/(2*f(c1)*hpr));
Propulsive_Efficiency(c1)=vpa((2*F_per_mo(c1)*vo)/(((1+f(c1))*(v9(c1)^2))-((vo^2))));
Overall_Efficiency(c1)=Thermal_Efficiency(c1)*Propulsive_Efficiency(c1);

        A(c1,1)=Pi_c(i);
        A(c1,2)=Tt4(k);
        A(c1,3)=F_per_mo(c1);
        A(c1,4)=Specific_FuelConsumption(c1);
        A(c1,5)=Thermal_Efficiency(c1);
        A(c1,6)=Propulsive_Efficiency(c1);
        A(c1,7)=Overall_Efficiency(c1);
        A(c1,8)=v9(c1);
        A(c1,9)=Nozzle_Type(c1);
        



end

end
 r=0;
 
 figure(1)
 xlabel('Specific Thrust'), ylabel('Specific Fuel Consumption'), zlabel('πc')

for z=1:length(Tt4):c1
    
   r=r+1;
    hold on
    plot3(A(z:z+length(Tt4)-1,3),A(z:z+length(Tt4)-1,4),A(z:z+length(Tt4)-1,1))
    caption=sprintf('Tt4=%d',Tt4(r));
    t=text(A(z,3),A(z,4), caption);
    t.Color='red';
    t.FontSize=6;
    
    
    [a1(r),b1(r)]=max(A(z:z+length(Tt4)-1,3));
    opt_pi_ST(r)=Pi_c(b1(r));
    
    [a2(r),b2(r)]=min(A(z:z+length(Tt4)-1,4));
    opt_pi_SF(r)=Pi_c(b2(r));
    

end

for c=1:length(Tt4)
    hold on
    plot3(A(c:length(Tt4):c1,3),A(c:length(Tt4):c1,4),A(c:length(Tt4):c1,1))
    
end

r=0;
figure(2)
for z=1:length(Tt4):c1
   
    r=r+1;
    hold on
    plot3(A(z:z+length(Tt4)-1,5),A(z:z+length(Tt4)-1,6),A(z:z+length(Tt4)-1,1))
    xlabel('Thermal Efficiency'), ylabel('Propulsive Efficiency'), zlabel('πc')
    
    [a3(r),b3(r)]=max(A(z:z+length(Tt4)-1,5));
    opt_pi_Therm(r)=Pi_c(b3(r));
    
    
    [a4(r),b4(r)]=max(A(z:z+length(Tt4)-1,6));
    opt_pi_Prop(r)=Pi_c(b4(r));
    
    [a5(r),b5(r)]=max(A(z:z+length(Tt4)-1,7));
    opt_pi_OA(r)=Pi_c(b5(r));
    

    
end
for c=1:length(Tt4)
    hold on
    plot3(A(c:length(Tt4):c1,5),A(c:length(Tt4):c1,6),A(c:length(Tt4):c1,1)) 
end

figure(3)
plot(a5,Tt4)
xlabel('Overall Efficiency'), ylabel('Turbine Inlet Temperature'), zlabel('πc')