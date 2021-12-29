clear;close all;clc
global R T_rxnr F_Me vo 

%For the purpose of the report 1 m^3/s was chosen 
prompt = 'Enter volumetric feed rate:'; 
vo = input(prompt);


 % m3 per s, basis value
R = 8.314; 
T_rxnr = 270; 


rho_Me = 813.3; %At T_rxnr
MM_Me = 32.04/1000; %kg/mol
v_Me = rho_Me / MM_Me; % mol/m^3

% Total flow rate of Methanol.
F_Me = vo*v_Me;


% Assume ideal gas (Pa)

P_rxnr=ones(1,30);
Max_H2O2_out=ones(1,30);

n=1;
for P_rxnrbar=1:30
    
P_rxnr(n) = P_rxnrbar * 100000; 
P_CO2 = 0.76 * P_rxnr(n); 
P_O2 = 0.21 * P_rxnr(n); 
P_H2 = 0.03 * P_rxnr(n); 

% Convert to concentration
C_H2Oo = 0;
C_H2O2o = 0;
C_H2o = P_H2/(R*T_rxnr); 
C_O2o = P_O2/(R*T_rxnr); 

% Enthalpy of feed (kJ/s)
H_H2O_in = 0; 
H_H2O2_in = 0; 
H_H2_in = H2(T_rxnr); 
H_O2_in = O2(T_rxnr); 
H_Me_in = Me(T_rxnr);
% Value (1) in td_ array is enthalpy, (2) is Cp.
H_in = vo*(C_H2o*H_H2_in(1) + C_O2o*H_O2_in(1) + C_H2Oo*H_H2O_in + ...
C_H2O2o*H_H2O2_in) + F_Me*H_Me_in(1); % kJ/s

% Set conditions and solve the ODE.
% H2O2, H2O, H2, O2
Co = [0 0 C_H2o C_O2o T_rxnr];
[W C] = ode23s(@diff_eqH2O2, [0 10], Co);
Max_H2O2_out(n) = max(C(:,1)); % find the most H2O2 production possible
k = find(C(:,1)==Max_H2O2_out(n)); % get its array index

C_H2O2_out = C(k,1);
C_H2O_out = C(k,2);
C_H2_out = C(k,3);
C_O2_out = C(k,4);
T_out = C(k,5);
C(k,:);

cat_amt = W(k); % mol
Production= Max_H2O2_out(n)*34*vo*(1/1000)*3600*8400/1000; %Prod tons per year

% Unpacking from function
H_H2O2_out = H2O2(T_out);
H_H2O_out = H2O(T_out);
H_H2_out = H2(T_out);
H_O2_out = O2(T_out);
H_Me_out = Me(T_out);

% Enthalpy leaving
H_H2O2_out = C_H2O2_out*(vo)*H_H2O2_out(1); 
H_H2O_out = C_H2O_out*(vo)*H_H2O_out(1); 
H_H2_out = C_H2_out*(vo)*H_H2_out(1); 
H_O2_out = C_O2_out*(vo)*H_O2_out(1); 
H_Me_out = F_Me*H_Me_out(1);
H_out = H_H2_out + H_O2_out + H_H2O_out + H_H2O2_out + H_Me_out;

Q = H_out - H_in; %Heat duty (kW)


PdMW=106.42; %Molecular weight of Pd 
SiMW=28.0855; %Molecular weight of Si 
OMW=16; %Molecular weight of O 

SiO2MW=SiMW+2*OMW; %Molecular weight of SiO2 
fracPd=0.04/100; %weight fraction Pd 
CatMW=fracPd*PdMW+(1-fracPd)*SiO2MW; %Avg molecular weight of catalyst 

Wkg=W*CatMW./1000/.0004; 
F=vo*C; 

n=n+1;
end

figure 
plot(Wkg, F(:,1), Wkg, F(:,2), Wkg, F(:,3)); 
title('Flow Rates vs Catalyst Weight'); 
xlabel('Catalyst Weight (kgcat)'); 
ylabel('Flow Rate (mol/s)'); 
legend('F_H_2_O_2', 'F_H_2_O', 'F_H_2'); 
S=F(:,1)./F(:,2); 

figure 
plot(Wkg,S,'r'); 
title('Selectivity vs Catalyst Weight'); 
xlabel('Catalyst Weight (kgcat)'); 
ylabel('Selectivity (S_H_2_O_2_/_H_2_O)'); 


P=P_rxnr/10000;

figure
plot(P,Max_H2O2_out)
title('H_2O_2 Concentration vs Pressure')
xlabel('Reactor Pressure [bar]')
ylabel('Concentration of H_2O_2 [mol/m^3]')


