function dW = diff_eqH2O2(W, C)
global R T_rxnr vo F_Me 

% Direct synthesis reaction.
E_H2O2 = 23840; 
A_H2O2 = 4.18e12 / (100)^6; 
k_H2O2= A_H2O2 * exp(-E_H2O2/(R*T_rxnr));

% Water formation reaction.
E_H2O = 45310;
A_H2O = 4.27e14 / 100^(9/2);
k_H2O = A_H2O * exp(-E_H2O/(R*T_rxnr));

% Decomposition reaction.
E_Decomp = 45060;
A_Decomp = 2.61e10 / 100^3;
kd = A_Decomp * exp(-E_Decomp/(R*T_rxnr));

% Hydrogenation reaction.
E_Hyd = 46580;
A_Hyd = 5.65e17 / 100^6;
k_Hyd = A_Hyd * exp(-E_Hyd/(R*T_rxnr));

% Reaction rates.
R_H2O2 = k_H2O2 * C(4) * C(3)^2;
R_H2O = k_H2O * C(3) * C(4)^(1/2);
R_Decomp = kd * C(1);
R_Hyd = k_Hyd * C(1) * C(3);

% Differential concentration change.
dC_H2O2 = (R_H2O2 - R_Decomp - R_Hyd)/vo;
dC_H2O = (R_H2O + R_Decomp + R_Hyd)/vo;
dC_H2 = (-R_Hyd - R_H2O2 - R_H2O)/vo;
dC_O2 = (-R_H2O2 - 0.5*R_H2O + 0.5*R_Decomp)/vo;

% Get enthalpy & Cp at current T.
td_h2o = H2O(T_rxnr);
td_h2o2 = H2O2(T_rxnr);
td_h2 = H2(T_rxnr);
td_o2 = O2(T_rxnr);
td_me = Me(T_rxnr);

% Enthalpy(1).
H_H2O = td_h2o(1);
H_H2O2 = td_h2o2(1);
H_H2 = td_h2(1);
H_O2 = td_o2(1);

% Cp(2)
Cp_H2O = td_h2o(2);
Cp_H2O2 = td_h2o2(2);
Cp_H2 = td_h2(2);
Cp_O2 = td_o2(2);
Cp_Me = td_me(2);

% Heats of reaction for each reaction.
Hr_H2O2 = -H_O2 - H_H2 + H_H2O2;

Hr_H2O = -H_H2 - 0.5*H_O2 + H_H2O;

Hr_Decomp = -H_H2O2 + H_H2O + 0.5*H_O2;

Hr_Hyd = -H_H2O2 - H_H2 + 2*H_H2O;

% Total heat generated = sum(-ri*dHrx,i)
Hgen = -1*(R_H2O2*Hr_H2O2 + R_H2O*Hr_H2O + R_Decomp*Hr_Decomp + R_Hyd*Hr_Hyd);
FCp = vo*(C(3)*Cp_H2 + C(4)*Cp_O2 + C(2)*Cp_H2O + ...
C(1)*Cp_H2O2) + F_Me*Cp_Me;

% Differential rate of temperature change. Assume adiabatic.
dT = Hgen/FCp;
% Pack & return results.
dW = [dC_H2O2; dC_H2O; dC_H2; dC_O2; dT];
end
