function Thermo = Me(T)
Hf_std = -238.4; 
A = 2.211;
B = 12.22 * 10^-3;
C = -3.45 * 10^-6;
H_out = Hf_std + 0.00831*(A*T + B*T^2/2 + C*T^3/3);
Cp_out = 0.00831*(A + B*T + C*T^2);
Thermo = [H_out, Cp_out];
end