function Thermo = H2O2(T)
H_vap = 48.5; 
H_std = -136.1; 
A = 34.2567;
B = 55.184;
C = -35.154;
D = 9.087;
t = T/1000;
Cp_out = A + B*t + C*t^2 + D*t^3;
H_out = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - H_vap + H_std;

Thermo = [H_out, Cp_out];
end