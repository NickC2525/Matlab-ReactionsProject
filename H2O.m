function Thermo = H2O(T)
H_std = -285.83; 
A = -203.61;
B = 1523.3;
C = -3196.4;
D = 2474.46;
E = 3.8553;
t = T/1000;
Cp_out = A + B*t + C*t^2 + D*t^3;
H_out = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 + H_std;

Thermo = [H_out, Cp_out];
end