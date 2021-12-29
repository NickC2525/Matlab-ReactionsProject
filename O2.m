function Thermo = O2(T)
A = 31.322;
B = -20.235;
C = 57.866;
D = -36.506;
t = T/1000;
Cp_out = A + B*t + C*t^2 + D*t^3;
H_out = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4;
Thermo = [H_out, Cp_out];
end