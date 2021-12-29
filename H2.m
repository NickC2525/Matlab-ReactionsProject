function Thermo = H2(T)
A = 33.07;
B = -11.36;
C = 11.43;
D = -2.77;
t = T/1000;
Cp_out = A + B*t + C*t^2 + D*t^3;
H_out = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4;
Thermo = [H_out, Cp_out];
end