function [Gamma] = BurnCalc(C,X,Time,TS)

T = C(X.MN,:);
x = find(T >= 317.15);
tburn = TS(x(1));

t = tburn:Time.dt:50;
T = T(x(1):end);
Integrand = (2*10^98)*exp(-12017/(T'-273.15));

Gamma = Time.dt * trapz(t,Integrand);

end
