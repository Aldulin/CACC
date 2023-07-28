% Ex12_1.m
% Define the linearized model parameters
u0=20; g=9.81; m=5000; f=0.015; Theta=0;
rho=1.202; A=3; Cd=0.5; uw=2;
% Calculate the equilibrium force, Fx0:
Fx0 = m*g*sin(Theta) + f*m*g + ...
0.5*rho*A*Cd*(u0+uw)^2;
% The time constant and dc gain are:
Tau=(m/(rho*A*Cd*(u0+uw)));
K=Tau/m;

sys = tf(K,[Tau 1]);

Ka = 10
TauA = 0.05
Ga = tf(Ka, [TauA 1 0])

Gm = Ga*sys

T0 = 16.27;
a = 9.16*10^(-3);
d = 1;
Kc = (4*d)/(pi*a);
%Kc=79.5;

%acordarea parametrilor PID
Md=2;
fid=45;
Kpm=Kc*Md*cos(fid);
Tim=(T0/pi)*((1+sin(fid))/cos(fid));
Tdm=0.25*Tim;
Kim=Kpm/Tim;
Kdm=Kpm*Tdm;


A = [0 -1 0 0; 0 -1/Tau 0 0; -1 0 0 0;0 0 1 0];
b = [0;K/Tau; 0 ; 0];
c = [1 0 0 0];
d = 0;
Ts = 0.01;
[Ad, bd, cd, dd] = c2dm(A,b,c,d,Ts,'zoh');

R = ctrb(Ad,bd);
inR = inv(R);
ht = inR(end,:);
sigma = 0.043;
tt = 20;
zeta = (-log(sigma))/(sqrt(pi^2+(log(sigma))^2));
wn = 4/(tt*zeta);
alpha1 = -2*exp(-zeta*wn*Ts)*cos(wn*Ts*sqrt(1-zeta^2));
alpha2 = exp(-2*zeta*wn*Ts);
alpha3 = exp(-4*zeta*wn*Ts);

Pcr = (Ad^2+alpha1*Ad+alpha2*eye(4))*(Ad-alpha3*eye(4))*(Ad-alpha3*eye(4));
ft = -ht*Pcr;



F1 = ft(1);
F2 = ft(2);
F3 = ft(3);
F4 = ft(4);

T = Ts;
