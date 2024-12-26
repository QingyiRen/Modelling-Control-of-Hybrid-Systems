%% Assignment Hybrid Systems
% Vivek Varma and Justine Kroese
clear all; close all; clc
%% Exercise 2.7
%% Variables
Ts=0.25;                % sampling time [h]

% Specification batteries
bat1.eta_c=0.9;         % charging efficiency battery 1
bat2.eta_c=0.95;        % charging efficiency battery 2
bat1.eta_d=0.8;         % discharging efficiency battery 1
bat2.eta_d=0.77;        % dicharging efficiency batttery 2
bat1.x_max=48;          % max stored energy battery 1 [kWh]
bat2.x_max=64;          % max stored energy battery 2 [kWh]
bat1.u_min=-3;
bat2.u_min=-4;
bat1.u_max=2;
bat2.u_max=3;

% Specification diesel generator
dies.x_min=20;         % [kg]
dies.x_max=120;        % [kg]
dies.u_max=15;         % [kW]
dies.Rf=0.4;           % [kg/h]

% ....
a1=1;
a2=1;
a3=1;
a4=1;
b1=1;
b2=1;
b3=1;
b4=1;
dies.a=[a1 a2 a3 a4];
dies.b=[b1 b2 b3 b4];
delta1=0;
delta2=0;
delta3=0;
delta4=0;
dies.delta_d=[delta1; delta2; delta3; delta4];
u1=3.75;
u2=6.5;
u3=10.25;

%% Description of the MLD matrices

% Build MLD matrices batteries
MLDbat1=build_matrices_bat(bat1,Ts);
MLDbat2=build_matrices_bat(bat2,Ts); 

% Build MLD matrices diesel generator
MLDdies.A=1;
MLDdies.B1=0;
MLDdies.B2=-Ts*dies.a;
MLDdies.B3=-Ts*dies.b;
MLDdies.B4=dies.Rf*Ts;
MLDdies.E1=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1]';
MLDdies.E2=[-1 1 -1 1 -1 1 -1 1 1 0 0 -1 1 0 0 -1 1 0 0 -1 1 0 0 -1 1 0 0]';
MLDdies.E3=[0 0 0 0; -u1+dies.u_max+eps 0 0 0;0 u1 0 0; 0 -u2+dies.u_max+eps 0 0;0 0 u2 0; 0 0 -u3+dies.u_max+eps 0;0 0 0 u3; 0 0 0 0;1 1 1 1; -u1 0 0 0;0 0 0 0; 0 0 0 0;dies.u_max 0 0 0; 0 -u2 0 0; 0 u1 0 0;0 0 0 0; 0 dies.u_max 0 0;0 0 -u3 0; 0 0 u2 0; 0 0 0 0;0 0 dies.u_max 0; 0 0 0 -15;0 0 0 u3;0 0 0 0; 0 0 0 dies.u_max;0 0 0 0;0 0 0 0 ];
MLDdies.E4=[0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 0; 1 0 0 0;-1 0 0 0; 1 0 0 0;-1 0 0 0; 0 1 0 0; 0 -1 0 0;0 1 0 0; 0 -1 0 0;0 0 1 0; 0 0 -1 0; 0 0 1 0;0 0 -1 0; 0 0 0 1;0 0 0 -1;0 0 0 1; 0 0 0 -1;0 0 0 0;0 0 0 0 ];
MLDdies.g5=[0 dies.u_max 0 dies.u_max 0 dies.u_max 0  15 1 0 0 0 dies.u_max 0 0 0 dies.u_max 0 1 0 dies.u_max 0 0 0 dies.u_max dies.x_max dies.x_max]';

% Combining the matrices
sys.A=blkdiag(MLDbat1.A,MLDbat2.A,MLDdies.A);
sys.B1=blkdiag(MLDbat1.B1,MLDbat2.B1,MLDdies.B1);
sys.B2=blkdiag(MLDbat1.B2,MLDbat2.B2,MLDdies.B2);
sys.B3=blkdiag(MLDbat1.B3,MLDbat2.B3,MLDdies.B3);
sys.B4=[MLDbat1.B4;MLDbat2.B4;MLDdies.B4];

sys.E1=blkdiag(MLDbat1.E1,MLDbat2.E1,MLDdies.E1);
sys.E2=blkdiag(MLDbat1.E2,MLDbat2.E2,MLDdies.E2);
sys.E3=blkdiag(MLDbat1.E3,MLDbat2.E3,MLDdies.E3);
sys.E4=blkdiag(MLDbat1.E4,MLDbat2.E4,MLDdies.E4);
sys.g5=[MLDbat1.g5;MLDbat2.g5;MLDdies.g5];


%% Exercise 2.7
%% Variables
Nb=2;                   % amount of batteries
Np=5;                   % prediction horizon
Nc=3;                  % control horizon

% weight factors
weight.Wb1=3;              % weight factor batteries
weight.Wb2=4;
weight.Wd=10;                  % weight factor diesel generator
weight.W_fuel=4;               % weight factor fuel
weight.We=0.4;                % weight factor ?

%weight.Ce=1;                   % benefit for importing energy
%weight.P_load=1;               % power of load?

M1=[];                             %M1
for i=1:Np
    M1=[M1;sys.A^i];
end

M3=[];                              %M3
sum=0;
for i=1:Np
    sum=sum+sys.A^(i-1)*sys.B4;
    M3=[M3; sum];
end

T1=[];
T2=[];
T3=[];
matrt1=[];
matrt2=[];
matrt3=[];
for i=1:Np
    matrt1=[sys.A^(Np-1)*sys.B1 matrt1];
    T1=[T1;matrt1 zeros(size(sys.B1, 1), size(sys.B1, 2)*(Np-i))];       %T1
    matrt2=[sys.A^(Np-1)*sys.B2 matrt2];
    T2=[T2;matrt2 zeros(size(sys.B2, 1), size(sys.B2, 2)*(Np-i))];      %T2
    matrt3=[sys.A^(Np-1)*sys.B3 matrt3];
    T3=[T3;matrt3 zeros(size(sys.B3, 1), size(sys.B3, 2)*(Np-i))];      %T3
end

Inu=eye(3);
Ku=[];                                                          %Ku
for i=1:Nc
    Ku=blkdiag(Ku, Inu);
end
for j=1:(Np-Nc)
    Ku=[Ku; zeros(3, 3*(Nc-1)) Inu];
end

Indelta=eye(6);
Kdelta=[];                                                  %Kdelta
for i=1:Np
    Kdelta=blkdiag(Kdelta, Indelta);
end

Inz=eye(6);
Kz=[];                                                      %Kz
for i=1:Np
    Kz=blkdiag(Kz, Inz);
end

M2=[T1*Ku T2*Kdelta T3*Kz];                                     %M2

E1hat=[];
E2hat=[];
E3hat=[];
E4hat=[];
for i=1:Np
    E1hat=blkdiag(E1hat, sys.E1);
    E2hat=blkdiag(E2hat, sys.E2);
    E3hat=blkdiag(E3hat, sys.E3);
    E4hat=blkdiag(E4hat, sys.E4);
end

E1hat=blkdiag(E1hat, eye(3));
E1hat=[E1hat; zeros(3, 3*Np) -eye(3)];                      %E1hat

E2hat=[E2hat; zeros(6, size(E2hat, 2))];                    %E2hat
E3hat=[E3hat; zeros(6, size(E3hat, 2))];                    %E3hat
E4hat=[E4hat; zeros(6, size(E4hat, 2))];                    %E4hat

g5hat=[];
for i=1:Np
    g5hat=[g5hat; sys.g5];
end
xupper=[48; 64; 120];                               %Given in the Question
xlower=[0 ; 0;20];

g5hat=[g5hat; xupper; -xlower];                             %g5hat

F1=E1hat*[zeros(3, size(M2, 2)); M2]+[E2hat*Ku E3hat*Kdelta E4hat*Kz]; %F1
F2=g5hat-E1hat*[zeros(3,1); M3];                                %F2
F3=-E1hat*[eye(3); M1];                                         %F3