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


a1=1.6318;
a2=-85.6286;
a3=111.8701 ;
a4=-207.2789;
b1=3.5292;
b2=20.8671;
b3=-9.1653;
b4=19.6968;
dies.a=[a1 a2 a3 a4];
dies.b=[b1 b2 b3 b4];
% delta1=0;
% delta2=0;
% delta3=0;
% delta4=0;
% dies.delta_d=[delta1; delta2; delta3; delta4];
u1=5;
u2=6.5;
u3=11;

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
MLDdies.E1=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1]';
MLDdies.E2=[-1 1 -1 1 -1 1 -1 1 1 0 0 0 -1 1 0 0 -1 1 0 0 -1 1 0 0 -1 1 0 0]';
MLDdies.E3=[0 0 0 0; -u1+dies.u_max+eps 0 0 0;0 u1 0 0; 0 -u2+dies.u_max+eps 0 0;0 0 u2 0; 0 0 -u3+dies.u_max+eps 0;0 0 0 u3; 0 0 0 0;-dies.u_max -dies.u_max -dies.u_max -dies.u_max;1 1 1 1; -u1 0 0 0;0 0 0 0; 0 0 0 0;dies.u_max 0 0 0; 0 -u2 0 0; 0 u1 0 0;0 0 0 0; 0 dies.u_max 0 0;0 0 -u3 0; 0 0 u2 0; 0 0 0 0;0 0 dies.u_max 0; 0 0 0 -15;0 0 0 u3;0 0 0 0; 0 0 0 dies.u_max;0 0 0 0;0 0 0 0 ];
MLDdies.E4=[0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 0;0 0 0 0; 1 0 0 0;-1 0 0 0; 1 0 0 0;-1 0 0 0; 0 1 0 0; 0 -1 0 0;0 1 0 0; 0 -1 0 0;0 0 1 0; 0 0 -1 0; 0 0 1 0;0 0 -1 0; 0 0 0 1;0 0 0 -1;0 0 0 1; 0 0 0 -1;0 0 0 0;0 0 0 0 ];
MLDdies.g5=[0 dies.u_max 0 dies.u_max 0 dies.u_max 0  15 0 1 0 0 0 dies.u_max 0 0 0 dies.u_max 0 1 0 dies.u_max 0 0 0 dies.u_max dies.x_max -dies.x_min]';

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


%% Exercise 2.8
%% Variables
Nb=2;                   % amount of batteries
Np=25;                   % prediction horizon
Nc=25;                  % control horizon

Nx=size(sys.A,2);
Nu=size(sys.B1,2);
Ndelta=size(sys.B2,2);
Nz=size(sys.B3,2);

Nrho=(Nb+1)*Np;

x0=[10;10;50];
xk=x0;
delta_min1=[0;0;0;0;0;0];
deltak_min1=delta_min1;

% weight factors
weight.Wb1=3;              % weight factor batteries
weight.Wb2=4;
weight.Wd=10;                  % weight factor diesel generator
weight.W_fuel=4;               % weight factor fuel
weight.We=0.4;                % weight factor ?

k=0:1:Np-1;
weight.Ce=50+50*sin(pi*Ts*k/12)';            % benefit for importing energy
weight.P_load=zeros(Np,1);
for i=1:Np-1
    if i>20 && i<51
        weight.P_load(i)=30+2*i;               % power of load?
    elseif i>=51
        weight.P_load(i)=45;
    end
end

MILPsys=MILP_model(sys,bat1,bat2,dies,Nb,Np,Nc,weight,xk,deltak_min1);

% solve the problem for the initial conditions
model.obj=MILPsys.c;
model.A=MILPsys.A;
model.rhs=MILPsys.b;
model.modelsense ='min';
model.vtype=repmat('C',size(MILPsys.c,2),1);
model.vtype(MILPsys.bin_part)='B';
model.lb =-inf(size(MILPsys.c,2),1);

optimum=gurobi(model);
opt.V=optimum.x;
opt.rho=opt.V(1:Nrho);
opt.u=opt.V(Nrho+1:(Nrho+Nu*Nc));
opt.delta=opt.V((Nrho+Nu*Nc+1):(Nrho+(Nu+Ndelta)*Nc));
opt.z=opt.V((Nrho+(Nu+Ndelta)*Nc+1):end);


%% Exercise 2.9
%%
%Simulation time
N=48/Ts;                

% initialisation
x29=zeros(3,Np+1);
u29=zeros(3,Np+1);
delta29=zeros(6,Np+1);
z29=zeros(6,Np+1);
x29(:,1)=[10;10;50];
constraints29=zeros(1,Np+1);
J29=zeros(1,Np);

for j=1:N
    t=j+k;                                  % for each step in the horizon we follow the whole prediction horizon
    weight.Ce=50+50*sin(pi*Ts*t/12)';           
    weight.P_load=zeros(Np,1);
   
    weight.P_load(t<=20)=0;
    weight.P_load(t>=21&t<=50)=30+2*t(t>= 21&t<=50);
    weight.P_load(t>=51)=45;
    
    if j>1
        [MILPsys,Cx,M2,M3]=MILP_model(sys,bat1,bat2,dies,Nb,Np,Nc,weight,x29(:,j),delta29(:,j-1));
    else
        [MILPsys29,Cx,M2,M3]=MILP_model(sys,bat1,bat2,dies,Nb,Np,Nc,weight,x29(:,j),delta_min1);
    end
    
    model=struct();
    model.obj=MILPsys29.c;
    model.A=MILPsys29.A;
    model.rhs=MILPsys29.b;
    model.modelsense ='min';
    model.vtype=repmat('C',size(MILPsys29.c,2),1);
    model.vtype(MILPsys29.bin_part)='B';
    model.lb =-inf(size(MILPsys29.c,2),1);                    %otherwise no negative values
    
    params.outputflag = 0;
    optimum29=gurobi(model,params);
    opt29.V=optimum29.x;
    opt29.rho=opt29.V(1:Nrho);
    opt29.u=opt29.V(Nrho+1:(Nrho+Nu*Nc));
    opt29.delta=opt29.V((Nrho+Nu*Nc+1):(Nrho+(Nu+Ndelta)*Nc));
    opt29.z=opt29.V((Nrho+(Nu+Ndelta)*Nc+1):end);

    u29(:,j)=opt29.u(1:Nu);
    delta29(:,j)=opt29.delta(1:Ndelta);
    z29(:,j)=opt29.z(1:Nz);

    x29(:,j+1) = sys.A*x29(:,j)+sys.B1*u29(:,j)+sys.B2*delta29(:,j)+sys.B3*z29(:,j)+sys.B4;
    % check constraints:
    constraints29(j)=all(sys.E1*x29(:,j)+sys.E2*u29(:,j)+sys.E3*delta29(:,j)+sys.E4*z29(:,j)-sys.g5<=eps);
    
    J29(j)=optimum29.objval+Cx*M2*x29(:,j)+Cx*M3+[weight.We weight.We weight.W_fuel]*x29(:,j)+weight.Ce'*weight.P_load;

end

figure
plot(J29)
grid on;
title('Optimal Operation Cost of the Microgrid for 48 Hours')
xlabel('Time step')
ylabel('Operational Cost')

figure
subplot(3,1,1)
plot(x29(1,:))
hold on;
plot([1:N], ones(1, N)*bat1.x_max)
hold on;
plot([1:N], ones(1, N)*0)
grid on;
title('Stored Energy in Battery 1 for 48 Hours')
xlabel('Time step')
ylabel('Stored Energy(kWh)')
legend('x_{b,1}', 'x_{b,1} max', 'x_{b,1} min')
subplot(3,1,2)
plot(x29(2,:))
hold on;
plot([1:N], ones(1, N)*bat2.x_max)
hold on;
plot([1:N], ones(1, N)*0)
grid on;
title('Stored Energy in Battery 2 for 48 Hours')
xlabel('Time step')
ylabel('Stored Energy(kWh)')
legend('x_{b,2}', 'x_{b,2} max', 'x_{b,2} min')
subplot(3,1,3)
plot(x29(3,:))
hold on;
plot([1:N], ones(1, N)*dies.x_max)
hold on;
plot([1:N], ones(1, N)*dies.x_min)
grid on;
title('Remaning Fuel in Diesel Generator for 48 Hours')
xlabel('Time step')
ylabel('Remaining Fuel(kg)')
legend('x_{d}', 'x_{d} max', 'x_{d} min')

figure
subplot(3,1,1)
plot(u29(1,:))
hold on;
plot([1:N], ones(1,N)*bat1.u_max)
hold on;
plot([1:N], ones(1, N)*bat1.u_min)
grid on;
title('Exchanged Power in Battery 1 for 48 Hours')
xlabel('Time step')
ylabel('Exchanged Power(kW)')
legend('u_{b,1}', 'u_{b,1} max', 'u_{b,1} min')
subplot(3,1,2)
plot(u29(2,:))
hold on;
plot([1:N], ones(1, N)*bat2.u_max)
hold on;
plot([1:N], ones(1, N)*bat2.u_min)
grid on;
title('Exchanged Power in Battery 2 for 48 Hours')
xlabel('Time step')
ylabel('Exchanged Power(kW)')
legend('u_{b,2}', 'u_{b,2} max', 'u_{b,2} min')
subplot(3,1,3)
plot(u29(3,:))
hold on;
plot([1:N], ones(1, N)*dies.u_max)
hold on;
plot([1:N], ones(1, N)*0)
grid on;
title('Generated Power in Diesel Generator for 48 Hours')
xlabel('Time step')
ylabel('Generated Power(kW)')
legend('u_{d}', 'u_{d} max', 'u_{d} min')

figure
subplot(3,1,1)
plot(delta29(1,:))
grid on;
title('Switching Signal in Battery 1 for 48 Hours')
xlabel('Time step')
ylabel('Switching Signal')
subplot(3,1,2)
plot(delta29(2,:))
grid on;
title('Switching Signal in Battery 2 for 48 Hours')
xlabel('Time step')
ylabel('Switching Signal')
subplot(3,1,3)
plot(delta29(3:end,:)')
grid on;
title('Switching Signal in Diesel Generator for 48 Hours')
xlabel('Time step')
ylabel('Switching Signal')
legend('Delta1', 'Delta2' , 'Delta3', 'Delta4')




%% 2.10
%%
%simmulation time
N=48/Ts;                

% initialisation
x210=zeros(3,Np+1);
u210=zeros(3,Np+1);
delta210=zeros(6,Np+1);
z210=zeros(6,Np+1);
x210(:,1)=[10;10;50];
constraints210=zeros(1,Np+1);
J210=zeros(1,Np);
T12h=0;

for j=1:N
    t=j+k;                                  % for each step in the horizon we follow the whole prediction horizon
    weight.Ce=50+50*sin(pi*Ts*t/12)';           
    weight.P_load=zeros(Np,1);
   
    weight.P_load(t<=20)=0;
    weight.P_load(t>=21&t<=50)=30+2*t(t>= 21&t<=50);
    weight.P_load(t>=51)=45;
    
    T12h=find(mod(t,48)==0);%[];
    
    
    if j>1
        [MILPsys210,Cx,M2,M3]=MILP_model210(sys,bat1,bat2,dies,Nb,Np,Nc,weight,x210(:,j),delta210(:,j-1),T12h);
    else 
        [MILPsys210,Cx,M2,M3]=MILP_model210(sys,bat1,bat2,dies,Nb,Np,Nc,weight,x210(:,j),delta_min1,T12h);
    end
    
    model=struct();
    model.obj=MILPsys210.c;
    model.A=MILPsys210.A;
    model.rhs=MILPsys210.b;
    model.modelsense ='min';
    model.vtype=repmat('C',size(MILPsys210.c,2),1);
    model.vtype(MILPsys210.bin_part)='B';
    model.lb =-inf(size(MILPsys210.c,2),1);                    %otherwise no negative values
    
    params.outputflag = 0;
    optimum210=gurobi(model,params);
    opt210.V=optimum210.x;
    opt210.rho=opt210.V(1:Nrho);
    opt210.u=opt210.V(Nrho+1:(Nrho+Nu*Nc));
    opt210.delta=opt210.V((Nrho+Nu*Nc+1):(Nrho+(Nu+Ndelta)*Nc));
    opt210.z=opt210.V((Nrho+(Nu+Ndelta)*Nc+1):end);

    u210(:,j)=opt210.u(1:Nu);
    delta210(:,j)=opt210.delta(1:Ndelta);
    z210(:,j)=opt210.z(1:Nz);

    x210(:,j+1) = sys.A*x210(:,j)+sys.B1*u210(:,j)+sys.B2*delta210(:,j)+sys.B3*z210(:,j)+sys.B4;
    % check constraints:
    constraints210(j)=all(sys.E1*x210(:,j)+sys.E2*u210(:,j)+sys.E3*delta210(:,j)+sys.E4*z210(:,j)-sys.g5<=eps);
    
    J210(j)=optimum210.objval+Cx*M2*x210(:,j)+Cx*M3+[weight.We weight.We weight.W_fuel]*x210(:,j)+weight.Ce'*weight.P_load;

end



figure
plot(J210)
grid on;
title('Optimal Operation Cost of the Microgrid with Additional Constraint for 48 Hours')
xlabel('Time step')
ylabel('Operational Cost')
 
figure
subplot(3,1,1)
plot(x210(1,:))
hold on;
plot([1:N], ones(1, N)*bat1.x_max)
hold on;
plot([1:N], ones(1, N)*0)
grid on;
title('Stored Energy in Battery 1 with Additional Constraint for 48 Hours')
xlabel('Time step')
ylabel('Stored Energy(kWh)')
legend('x_{b,1}', 'x_{b,1} max', 'x_{b,1} min')
subplot(3,1,2)
plot(x210(2,:))
hold on;
plot([1:N], ones(1, N)*bat2.x_max)
hold on;
plot([1:N], ones(1, N)*0)
grid on;
title('Stored Energy in Battery 2 with Additional Constraint for 48 Hours')
xlabel('Time step')
ylabel('Stored Energy(kWh)')
legend('x_{b,2}', 'x_{b,2} max', 'x_{b,2} min')

subplot(3,1,3)
plot(x210(3,:))
hold on;
plot([1:N], ones(1, N)*dies.x_max)
hold on;
plot([1:N], ones(1, N)*dies.x_min)
grid on;
title('Remaning Fuel in Diesel Generator with Additional Constraint for 48 Hours')
xlabel('Time step')
ylabel('Remaining Fuel(kg)')
legend('x_{d}', 'x_{d} max', 'x_{d} min')

figure
subplot(3,1,1)
plot(u210(1,:))
hold on;
plot([1:N], ones(1, N)*bat1.u_max)
hold on;
plot([1:N], ones(1, N)*bat1.u_min)
grid on;
title('Exchanged Power in Battery 1 with Additional Constraint for 48 Hours')
xlabel('Time step')
ylabel('Exchanged Power(kW)')
legend('u_{b,1}', 'u_{b,1} max', 'u_{b,1} min')
subplot(3,1,2)
plot(u210(2,:))
hold on;
plot([1:N], ones(1, N)*bat2.u_max)
hold on;
plot([1:N], ones(1, N)*bat2.u_min)
grid on;
title('Exchanged Power in Battery 2 with Additional Constraint for 48 Hours')
xlabel('Time step')
ylabel('Exchanged Power(kW)')
legend('u_{b,2}', 'u_{b,1} max', 'u_{b,2} min')
subplot(3,1,3)
plot(u210(3,:))
hold on;
plot([1:N], ones(1, N)*dies.u_max)
hold on;
plot([1:N], ones(1, N)*0)
grid on;
title('Generated Power in Diesel Generator with Additional Constraint for 48 Hours')
xlabel('Time step')
ylabel('Generated Power(kW)')
legend('u_{d}', 'u_{d} max', 'u_{d} min')


figure
subplot(3,1,1)
plot(delta210(1,:))
grid on;
title('Switching Signal in Battery 1 for 48 Hours')
xlabel('Time step')
ylabel('Switching Signal')
subplot(3,1,2)
plot(delta210(2,:))
grid on;
title('Switching Signal in Battery 2 for 48 Hours')
xlabel('Time step')
ylabel('Switching Signal')
subplot(3,1,3)
plot(delta210(3:end,:)')
grid on;
title('Switching Signal in Diesel Generator for 48 Hours')
xlabel('Time step')
ylabel('Switching Signal')
legend('Delta1', 'Delta2' , 'Delta3', 'Delta4')

