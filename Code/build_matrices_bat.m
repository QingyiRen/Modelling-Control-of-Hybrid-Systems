function f=build_matrices_bat(bat,Ts)
A=1;
B1=-bat.eta_d*Ts;
B2=0;
B3=(bat.eta_d-bat.eta_c)*Ts;
B4=0;

E1=[0 0 0 0 0 0 1 -1]';
E2=[-1 1 0 0 -1 1 0 0]';
E3=[bat.u_min-eps bat.u_max -bat.u_max bat.u_min -bat.u_min bat.u_max 0 0]';
E4=[0 0 1 -1 1 -1 0 0]';
g5=[-eps bat.u_max 0 0 -bat.u_min bat.u_max bat.x_max 0]';

% Return as structure
f = struct('A',A,'B1',B1,'B2',B2,'B3',B3,'B4',B4,'E1',E1,'E2',E2,'E3',E3,'E4',E4,'g5',g5);
end