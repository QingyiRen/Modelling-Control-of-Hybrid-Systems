%Justine Kroese and Vivek Varma 
clc;
clear all;
close all;
%% 2.3
syms u a1 b1 a2 b2 a3 b3 a4 b4; 
u1=5;
u2=6.5;
u3=11;
%Non-linear functions
f1=(u^2)+4;
f2=4*u;
f3=(-9.44*u^3)+(166.06*u^2)-(948.22*u)+1790.28;
f4=-11.78*u+132.44;
f5=(4.01*u^2)-(20.94*4.01*u)+17.79+(4.01*10.47^2);
%Optimisation objective function
fun=@(x) int((f1-(x(1)+x(2)*u))^2,0, 2)+int((f2-(x(1)+x(2)*u))^2, 2, 5)+int((f3-(x(3)+x(4)*u))^2, 5, 6.5)+int((f3-(x(5)+x(6)*u))^2, 6.5, 7)+int((f4-(x(5)+x(6)*u))^2, 7, 9)+int((f5-(x(5)+x(6)*u))^2, 9, 11)+int((f5-(x(7)+x(8)*u))^2, 11, 15);
guess = [1 1 1 1 1 1 1 1];
options = optimoptions(@ga,'PopulationSize', 700, 'MaxGenerations', 1600,'EliteCount', 25, 'Display', 'iter'); 
[goptidiscrete, optvaldiscrete] =ga(fun, 8, [], [], [], [], [], [], [], [], options);
%% 2.3 Plots
pc=piecewise(0<=u & u<2, (u^2)+4, 2<=u & u<5, 4*u, 5<=u & u < 7, (-9.44*u^3)+(166.06*u^2)-(948.22*u)+1790.28, 7<=u & u< 9, -11.78*u+132.44, 9<=u & u<=15, (4.01*u^2)-(20.94*4.01*u)+17.79+(4.01*10.47^2));
figure
fplot(pc, [0 15])
hold on;
grid on;
pc2=piecewise(0<=u & u<5, goptidiscrete(1)+goptidiscrete(2)*u, 5<=u & u<6.5, goptidiscrete(3)+goptidiscrete(4)*u, 6.5<=u & u < 11, goptidiscrete(5)+goptidiscrete(6)*u, 11<=u & u<= 15, goptidiscrete(7)+goptidiscrete(8)*u);
fplot(pc2, [0 15])
xlabel('$u_d(k)$', 'Interpreter','latex')
ylabel('$f(u_d(k))/\hat{f}(u_d(k))$', 'Interpreter','latex')
legend('$f(u_d(k))$','$\hat{f}(u_d(k))$', 'Interpreter','latex')
title('Non-linear model vs PWA model of Diesel Generator')
%% 2.4
%Optimisation objective function
fun=@(x) int((pc-piecewise(0<=u & u<x(3), x(1)+x(2)*u, x(3)<=u & u<x(6), x(4)+x(5)*u, x(6)<=u & u<x(9), x(7)+x(8)*u, x(9)<=u & u<=15, x(10)+x(11)*u))^2, 0, 15)
%Initial guess=Solution of optisation in 2.3
x0=[1.6317452147803 3.5391964828094 5 -85.6286320070643 20.8671127088603 6.5 111.8700873269495 -9.1652452629377 11 -207.2789168622208 19.6967714691146];
options = optimoptions(@simulannealbnd,'Display', 'iter');
%Since the values of u1, u2, u3 have to be bounded
lb=[-Inf; -Inf; 0; -Inf; -Inf; 0; -Inf;-Inf;0;-Inf;-Inf];
ub=[Inf; Inf; 15; Inf; Inf; 15; Inf;Inf;15;Inf;Inf];
[gopti, optval]=simulannealbnd(fun,x0,lb,ub,options) 
%% 2.4 Plots
pc=piecewise(0<=u & u<2, (u^2)+4, 2<=u & u<5, 4*u, 5<=u & u < 7, (-9.44*u^3)+(166.06*u^2)-(948.22*u)+1790.28, 7<=u & u< 9, -11.78*u+132.44, 9<=u & u<=15, (4.01*u^2)-(20.94*4.01*u)+17.79+(4.01*10.47^2));
figure
fplot(pc, [0 15])
hold on;
grid on;
pc2=piecewise(0<=u & u<gopti(3), gopti(1)+gopti(2)*u, gopti(3)<=u & u<gopti(6), gopti(4)+gopti(5)*u, gopti(6)<=u & u < gopti(9), gopti(7)+gopti(8)*u, gopti(9)<=u & u<= 15, gopti(10)+gopti(11)*u);
fplot(pc2, [0 15])
xlabel('$u_d(k)$', 'Interpreter','latex')
ylabel('$f(u_d(k))/\hat{f}(u_d(k))$', 'Interpreter','latex')
legend('$f(u_d(k))$','$\hat{f}(u_d(k))$', 'Interpreter','latex')
title('Non-linear model vs PWA model of Diesel Generator with variable PWA limits')