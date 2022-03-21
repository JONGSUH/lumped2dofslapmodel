%% z0_2DOF_SlapMain
clear all; clc;  close all
%%
global M dt 
%% 2DOF
m1=1; m2=.5; 
M = diag([m1 m2]); 
% k1=1; k2=0.5*k1; k3=1.5*k1; % dU=~0
k1=1; k2=0.5*k1; k3=0.5*k1; % dU=0
K=[k1+k2 -k2; -k2 k2+k3];
C = zeros(size(M));
fs = 100; dt = 1/fs;
T=1000; tt = 0:dt:T-1/fs;

Initial_gap = 0.005; 

[ve,va] = eig(inv(M)*K); ve1 = ve(:,1); ve2=ve(:,2); 
ref_ve = ve2;

fn =sqrt(diag(va))/(2*pi)'
v0 = 0.00; d0 = 0.0500;
v_n0 = -v0*ref_ve; d_n0 = -d0*ref_ve;

fre0=0.01; fre1=20; k=(fre1/fre0)^(1/T); phi_0=0; 
%%
% AnalType = 'Penalty'; 
AnalType ='AugLag';
F=zeros(2,length(tt));
A=0.0000; % Free vibration
F(2,:) = -A*sin(phi_0 + 2*pi*fre0*(k.^(tt) - 1)/log(k));
[displ,velo] = z1_LumpedModelSlap(M,C,K,Initial_gap,v_n0,d_n0,F,tt,AnalType);