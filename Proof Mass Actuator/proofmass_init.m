%Shishir Khanal
%MCE 6642
%Feedback control of proof mass actuator system
%--------------------------------------
clc;clear; close all;
% This script generates the linearized A,B,C,D 
% matrices for CE 5
%--------------------------------------
M = 1.3608;
m = 0.096;
J = 0.0002175;
e = 0.0592;
k = 186.3;
thetat = pi/4;
d = (M+m)*(J + m*e^2)-(m*e*cos(thetat))^2;
%Linearized State Space Realization
D = 0;
A = computeA(m,J,e,thetat,k,d);
B = computeB(M,m,e,thetat,d);
% Original C = [1 0 0 0]
% Controllable but not observable
%Hence C is modified to multiple output to make the system observable
C = computeC();
%--------------------------------------
%In order to implement Luenberger, the system pair (A,C) must be
%observable
%Check if the system is controllable & observable
checkControlObsv(A,B,C)
%--------------------------------------
%desired eigenvalues for the gain and Luenberger observer
fprintf('\n ----------------------------------------------------------------------------------------------------\n')
desp_K = [-2+2i -2-2i -4+4i -4-4i]
desp_L = [-10-2i -10+2i -11-3i -11+3i]
[K,L] = computeKnL(A,B,C, desp_K, desp_L)
fprintf('\n Matrix K provides the gain matrix L provides the damping ratio required to implement the observer\n')
fprintf('\n ----------------------------------------------------------------------------------------------------\n')
%--------------------------------------
% Now build the model 8.5 pg. 326, with the alteration that we want both
% the true observation and the estimated observation:
[Aobc,Bobc,Cobc] = estimateed_Luenberger(A,B,C,K,L)
[n,p] = size(B);  
[m,n] = size(Cobc);
Dobc = zeros(m,p)
fprintf('\n ----------------------------------------------------------------------------------------------------\n')

%--------------------------------------
%Set up additional variables for the simulink model

%Use sine as a imput signal
%set the amplitude & phase
amp1 = 2; w1 = pi; 
%period = 2
% choose tfinal to get a good idea of what the system is doing:
tfinal = 10;

% to try random initial values: 
 x0_obc = [randn(n/2,1); zeros(n/2,1)];

 % If you want to simulate the uncompensated (original) system:
 x0 = x0_obc(1:n/2);

%---------------------------------------------------------
% Functions to make computations
function A = computeA(m,J,e,thetat,k,d)
A = [ 0 1 0 0; -k*(J+m*e^2)/d 0 0 0; 
      0 0 0 1; k*m*e*cos(thetat)/d 0 0 0];
end

function B = computeB(M,m,e,thetat,d)
B = [ 0 -m*e*cos(thetat)/d 0 (M+m)/d]';
end

function C = computeC()
C = [1 0 0 0;0 0 1 0];
end

function checkControlObsv(A,B,C)
fprintf('\n ----------------------------------------------------------------------------------------------------\n')
[v,d] = eig(A)
fprintf('\n The system has repeated eigenvalues and hence repeated eigenvectors. Hence, look at Jordan matrix\n')
[V,J]=jordan(A)
fprintf('\n Jordan blocks: [0 1; 0 0] & [-11.4i 0;0 +11.4i]\n')
fprintf('\n The angle response is expected to be oscillating\n')
%Checking for controability
a = (V^-1)*B
fprintf('\n The dominant matrix element n = 2 is not 0. Hence the system is controllable.\n')
%checking for observability
b = C*V
fprintf('\n The vector at the col = 1 is non-zero. Hence, the system is observable.\n')
fprintf('\n Hence, we can design a Luenberger observer for our system\n')
fprintf('\n ----------------------------------------------------------------------------------------------------\n')
end

function [K,L] = computeKnL(A,B,C, desp_K, desp_L)
K = place(A,B,desp_K);
L = place(A',C',desp_L)';
end

function [Aobc,Bobc,Cobc] = estimateed_Luenberger(A,B,C,K,L)
Aobc = [ A,  -B*K; L*C , A-B*K-L*C];
Bobc = [B;B];
Cobc = [C , zeros(size(C)); zeros(size(C)), C];
end
%--------------------------------------------------------