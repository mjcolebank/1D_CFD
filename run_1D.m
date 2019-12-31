%=============================================================*
%                                                             *
% run_1D.m                                                    *
% Version: 1.1 (created on 31 Dec. 2019)                      *
% AUTHORS: M.J. Colebank,M.U. Qureshi,  M.S. Olufsen          *
%          Department of Mathematics                          *
%          North Carolina State University, Raleigh, USA      *
% DATE UPDATED: 31 Dec. 2019.                                 *
%                                                             *
% DESCRIPTION: This script creates the inteface with the C++  *
% code for the 1D fluid dynamics model. The script creates    * 
% and runs the executable by passing selected model           *
% parameters and plots the hemodynamic waveforms              *
%                                                             *
%=============================================================*



% Model usines a simple linear elastic wall model and a three element Windkessl outflow
% boundary condition. The system is driven by imposing an inflow profile.

clc; clear all; close all;
!chmod +x sor06
!make clean
!make
format shortg;
%% Define the parameters
% The stiffness is assumed to follow an exponential curve, i.e.
% Eh/r0 = k1*exp(k2*r0) + k3
% where E, h, and r0 are the Youngs Modulus, wall thickness, and reference
% radius, respectively. To set a single stiffness value for all the
% vessels, set k1=0;
k1 = 1e+5;
k2 = -25;
k3 = 9e+4;
%% Create an Inflow waveform for pulse wave simlation in a single vessel
% Comment this section out if you have a flow waveform you wish to provide
% as an input to the model
Qmin = 0;    % minumum flow rate during diastole (ml/s)
Qmax = 500;  % maximum flow rate during systole (ml/s)


t0 = 0.0; % initial time (s)
T  = 0.85; % Length of the cardiac cycle (s)

Tm  = 0.15; % Peak systole time (s)
Td  = 0.4;  % End systole time (s)

tt = Td-Tm; % time difference between peak and end of systole

[Qin,tQ] = InFlow(Qmin,Qmax,t0,Tm,tt,T); % Generates a desired inflow profile and saves it to a file Qin_*.dat
% OR load in a flow waveform
% Qin = load('MYFLOWPROFILE.dat');
fname = strcat('Qin.dat');
dlmwrite(fname,Qin');

%% Load in a pressure profile
% This code was developed for parameter estimation purposes. If no pressure
% profile is available, simply provide a systolic, diastolic, and mean
% pressure value below
% Pdat = load('MYPRESSUREPROFILE.dat');
% Psys = max(Pdat); Pdia = min(Pdat); Pmean = mean(Pdat);

% Or hard code the systolic, diastolic, and mean value
Psys = 30; Pdia = 8; Pmean = (Psys+2*Pdia)/3;
Pdat = [Psys, Pmean, Pdia];
%% Specify vessel dimensions
% Written for main, left and right pulmonary artery
% Note that this can also be adapted to read in text files if you choose.
% Length and radius arrays start with vessel 1 and move up till vessel N
L    = [2.87,4.74,6.23];

% Rin = Inlet radius, Rout = Outlet radius
% The code can handle tapering if that information is available; simply
% change Rout to Rin*taper or hardcode the outlet radius
Rin  = [1.64,1.55,1.36];
Rout = [1.64,1.45, 1.20];


dim_mat = [L; Rin; Rout]';
%% Define the connectivity matrix
% order:  [parent d1 d2
%          d1    d3 d4 
%          d2    d5 d6 etc ...]
% NOTE: it is important that the largest vessel 'ids' are located at the
% end of the connectivity matrix. For example, 
%
%   conn = [1 2 3; 3 6 7; 2 4 5]
%
% will NOT work because of how the C++ code declares the artery class. Once
% this connectivity matrix is passed to the C++ code, it will take each
% "vessel ID" in the connectivity matrix and associate it with its
% parent/daughters and with its length and radius specified in "dim_mat".

connectivity = [1 2 3];
terminal     = [2 3];

% Write all this information to file. We subtract one to keep up with C++
% indexing.
dlmwrite('connectivity.txt',connectivity-1,'\t');
dlmwrite('terminal_vessels.txt',terminal-1,'\t');
dlmwrite('Dimensions.txt',dim_mat,'\t');

num_ves = length(L);
num_term = length(terminal);
num_pts  = 6;
%% Specify Windkessel boundary
% Set IMP_FLAG=1 if you want to use the "characteristic impedance" for the
% proximal resistance value. Note that this can negative resistance using
% the method here for Windkessel parameterization.
IMP_FLAG = 0;
kvals = [k1 k2 k3];
[BC_matrix] = get_nominal_WK(IMP_FLAG,connectivity,terminal,L,Rout,Qin,Pdat,kvals,tQ);

%% call the model
pars = [k1,k2,k3,T,num_ves,num_term,num_pts];
pars_str = mat2str(pars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unix(sprintf('./sor06  %s',pars_str(2:end-1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting Simulated data
%


for i=1:3
    data = load (strcat('art',num2str(i),'_1.2d'));
    
    [time,x,p,q,A,C] = gnuplot(data);
    
    P_p = p(:,1);  Q_p = q(:,1);   A_p = A(:,1); C_p = C(:,1);
    P_m = p(:,floor((end)/2));  Q_m = q(:,floor((end)/2));   A_m = A(:,floor((end)/2)); C_m = C(:,floor((end)/2));
    P_d = p(:,end); Q_d = q(:,end); A_d = A(:,end); C_d = C(:,end);
    t = time(:,1)-time(1,1);
    
    
    figure(1+i*10); hold on;
    plot(t,P_p,':','linewidth',3);
    plot(t,P_m,':','linewidth',3);
    plot(t,P_d,':','linewidth',3)
    hold off;
    set(gca, 'fontsize',30);grid on;
    xlabel 'Time (s)';
    ylabel 'Pressure (mmHg)';
    legend('Prox','Mid','Dist');
    
    figure(2+i*10); hold on;
    plot(t,Q_p,':','linewidth',3);
    plot(t,Q_m,':','linewidth',3);
    plot(t,Q_d,':','linewidth',3)
    set(gca, 'fontsize',30);grid on;
    xlabel 'Time (s)';
    ylabel 'Flow (ml/s)';
    legend('Prox','Mid','Dist');
    
    
    % Can plot a 3D representation
    figure(1+i*100); hold on;
    mesh(x,t,p);
    
    figure(2+i*100); hold on;
    mesh(x,t,q);
end
figure(12); hold on;
plot(tQ,Qin,'k','LineWidth',3);
set(gca, 'fontsize',30);grid on;
xlabel 'Time (s)';
ylabel 'Flow (ml/s)';
legend('Prox','Mid','Dist','Data');

