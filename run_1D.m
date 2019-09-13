%%

% Driver file for a single vessel model of pulse wave propagation.
% Model usines a simple linear elastic wall model and a three element Windkessl outflow
% boundary condition. The system is driven by imposing an inflow profile.

clc; clear all; close all;
!rm art*.2d
!rm Zhat
!make clean
!make
format shortg;
%% Define the parameters
f1 = 1e+5;
f2 = -25;
f3 = 9e+4;
HB = 10;
%% Create an Inflow waveform for pulse wave simlation in a single vessel

Qmin = 0;    % minumum flow rate during diastole (ml/s)
Qmax = 500;  % maximum flow rate during systole (ml/s)


t0 = 0.0; % initial time (s)
T = 0.85; % Length of the cardiac cycle (s)

Tm  = 0.15; % Peak systole time (s)
Td  = 0.4;  % End systole time (s)

tt = Td-Tm; % time difference between peak and end of systole

[Qin] = InFlow(Qmin,Qmax,t0,Tm,tt,T); % Generates a desired inflow profile and saves it to a file Qin_*.dat
%% Make the network
% Written for main, left and right pulmonary artery
L    = [2.87,4.74,6.23];
Rin  = [1.64,1.55,1.36];
Rout = [1.64,1.55,1.36];
dim_mat = [L; Rin; Rout]';

% order:  [parent d1 d2
%         d1    d3 d4 ...]
connectivity = [1 2 3];
terminal     = [2 3];

dlmwrite('connectivity.txt',connectivity-1,'\t');
dlmwrite('terminal_vessels.txt',terminal-1,'\t');
dlmwrite('Dimensions.txt',dim_mat,'\t');

num_ves = length(L);
num_term = length(terminal);


IMP_FLAG = 0;
[BC_matrix] = get_nominal_WK(IMP_FLAG,connectivity,terminal,L,Rin,Qin,f3);

%% call the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unix(sprintf('./sor06  %f %f %f %d %d %d',f1,f2,f3,HB,num_ves,num_term));
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
    
    Qmean_sim = [mean(Q_p) mean(Q_d)]
    Pmean_sim = [mean(P_p) mean(P_d)]
    
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
end
figure(12); hold on;
plot(linspace(0,1,length(Qin)),Qin,'k','LineWidth',3);
set(gca, 'fontsize',30);grid on;
xlabel 'Time (s)';
ylabel 'Flow (ml/s)';
legend('Prox','Mid','Dist','Data');

