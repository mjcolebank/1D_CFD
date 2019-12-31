%=============================================================*
%                                                             *
% get_nominal_WK.m                                            *
% Version: 1.1 (created on 31 Dec. 2019)                      *
% AUTHORS: M.J. Colebank, M.U. Qureshi, M.S. Olufsen          *
%          Department of Mathematics                          *
%          North Carolina State University, Raleigh, USA      *
% DATE UPDATED: 31 Dec. 2019.                                 *
%                                                             *
% DESCRIPTION: This function takes in various inputs and      *
% calculates the nominal Windkessel boundary conditions for   *
% each terminal vessel. Subroutines to calcuate the diastolic *
% decay constant are also included.                           *
%=============================================================*


function [BC_matrix] = get_nominal_WK(IMP_FLAG,conn,term,L,R,Qdat,Pdat,kvals,t)
%% Initialize all the Windkessel variables
num_vessels = length(L);
r4l     = zeros(num_vessels,1);
Q       = zeros(num_vessels,1);
Rtotal  = zeros(num_vessels,1);
R1      = zeros(num_vessels,1);
R2      = zeros(num_vessels,1);
CT      = zeros(num_vessels,1);
Zc      = zeros(num_vessels,1);


%% Perform a nondimensionalization on all the input parameters
% To calculate nondimensional mean pressure and tau parameter, we need to
% consider the two different scenarios when 1) we have dynamic pressure
% data, or 2) we have only static values.
if length(Pdat)>3
    % Determine the Tau parameter in the Windkessel using the available
    % time-series data.
    P_spread = mean(Pdat);
    [tau,~] = TimeConstant_timeseries(Pdat,Qdat,t);
else
    % Determine the Tau parameter in the Windkessel using the available
    % static data values.
    P_spread= Pdat(2);
    [tau,~] = TimeConstant_static(Pdat,Qdat,t);
end

% Nondimensional all constants
Q_spread = mean(Qdat);
cf = 1332.22;                % Conversion factor from mmHg to g/cm/s^2
g = 981;                     % Gravitational constant (g/cm^2)
ratio = 0.2;                 % Windkessel Ratio (R1/RT) when not using Zc impedance
rho = 1.055;                 % Density of blood, assumed constant (g/cm^3)
Lr = 1.0;                    % Nondimensional length
qc = 10*Lr^2;                % Nondimensional flow
Pnd = P_spread*(cf/rho/g/Lr);% Nondimensional mean pressure
Qnd = Q_spread/qc;           % Nondimensional mean flow
tc = tau*qc/Lr^3;            % Nondimensional tau parameter for WK

% Define the radius dependent stiffness and nondimensionalize it.
stiff = @(r) (kvals(1).*exp(kvals(2).*r) + kvals(3)).*(4/3).*(rho/g/Lr);
%%
% To calculate the nominal Windkessel resistance values, we follow a
% Poiseuille relationship. More on this can be found in publications by
% Colebank, Qureshi, and Olufsen (2018-2019)
for i=1:num_vessels
    r4l(i) = (R(i).^4)./L(i);
end
conn_inc = 1;
% If any parameters are zero or INF, its because those vessels are unused
for i=conn(:,1)'
    if i==1
        Q(i) = Qnd;
    end
    if any(conn(:,1) == i) && size(conn,2) > 1
        d1 = conn(conn_inc,2); %Assign Daughters
        d2 = conn(conn_inc,3); %Assign Daughters
        Q(d1) = Q(i)*r4l(d1)./(r4l(d1) + r4l(d2)); % Distribute flow according
        Q(d2) = Q(i)*r4l(d2)./(r4l(d1) + r4l(d2)); % to Poiseuille relation
        conn_inc = conn_inc + 1;
    end
end

% Here, we calculate RT = P/Q, and distribute the resistance. Note that we
% set a 20/80 rule for R1/R2 when characteristic impedance is not used. If
% you do use the Zc values, make sure they are NOT negative. If the
% stiffness values are too large, Zc can become greater than RT (hence
% making R2<0). Changing alpha above will change the 20/80 rule.
for i=1:num_vessels
    Rtotal(i) = Pnd./Q(i);
    R1(i) = round(ratio.*Rtotal(i),4);
    R2(i) = round(Rtotal(i)-R1(i),4);
    CT(i) = tc./Rtotal(i);
    Zc(i) = (1./(pi.*R(i).^2)).*sqrt(rho.*stiff(R(i))./2);
end

% Set up all the boundary condition structures
param = zeros(num_vessels.*2 + length(term).*3,1);
inneri=1;
for i=1:2:2*num_vessels
    param(i) = L(inneri);
    param(i+1) = R(inneri);
    inneri = inneri + 1;
end
inneri = 1;
BC_matrix = zeros(length(term),3);

% Decide whether to use alpha/(1-alpha) rule or use Zc.
if IMP_FLAG == 1
    R1 = Zc;
    R2 = Rtotal - Zc;
end
for i=2*num_vessels+1:3:length(param)
    param(i)   = round(R1(term(inneri)),2);
    param(i+1) = round(R2(term(inneri)),2);
    param(i+2) = round(CT(term(inneri)),6);
    BC_matrix(inneri,:) = param(i:i+2,1);
    inneri = inneri+1;
end

% Write to file for the C++ code.
format longg
dlmwrite('Windkessel_Parameters.txt',BC_matrix,'delimiter','\t','precision','%10.6f');


end

% This function will take the time-series pressure and flow data and
% estimate the Tau parameter. This considers two different scenarios, one
% where the capillary pressure is zero (tau1) and one where the capillary
% pressure is nonzero (tau2). For the time being, we only use tau1 in the
% code.
%
% Originally written by M. Umar Qureshi
% Modified by MJ Colebank.
function [tau,tau2] = TimeConstant_timeseries(p,q,t)

% global p q N t t0 p0 tdia pdia nc n;
N = length(t);

dq = diff(q);
dt = diff(t);
dqdt = dq./dt;
nc = 3; 
m = floor((N-1)/nc);
id = find(dqdt==0);

for i = 1:length(id)
    idtest = id(i);
    if(idtest > m)
        id1(i)=idtest;
    else
    end
end

id1 = id1(id1~=0);
ids = id1(1)-n;

pdia = p(ids:end);
tdia = t(ids:end);

t0 = tdia(1);
p0 = pdia(1);

%% Define the two different exponentials

tau_opt1 = @(x) sum((pdia - (p0*exp(-(tdia-t0)/x))).^2);
tau_opt2 = @(x) sum((pdia - (x(2)+(p0-x(2))*exp(-(tdia-t0)/x(1)))).^2);


%%
options = optimset('MaxFunEvals',1e8,'TolX',1e-8,'MaxIter',15000);

x0 = 0.1;
[xopt1,~,~,~] = fminsearch(tau_opt,x0,options);
myfit1 = p0*exp(-(t-t0)/xopt1);


x0 = [0.1 0];
[xopt2,~,~,~] = fminsearch(tau_opt2,x0,options);
myfit2 = xopt2(2)+(p0-xopt2(2))*exp(-(t-t0)/xopt2(1));


figure; clf;
plot(t,p,'-r',t,myfit1,'-k',t,myfit2,'--c');hold on;plot(t0,p0,'ob');
plot(t,min(p)+(max(p)-min(p))*(q-min(q))/(max(q)-min(q)),'b');
line([t0 t0],[5 p0],'LineStyle','--','Color',[0 0 0]);
set(gca,'ylim',[5 20]);

tau = xopt1;
tau2 = xopt2(1);
end


% This function will take the static pressure and dynamic flow data and
% estimate the Tau parameter. This considers two different scenarios, one
% where the capillary pressure is zero (tau1) and one where the capillary
% pressure is nonzero (tau2). For the time being, we only use tau1 in the
% code.
%
% Originally written by M. Umar Qureshi
% Modified by MJ Colebank.
function [tau,tau2] = TimeConstant_static(p,q,t)

% global p q N t t0 p0 tdia pdia nc n;
N = length(t);

dq = diff(q);
dt = diff(t);
dqdt = dq./dt;
nc = 3; n=0.0;%NOT SURE WHY: UMAR SET THIS

m = floor((N-1)/nc);
id = find(dqdt==0);

for i = 1:length(id)
    idtest = id(i);
    if(idtest > m)
        id1(i)=idtest;
    else
    end
end

id1 = id1(id1~=0);
ids = id1(1)-n;

% Try to fit the exponential between the mean pressure and the diastolic
% value
pdia = p(2:3);
tdia = t([ids(1) end]);

t0 = tdia(1);
p0 = pdia(1);

%% Define the two different exponentials

tau_opt1 = @(x) sum((pdia - (p0*exp(-(tdia-t0)/x))).^2);
tau_opt2 = @(x) sum((pdia - (x(2)+(p0-x(2))*exp(-(tdia-t0)/x(1)))).^2);


%%
options = optimset('MaxFunEvals',1e8,'TolX',1e-8,'MaxIter',15000);

x0 = 0.1;
[xopt1,~,~,~] = fminsearch(tau_opt1,x0,options);
myfit1 = p0*exp(-(t-t0)/xopt1);


x0 = [0.1 0];
[xopt2,~,~,~] = fminsearch(tau_opt2,x0,options);
myfit2 = xopt2(2)+(p0-xopt2(2))*exp(-(t-t0)/xopt2(1));

figure; clf;
plot(tdia,pdia,'-r',t,myfit1,'-k',t,myfit2,'--c');hold on;plot(t0,p0,'ob');
% plot(t,min(p)+(max(p)-min(p))*(q-min(q))/(max(q)-min(q)),'b');
line([t0 t0],[5 p0],'LineStyle','--','Color',[0 0 0]);
set(gca,'ylim',[5 20]);

tau = xopt1;
tau2 = xopt2(1);
end
