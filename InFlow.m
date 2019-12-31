%=============================================================*
%                                                             *
% InFlow.m                                                    *
% Version: 1.1 (created on 31 Dec. 2019)                      *
% AUTHORS: M.J. Colebank, M.U. Qureshi, M.S. Olufsen          *
%          Department of Mathematics                          *
%          North Carolina State University, Raleigh, USA      *
% DATE UPDATED: 31 Dec. 2019.                                 *
%                                                             *
% DESCRIPTION: This functions takes in a max and min flow     *
% value along with timing parameters to create a flow profile *
% which will drive the 1D model. This is only used if there   *
% is not a MRI/Echo obtained flow profile.
%=============================================================*
%%
function [Q,t] = InFlow(Qmin,Qmax, t0,Tm,tau,T )

% Fluids model typically deals with 2^N data points in inflow
N = 8192;

deltaT = (T-t0)/N;
t = linspace(0,T,N+1);
Q = zeros(1,N+1);

% Create the flow profile
for i = 1:N+1
    
    ti = i*deltaT;
    
    if (t0<=ti)&&(ti<=Tm)
        Q(i) = 0.5*(Qmax - Qmin)*(1 - cos(pi*ti/Tm)) + Qmin;
    elseif (Tm<=ti)&&(ti<=tau+Tm)
        Q(i) = 0.5*(Qmax - Qmin)*(1 + cos(pi*(ti-Tm)/tau)) + Qmin;
    elseif (tau+Tm<=ti)&&(ti<=T)
        Q(i) = Qmin;
        
    end
    
end
end

