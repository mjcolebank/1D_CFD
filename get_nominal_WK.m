function [BC_matrix] = get_nominal_WK(IMP_FLAG,conn,term,L,R,Qdat,Ehr0)
%%
num_vessels = length(L);
r4l     = zeros(num_vessels,1);
Q       = zeros(num_vessels,1);
Rtotal  = zeros(num_vessels,1);
R1      = zeros(num_vessels,1);
R2      = zeros(num_vessels,1);
CT      = zeros(num_vessels,1);
Zc      = zeros(num_vessels,1);


%% Define parameter scaling from WKnominal
Q_spread = mean(Qdat);
P_spread = 14; %Mean pressure for healthy case
con_factor = 1332.22;
g = 981;
ratio = 0.2;
tau = 0.8833; % This parameter should be measured from flow/pressure waveform
rho = 1.055;
Lr = 1.0;
qc = 10*Lr^2;
Pnd = P_spread*(con_factor/rho/g/Lr);
Qnd = Q_spread/qc;
tc = tau*qc/Lr^3;

for i=1:num_vessels
    r4l(i) = (R(i).^4)./L(i);
end
conn_inc = 1;
%%If any parameters are zero or INF, its because those vessels are unused
for i=conn(:,1)'
    if i==1
        Q(i) = Qnd;
    end
    if any(conn(:,1) == i) && size(conn,2) > 1
        d1 = conn(conn_inc,2); %Assign Daughters
        d2 = conn(conn_inc,3);
        Q(d1) = Q(i)*r4l(d1)./(r4l(d1) + r4l(d2));
        Q(d2) = Q(i)*r4l(d2)./(r4l(d1) + r4l(d2));
        conn_inc = conn_inc + 1;
    end
end
for i=1:num_vessels
    Rtotal(i) = Pnd./Q(i);
    R1(i) = round(ratio.*Rtotal(i),4);
    R2(i) = round(Rtotal(i)-R1(i),4);
    CT(i) = tc./Rtotal(i);
    Zc(i) = (1./(pi.*R(i).^2)).*sqrt(2.*rho.*Ehr0.*(rho/g/Lr)./3);
end
param = zeros(num_vessels.*2 + length(term).*3,1);
inneri=1;
for i=1:2:2*num_vessels
    param(i) = L(inneri);
    param(i+1) = R(inneri);
    inneri = inneri + 1;
end
inneri = 1;
BC_matrix = zeros(length(term),3);
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
format longg
dlmwrite('Windkessel_Parameters.txt',BC_matrix,'delimiter','\t','precision','%10.6f');


end