% Alex Howard
% MECH 4175
% HW 5 Prob 2
clc
clear all
close all

q = -1; % lb/in
w = 12; %in
h = 1; % in
L = 100; % in
E = 10e7; % psi
I = 1/12*w*h^3;

% Element matrices
k1 = E*I/L^3* [12    6*L     -12    6*L;
               6*L   4*L^2   -6*L   2*L^2;
               -12   -6*L    12     -6*L;
               6*L   2*L^2   -6*L   4*L^2];
k2 = k1;

% Load vector for element 1
r_q1 = q*L/2*[1; L/6; 1; -L/6];

% Build global K matrix
kglobal = zeros(6)
dof1 = [1 2 3 4];
dof2 = [3 4 5 6];

for i=1:4
    for j=1:4
        kglobal(dof1(i),dof1(j))=kglobal(dof1(i),dof1(j))+k1(i,j);
        kglobal(dof2(i),dof2(j))=kglobal(dof2(i),dof2(j))+k2(i,j);
    end
end

% Build global force vector
%   no point loads
r = zeros(6,1); 
for i = 1:4
    r(dof1(i)) = r(dof1(i))+r_q1(i); % add in element 1 load vector
end

% Apply essential boundary conditions
%   v_1 = 0, \theta_1 = 0, v_3 = 0, \theta_3 = 0
%   rmod remains the same b/c all ebc's = 0
kmod = kglobal;
rmod = r;

ebc = [1 2 5 6]; % EBC's at these locations
mu = max(max(kglobal))*1e8;

for i = 1:4
    kmod(ebc(i),ebc(i)) = kmod(ebc(i),ebc(i)) + mu;
end

d = kmod\rmod % solve the system

%% find and plot displacement, moment, whatever
s = 0:1:100;

% % Element 1
s1_3 = -2/L^3*d(3)+1/L^2*d(4);      % s^3 Terms
s1_2 = 3/L^2*d(3)-1/L*d(4);         % s^2
v_1 = s1_3*s.^3+s1_2*s.^2;          % displacement
M_1 = E*I*(6*s1_3.*s+2*s1_2);       % Bending Moment
V_1 = E*I*6*s1_3*ones(1,length(s)); % Shear Force


% Repeat for element 2
s2_3 = 2/L^3*d(3)+1/L^2*d(4);
s2_2 = -(3/L^2*d(3)+2/L*d(4));
s2_1 = d(4);
s2_0 = d(3);
v_2 = s2_3*s.^3 + s2_2*s.^2 + s2_1.*s + s2_0;
M_2 = E*I*(6*s2_3.*s+2*s2_2);
V_2 = E*I*6*s2_3*ones(1,length(s));

% plotz
plot(s,v_1,'r')
hold on
plot(s+L,v_2,'b')
title('Displacement Plot')
xlabel('Location ''s'' [in]')
ylabel('displacement [in]')
legend('Element 1','Element 2')

figure(2);
plot(s,M_1,'r');
hold on
plot(s+L,M_2,'b');
plot([100 100],[M_1(end) M_2(1)],'k--');
title('Bending Moment Plot')
xlabel('Location ''s'' [in]')
ylabel('Bending Moment [lbf*in]')
legend('Element 1','Element 2','Discontinuity')

figure(3)
plot(s,V_1,'r')
hold on
plot(s+L,V_2,'b')
plot([100 100],[V_1(end) V_2(1)],'k--');
title('Shear Force Plot')
xlabel('Location ''s'' [in]')
ylabel('Shear [lbf]')
legend('Element 1','Element 2','Discontinuity')