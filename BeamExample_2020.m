% BeamExample_2020.m
%--------------------------------------------------------------------------

E = 210e6; % Young's modulus in kN/m^2
L = 2; % Length of each element in m
I = 4e-4; % Moment of inertia in m^4
q = -10; % Distributed load in kN/m.
F = -18; % Concentrated force in kN

% Set up a general beam element stiffness matrix
ke = [12 6*L -12 6*L;
      6*L 4*L^2 -6*L 2*L^2;
      -12 -6*L 12 -6*L;
      6*L 2*L^2 -6*L 4*L^2];

% Element-specific matrices
k1 = ((2*E*I)/L^3)*ke
k2 = k1
k3 = ((E*I)/(L^3))*ke

pause;

% Element load vector for element 3
rq3 = ((q*L)/2)*[1;L/6;1;-L/6]

pause;

% Assemble global system
kglobal = zeros(8,8);
dof1 = [1 2 3 4];
dof2 = [3 4 5 6];
dof3 = [5 6 7 8];
for i=1:4
    for j=1:4
        kglobal(dof1(i),dof1(j))=kglobal(dof1(i),dof1(j))+k1(i,j);
        kglobal(dof2(i),dof2(j))=kglobal(dof2(i),dof2(j))+k2(i,j);
        kglobal(dof3(i),dof3(j))=kglobal(dof3(i),dof3(j))+k3(i,j);
    end
end
kglobal

pause;

% Assemble global force vector
r = zeros(8,1);
r(3) = F; % Concentrated force at node 2
for i=1:4 % For loop for entering rq3 into global vector
    r(dof3(i))=r(dof3(i)) + rq3(i);
end
r

pause;

% Apply essential boundary conditions using the penalty method
% v1 = 0, theta1 = 0, and v4 = 0

kMod = kglobal;
rMod = r;

mu = max(max(abs(kMod)))*10^8

% Apply v1 = 0
kMod(1,1) = kMod(1,1) + mu;

% Apply theta1 = 0
kMod(2,2) = kMod(2,2)+mu;

% Apply v4 = 0
kMod(7,7) = kMod(7,7)+mu;

kMod
rMod

pause;

% Solve for unknown displacements and rotations
d = kMod\rMod

pause;

s = linspace(0,2); % s values go from 0 to 2 for each element

% Get v(s), V(s), and M(s) for element 1
s31 = d(3)*-0.25+d(4)*0.25; % Coefficient for s^3 terms
s21 = d(3)*(3/4)+d(4)*(-1/2); % Coefficient for s^2 terms
v1 = s31*s.^3+s21*s.^2; % Displacement along element 1
M1 = 2*E*I*(6*s31*s + 2*s21); % Moment in element 1
% Note the 2EI term in elements 1 and 2
V1 = 2*E*I*6*s31*ones(1,length(s)); % Shear in element 1
 
% Get v(s), V(s), and M(s) for element 2
s32 = d(3)*0.25+d(4)*0.25+d(5)*-0.25+d(6)*0.25;% Coefficient for s^3 terms
s22 = d(3)*(-3/4)+d(4)*(-1)+d(5)*(3/4)+d(6)*(-0.5); % Coefficient for s^2 terms
s12 = d(4)*1; % Coefficient for s terms
s02 = d(3)*1; % Constant
v2 = s32*s.^3+s22*s.^2+s12*s+s02; % Displacement along element 2
M2 = 2*E*I*(6*s32*s+2*s22); % Moment in element 2
V2 = 2*E*I*6*s32*ones(1,length(s)); % Shear in element 2

% Get v(s), V(s), and M(s) for element 3
s33 = d(5)*(1/4)+d(6)*(1/4)+d(8)*(1/4); % Coefficient for s^3 terms
s23 = d(5)*(-3/4)+d(6)*(-1)+d(8)*(-1/2); % Coefficient for s^2 terms
s13 = d(6)*1;% Coefficient for s terms
s03 = d(5)*1; % Constant
v3 = s33*s.^3+s23*s.^2+s13*s+s03; % Displacement along element 3
M3 = E*I*(6*s33*s + 2*s23); % Moment in element 3
V3 = E*I*6*s33*ones(1,length(s)); % Shear in element 3

% Plot displacement
plot(s,v1);
hold on;
plot(s+2,v2); % Shift to plot in correct place along the beam
plot(s+4,v3); % Shift to plot in correct place along the beam
plot(linspace(0,6),0,'--'); % Plot x-axis
xlabel('Distance Along Beam (m)');
ylabel('Vertical Displacement (m)');

pause;

% Plot moment
figure(2);
plot(s,M1);
hold on;
plot(s+2,M2); % Shift to plot in correct place along the beam
plot([4 4],[M2(length(M2)) M3(1)]); % Plot discontinuity
plot(s+4,M3); % Shift to plot in correct place along the beam
plot(linspace(0,6),0,'--'); % Plot x-axis
xlabel('Distance Along Beam (m)');
ylabel('Moment (kN*m)');

pause;

% Plot shear
figure(3);
plot(s,V1);
hold on;
plot([2 2],[V1(length(V1)) V2(1)]); % Plot discontinuity
plot(s+2,V2);% Shift to plot in correct place along the beam
plot([4 4],[V2(length(V2)) V3(1)]); % Plot discontinuity
plot(s+4,V3);% Shift to plot in correct place along the beam
plot(linspace(0,6),0,'--'); % Plot x-axis
xlabel('Distance Along Beam (m)');
ylabel('Shear (kN)');



