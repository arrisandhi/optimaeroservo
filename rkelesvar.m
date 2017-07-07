function [Rfs, Cfs, dR1, dC1]=rkelesvar(R0,C0,gamma,kappa,e0,ds,dR0,dC0,dgamma,dkappa)
% Direct calculation of final rotation and displacement based on
% constant curvature and strain, without integration
% also computes variations from dR0, dC0, dgamma and dkappa
% Input:
%   R0      : Starting rotation 
%   C0      : Starting coordinate
%   gamma   : Strain mode in global axes
%   kappa   : Curvature mode in global axes
%   e0      : Running direction in global axes
%   ds      : Length of cell
%   dR0     : Starting rotation variation
%   dC0     : Starting coordinate varition
%   dgamma  : Strin mode variation
%   dkappa  : Curvature mode variation
%
% Output:
%   Rfs     : Displacement at the next node
%   Cfs     : Curvature at the next node
%   dR1     : Variation of displacement
%   dC1     : Variation of curvature


I = eye(3);
psi = kappa*ds;
phi = norm(psi);
if abs(phi) < 1e-6
    H0 = (I+(1-phi^2/6)*tilde(psi)+(0.5-phi^2/24)*tilde(psi)*tilde(psi));
    H1 = (I+(0.5-phi^2/24)*tilde(psi)+(1/6-phi^2/120)*tilde(psi)*tilde(psi))*ds;
else
    % Refer to Eqn. 9 in Structural and Aerodynamic Models in Nonlinear Flight 
    % Dynamics of Very Flexible Aircraft
    H0 = (I+sin(phi)/phi*tilde(psi)+(1-cos(phi))/phi^2*tilde(psi)*tilde(psi));
    H1 = (I+(1-cos(phi))/phi^2*tilde(psi)+(phi-sin(phi))/phi^3*tilde(psi)*tilde(psi))*ds;
end

% Refer to Eqn. 8 in Structural and Aerodynamic Models in Nonlinear Flight 
% Dynamics of Very Flexible Aircraft
Cfs = C0*H0;
Rfs = R0+C0*H1*(e0+gamma);

dH0 = zeros(3,3,3);
dH1 = zeros(3,3,3);
for i = 1:3
    if abs(phi) < 1e-6
        dH0(:,:,i) = ((1-phi^2/6)*tilde(I(:,i))+(-1/3+phi^2/30)*psi(i)*tilde(psi)+(0.5-phi^2/24)*(tilde(psi)*tilde(I(:,i))+tilde(I(:,i))*tilde(psi))+(-1/12+phi^2/180)*psi(i)*tilde(psi)*tilde(psi))*ds;
        dH1(:,:,i) = ((0.5-phi^2/24)*tilde(I(:,i))+(-1/12+phi^2/180)*psi(i)*tilde(psi)+(1/6-phi^2/120)*(tilde(psi)*tilde(I(:,i))+tilde(I(:,i))*tilde(psi))+(-1/60+phi^2/1260)*psi(i)*tilde(psi)*tilde(psi))*ds^2;
    else
        dH0(:,:,i) = (sin(phi)/phi*tilde(I(:,i))+(phi*cos(phi)-sin(phi))*psi(i)*tilde(psi)/phi^3+(1-cos(phi))/phi^2*(tilde(psi)*tilde(I(:,i))+tilde(I(:,i))*tilde(psi))+(phi*sin(phi)+2*(cos(phi)-1))/phi^4*psi(i)*tilde(psi)*tilde(psi))*ds;
        dH1(:,:,i) = ((1-cos(phi))/phi^2*tilde(I(:,i))+(phi*sin(phi)+2*(cos(phi)-1))/phi^4*psi(i)*tilde(psi)+(phi-sin(phi))/phi^3*(tilde(psi)*tilde(I(:,i))+tilde(I(:,i))*tilde(psi))+(3*sin(phi)-phi*cos(phi)-2*phi)/phi^5*psi(i)*tilde(psi)*tilde(psi))*ds^2;
    end
end

dH0g = zeros(3,3);
dH1g = zeros(3,3);
for i = 1:3
    dH0g = dH0g+dH0(:,:,i)*dkappa(i);
    dH1g = dH1g+dH1(:,:,i)*dkappa(i);
end

% Refer to Eqn. 8 in Structural and Aerodynamic Models in Nonlinear Flight 
% Dynamics of Very Flexible Aircraft. Chain rule!
dC1 = C0*dH0g+dC0*H0;
dR1 = dR0+dC0*H1*(e0+gamma)+C0*dH1g*(e0+gamma)+C0*H1*dgamma;

return