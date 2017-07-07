function [ Rout, Cout ] = dispintt( x, Rinitial, Cinitial )
% Integrates displacement from intrinsic solution with t. uses x, cmat, C0 
% and gammainit as starting conditions and also n and dl as settings

% x is now velocity instead of strains, 2 columns with start and end values
% rkele remains unchanged

% global settings
global dtt

% These are to be fed to the RK solver, not fed through from main program
global gamma0
global gamma1
global kp0
global kp1

global tcurrent

    % Initialise
           
    Ra=zeros(3,2);
    % gamma and kappa is in LOCAL frame
    Ra=Rinitial;
    % C0 at base is identity
    CaB=Cinitial;
    % only run once
    g0=x(:,1);
    g1=x(:,2);    
    gamma0=g0(1:3);
    gamma1=g1(1:3);
    kp0=g0(4:6);
    kp1=g1(4:6);
    % define c matrix
    Cinit=CaB(1:3,1:3);
    % define initial u
    Rinit=Ra;
    % finalise input y
    yvec=[Rinit;Cinit(1,:)';Cinit(2,:)';Cinit(3,:)'];
    % start R-K
    [tout,vout]=ode45(@rkelet,[0,dtt],yvec);
    yfin=vout(end,:)';
    CaB(:,:)=[yfin(4:6)';yfin(7:9)';yfin(10:12)'];
    VAtmp=mattoangle(CaB(:,:));
    CaB(:,:)=angletomat(VAtmp);
    Rout=yfin(1:3);
    % Cout=(CaB(:,:))';
    Cout=(CaB(:,:));


end

