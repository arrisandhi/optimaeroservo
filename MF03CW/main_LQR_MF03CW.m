function [u,x,y,Klqr] = main_LQR_MF03CW(ssLinRed,StateWORBRed,SimDef,GustDef,CW)
% MPC set point tracking including input rate cost. This code performs tracking of piecewise constant
% reference signal(uref, xref). Cost function is a function of the
% difference between the value and the reference signal, and also
% difference between the current input value and the predicted input value.

% The arrangement of quadratic programming matrix is as follow:
% G = diag([RSQ RSQ RSQ ... RSP])
% other QP matrix follow the arrangement of G

%% PLANT MODEL
% Set discrete time MPC controller dynamics
SrA = ssLinRed.a;
SrB = ssLinRed.b;
SrC = ssLinRed.c;
SrD = ssLinRed.d;
% 
% SrBG = SrB(:,end);    % Gust input vector
% SrB = SrB(:,1:end-1);   % Exclude the gust input
% SrD = SrD(:,1:end-1);   % Exclude the gust input 

dt = SimDef.dt;           % Time step
tSim = SimDef.tSim;    % Simulation time
[Ad,Bd,Cd,~]=c2dm(SrA,SrB,SrC,SrD,dt); % Continous to discrete system

BdG = Bd(:,end);
Bd = Bd(:,1:end-1);

n = size(Ad,1);     % Ad row size 
m = size(Bd,2);     % Bd column size
p = size(Cd,1);     % Cd row size

%% SET UP QP
Q = CW(1)*diag(ones(1,n));      % State weight matrix
R = CW(2)*diag(ones(1,m));      % Input weight matrix  

%% SET UP DLQR
[Klqr,S,cl_poles]=dlqr(Ad,Bd,Q,R);


ul = -0.5*ones(m,1);
uh = 0.5*ones(m,1);

%% SIMULATION LOOP
x(:,1) = StateWORBRed;

u = zeros(m,length(tSim));
y = zeros(p,length(tSim));

Vgust_alltime = GustDef.Vgust_alltime;

for it = 1:numel(tSim);
    %% DLQR
    u(:,it) = -Klqr*x(:,it);
    
    % Introduce control saturation
    if u(1,it)> uh(1)
       u(1,it) = uh(1);
    elseif u(1,it) < ul(1)
       u(1,it) = ul(1);
    end
    
%     if u(2,it)> uh(2)
%        u(2,it) = uh(2);
%     elseif u(2,it) < ul(2)
%        u(2,it) = ul(2);
%     end
    
    x(:,it+1)=Ad*x(:,it)+Bd*u(:,it)+BdG*Vgust_alltime(3,it);
    y(:,it)=Cd*x(:,it);
   
end

    
    
 