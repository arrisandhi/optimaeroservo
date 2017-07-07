function [f0mat,fsmat,f0matr,fsmatr] = beamsectioncoeff_Arri (totals, compmat)
% Input :
% totals    : Length of beam section
% compmat   : Compliance matrix
%
% Output :
% f0mat     : Base force
% fsmat     : External applied force which will cause the said displacement
%             at tip
% f0matr    : Base force in reverse direction
% fsmatr    : External applied force in reverse direction

EMAT=[[     0     0     0      0      0       0     ];...
      [     0     0     0      0      0       0     ];...
      [     0     0     0      0      0       0     ];...
      [     0     0     0      0      0       0     ];...
      [     0     0    -1      0      0       0     ];...
      [     0     1     0      0      0       0     ]];
  

%========================= Define force at BASE ==========================%
% x2' + Ex2 = 0         (Eqn. 4.28a)
% COLUMN vectors representing x2_0
% Unity disturbance in each variable
x20 = eye(6);                             
x2s = -EMAT * x20;                        % Gradient of x2 

% COLUMN vectors representing x2_d
%       x2 =  x20 + s*x2s
% gm = cx2 = cx20 + s*cx2s
cx20 = compmat * x20;
cx2s = compmat * x2s;

% COLUMN vectors
% rotation = integral(gmn)ds
% rt = cx20*s + 1/2*cx2d*s*s
% rt =  rt0*s +      rts*s*s
% wrt s= 0;
rt0 = cx20(4:6,:);
rts = cx2s(4:6,:)/2;

% disp = integral(gm)ds + integral(rt)ds
% disp = s*cx20 + s*s*(cx2d/2+rt0~/2) + s*s*s*(rts~/3)
% disp_gm = s*cx20 + cx2d*s*s/2
% disp_rt = rt0~*s*s/2 + rts~*s*s*s/3                 
dp0  = cx20(1:3,:);

% rotation in z direction contributes to displacement in y direction
% rotation in y direction contributes to displacement in z direction
dps  = cx2s(1:3,:)/2+[zeros(1,6);rt0(3,:);-rt0(2,:)]/2;   
dpss = [zeros(1,6);rts(3,:);-rts(2,:)]/3;                 

% define totals
dpt = totals*dp0 + totals*totals*dps + totals*totals*totals*dpss;
rtt = totals*rt0 + totals*totals*rts;
totaldpmats = [dpt;rtt];

% tinv defines displacement to force
tinv = inv(totaldpmats);

% tinv*tipdisp = baseforce
% However base used to balance tip, therefore minus
f0mat = -tinv;

% tip force uses base force and gradient E
fsmat = tinv - EMAT*tinv*totals;
% These are EXTERNAL APPLIED FORCES which will cause the SAID DISPLACEMENT
% at the TIP

% for a negative totals, in case displacement is applied at the other side,
% ie s is -total
% while still fixed at s=0
dptr = -totals*dp0 + totals*totals*dps - totals*totals*totals*dpss;
rttr = -totals*rt0 + totals*totals*rts;
totaldpmatsr = [dptr;rttr];

% tinv defines displacement to force
tinvr = inv(totaldpmatsr);

% Note this time the reversed order of start and finish (ie start at -total
% and finish at 0)
fsmatr = tinvr;
f0matr = -(tinvr+EMAT*tinvr*totals);
% end