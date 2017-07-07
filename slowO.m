function [a11h, v1, s1, t1] = slowO(SSS, cut)
%SLOWFAST State-space slow-fast decomposition.
%
% [GS,GF] = SLOWFAST(G,CUT)  produces a decomposition of LTI object G
% into the sum of a slow part GS and a fast part GF; CUT is the order of GS. 
%
%                 G(s) = Gs(s) + Gf(s)

%old help
%SLOWFAST State-space slow-fast decomposition.
%
% [SS_1,SS_2] = SLOWFAST(SS_,CUT) or
% [A11H,B1H,C1H,D1H,A22H,B2H,C2H,D2H] = SLOWFAST(A,B,C,D,CUT) produces
%     a decomposition of G(s) into the sum of a slow part Gs(s) and a
%     fast part Gf(s)
%
%                 G(s) = Gs(s) + Gf(s)
%     where
%                 Gs(s):= ss_1 = mksys(a11h,b1h,c1h,d1h);
%                 cut = no. of modes in Gs(s)
%                 Gf(s):= ss_2 = mksys(a22h,b2h,c2h,d2h);
%
%   The regular state-space can be recovered by "branch".

% R. Y. Chiang & M. G. Safonov 4/4/86
% Copyright 1988-2004 The MathWorks, Inc.
%       $Revision: 1.1.6.2 $
% All Rights Reserved.
% -----------------------------------------------------------------------

disp('  ');
disp('        - - Working on slow/fast decomposition - -');

% nag1=nargin;
% [emsg,nag1,xsflag,Ts,a,b,c,d,cut]=mkargs5x('ss',varargin); error(emsg);
a=SSS.a;
b=SSS.b;
c=SSS.c;
d=SSS.d;


%
% ----- Real Schur Decomposition :
%
[ra,ca] = size(a);
[u,t,m] = blkrsch(a,6,cut);
% T=V'AV
% a11h=V'AV(1:cut first diag)=v1Av2
% v1=V(:,1:cut), less cols than V
%
a11h = t(1:cut,1:cut);
a22h = t(cut+1:ra,cut+1:ra);
a12h = t(1:cut,cut+1:ra);
v1 = u(:,1:cut); v2 = u(:,cut+1:ra);
%
x = lyap(a11h,-a22h,a12h);
%
t1 = v1;
t2 = v1 * x + v2;
s1 = v1' - x * v2';
s2 = v2';
b1h = s1 * b;
b2h = s2 * b;
c1h = c * t1;
c2h = c * t2;
d1h = 0.5*d;
d2h = 0.5*d;
% B1=s1*B
% C1=C*t1

   a11h = ss(a11h,b1h,c1h,d1h);
   % b1h = mksys(a22h,b2h,c2h,d2h);
end
%
% ------ End of SLOWFAST.M ---- RYC/MGS ---- %