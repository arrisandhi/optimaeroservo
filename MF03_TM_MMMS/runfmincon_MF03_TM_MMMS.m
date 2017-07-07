function [history,searchdir] = runfmincon_MF03_TM_MMMS(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,X0,BaseRA)
 fprintf(1, 'Run optimisation process.\n');
 
% Set up shared variables with OUTFUN
history.x    = [];
history.fval = [];
searchdir    = [];

% Fmincon configuration
LB = [-10000000000000;-10000000000000;-10000000000000];          % Lower bound of design vector
UB = [10000000000000;10000000000000;10000000000000];         % Upper bound of design vector

Fun = @(X)main_optimCost_MF03_TM_MMMS(X,BeamDef,BaseRA);
NLConstraintFun = @(X)confun_MF03_TM_MMMS(X,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,maxMomOri);
options = optimoptions(@fmincon,'OutputFcn',@outfun,'UseParallel',true,'Display','iter','Algorithm','active-set');

Atemp = [BeamDef.NodeL(1:BeamDef.Ndivs).^2, BeamDef.NodeL(1:BeamDef.Ndivs), ones(BeamDef.Ndivs,1)];
A = [-Atemp; Atemp];
b = [-5*ones(BeamDef.Ndivs,1);100*ones(BeamDef.Ndivs,1)];

X = fmincon(Fun,X0,A,b,[],[],LB,UB,NLConstraintFun,options);

 function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
%              hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval]
           history.x = [history.x, x]
           fprintf('%d\n',history.fval);
           fprintf('%d\n',history.x);
         % Concatenate current search direction with 
         % searchdir.
%            searchdir = [searchdir;... 
%                         optimValues.searchdirection'];
%            plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,x(2),... 
%                 num2str(optimValues.iteration));
%            title('Sequence of Points Computed by fmincon');
         case 'done'
%              hold off
         otherwise
     end
 end
end