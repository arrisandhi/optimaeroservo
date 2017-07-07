function [history,searchdir] = runfmincon_MA_TM_MMax(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,X0,BaseRA)
 fprintf(1, 'Run optimisation process.\n');
 
% Set up shared variables with OUTFUN
history.x    = [];
history.fval = [];
searchdir    = [];


% Fmincon configuration
Ndivs = BeamDef.Ndivs;
LB = 5*ones(Ndivs,1);
UB = 100*ones(Ndivs,1);

Fun = @(X)main_optimCost_MA_TM_MMax(X,BeamDef,BaseRA);
NLConstraintFun = @(X)confun_MA_TM_MMax(X,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,maxMomOri);
options = optimoptions(@fmincon,'OutputFcn',@outfun,'UseParallel',true,'Display','iter','Algorithm','active-set');

X = fmincon(Fun,X0,[],[],[],[],LB,UB,NLConstraintFun,options);

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