function zDispWORBRed = calc_disprot(u,x,y,BeamDef,SimDef,trimstate,StateWORBRed,StateWORB,Nmodes,BeamSeg,Mode,Ti,RVnum)

%% CALCULATE ROTATION AND DISPLACEMENT ALONG SIMULATION TIME (REDUCED SYSTEM)
NumNode = BeamDef.NumNode; 
NumIndices = BeamDef.NumIndices; 
Nseg = BeamDef.Nseg; 
tSim = SimDef.tSim;
DispWORBRed = zeros(3,NumNode,numel(tSim));
RotWORBRed = zeros(3,3,NumNode,numel(tSim));

nInterest = 1;
ttime = 1;
CurrentState = trimstate;
for it = 1:numel(tSim)-50
    FullStateWORBRed = x(:,it);
    FullStateWORB = [FullStateWORBRed;zeros(numel(StateWORB)-numel(StateWORBRed),1)];
    FullStateWORB = Ti*FullStateWORB;
    CurrentState(1:RVnum,1)                         = FullStateWORB(1:RVnum);
    CurrentState((1:RVnum)+Nmodes+6,1)              = FullStateWORB(RVnum+(1:RVnum));
    CurrentState((1:RVnum)+2*Nmodes+6,1)            = FullStateWORB(2*RVnum+(1:RVnum));
    CurrentState((1:RVnum)+2*Nmodes+6+Nmodes+6,1)   = FullStateWORB(3*RVnum+(1:RVnum));

    [ModalRE, ModalCE] = dispintegrationvar_Arri(CurrentState, nInterest, ttime, NumNode, Nmodes, NumIndices, BeamSeg, Nseg, Mode);
    DispWORBRed(:,:,it) = ModalRE;
    RotWORBRed(:,:,:,it) = ModalCE;
    tSim(it) 
end
zDispWORBRed(1:numel(tSim)) = DispWORBRed(3,21,1:numel(tSim));

%% PLOT
% subplot(5,1,1);
% plot(tSim,y(1,:),'LineWidth',1.5); hold on;
% xlabel('Time (s)'); ylabel('M_y at root (Nm)');
% 
% subplot(5,1,2);
% plot(tSim,y(3,:),'LineWidth',1.5); hold on;
% xlabel('Time (s)'); ylabel('F_z at tip (N)');
% 
% subplot(5,1,3);
% plot(tSim,y(2,:),'LineWidth',1.5); hold on;
% xlabel('Time (s)'); ylabel('V_z at tip (m/s)');
% 
% % Plot z-displacement of the tip
% subplot(5,1,4);
% plot(tSim,zDispWORBRed,'LineWidth',1.5); hold on;
% xlabel('Time (s)'); ylabel('Disp_z at tip (m)');
% 
% subplot(5,1,5);
% plot(tSim,u(1,:),'LineWidth',1.5); hold on;
% xlabel('Time (s)'); ylabel('Control Deflection');


