%% PLOT
tSim = SimDef.tSim;
% figure;
subplot(3,1,1);
plot(tSim,y_Goland(1,:),'LineWidth',1.5); hold on;
xlabel('Time [s]'); ylabel('M_{y-root} [Nm]');

subplot(3,1,2);
plot(tSim,zDispWORBRed_Goland,'LineWidth',1.5); hold on;
xlabel('Time [s]'); ylabel('Disp_{z-tip} [m]');

subplot(3,1,3);
plot(tSim,u_Goland(1,:),'LineWidth',1.5); hold on;
xlabel('Time [s]'); ylabel('\delta_{ail} [rad]');