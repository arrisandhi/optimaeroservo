function plotdisprot(iEval,ModalRE_all, ModalCE_all, timeRange, maxZ, maxTheta)
% This function returns plot of displacement and rotation in 2d figure
% INPUT:
% ModalRE_all & 3xNumNodexnumel(timeRange)   & Database of displacement along the time(from time marching code)
% ModalCE_all & 3x3xNumNodexnumel(timeRange) & Database of rotation along the time(from time marching code)
% timeRange   & 1xnumel(timeRange)           & Time range
% maxZ        & 1x1                          & Maximum value of displacement. Used for limiting axis plot
% maxTheta    & 1x1                          & Maximum value of rotation. Used for limiting axis plot

% Initialise video
writerObj = VideoWriter(sprintf('pole%d.avi',iEval));
open(writerObj);

% Video 
for frame = 1:numel(timeRange)
    xDisp = ModalRE_all(1,:,frame);         
    zDisp = ModalRE_all(3,:,frame);
    theta = atan(ModalCE_all(3,2,:,frame)./sqrt(ModalCE_all(2,2,:,frame).^2+ModalCE_all(1,2,:,frame).^2));
    theta = theta(:);
    h = figure;
    set(h,'Visible','off');
    subplot(2,1,1);
    plot(xDisp,zDisp);
    xlabel('x'); ylabel('z disp');  
    xlim([0 7])
    ylim([-maxZ maxZ])

    subplot(2,1,2);
    plot(xDisp,theta);
    xlabel('x'); ylabel('\theta');
    xlim([0 7])
    ylim([-maxTheta maxTheta])

    mov(frame) = getframe(gcf);
    writeVideo(writerObj,mov(frame));
end
close(writerObj);