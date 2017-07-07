function plot3Dwing(iEval, ModalRE_all, ModalCE_all, aebvalue, timeRange, Nseg, BeamSeg, maxZ)
% This function returns plot of displacement and rotation in 2d figure
% INPUT:
% ModalRE_all & 3xNumNodexnumel(timeRange)   & Database of displacement along the time(from time marching code)
% ModalCE_all & 3x3xNumNodexnumel(timeRange) & Database of rotation along the time(from time marching code)
% aebvalue    & 1x1                          & Half chord length
% timeRange   & 1xnumel(timeRange)           & Time range
% Nseg        & 1x1                          & Number of segment. same with number of beam
% BeamSeg     & struct                       & Datastructure of BeamSeg
% maxZ        & 1x1                          & Maximum value of displacement. Used for limiting axis plot


% Initialise video
writerObj = VideoWriter(sprintf('3Dpole%d.avi',iEval));
open(writerObj);

% Leading edge vector
vecLE = [0;aebvalue;0];

% Trailing edge vector
vecTE = [0;-aebvalue;0];

% figure
for nt = 1
    MR = ModalRE_all(:,:,nt);
    MC = ModalCE_all(:,:,:,nt);
    
    for i = 1:Nseg
        for j = 1:BeamSeg(i).NodeNo-1
            curnode = BeamSeg(i).NodeOrder(j);
            nxtnode = BeamSeg(i).NodeOrder(j+1);
            Rcur = MR(:,curnode);
            Rnxt = MR(:,nxtnode);
            Ccur = MC(:,:,curnode);
            Cnxt = MC(:,:,nxtnode);
            vLEcur = Ccur*vecLE;
            vTEcur = Ccur*vecTE;
            vLEnxt = Cnxt*vecLE;
            vTEnxt = Cnxt*vecTE;
            n1 = Rcur+vLEcur;
            n2 = Rcur+vTEcur;
            n3 = Rnxt+vTEnxt;
            n4 = Rnxt+vLEnxt;
            v1 = (n1-n4)+(n2-n3);
            v2 = (n1-n2)+(n4-n3);
            normdir = cross(v2,v1);
            normdir = normdir/norm(normdir);
            Nmat = [n1,n2,n3,n4];
            vX = Nmat(1,:)';
            vY = Nmat(2,:)';
            vZ = Nmat(3,:)';
            zlim([-maxZ+0.1*maxZ maxZ+0.1*maxZ])
            fill3(vX,vY,vZ,[0,0,1]);
            hold on;
        end
    end
    axis equal 
    zlim([-maxZ+0.1*maxZ maxZ+0.1*maxZ]);
    drawnow
    mov(nt) = getframe(gcf);
    writeVideo(writerObj,mov(nt));
%     clf;
end
