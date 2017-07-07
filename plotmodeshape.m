function plotmodeshape(BeamSeg,Vreal,iMode)
plot(BeamSeg.NodeL,reshape(Vreal(:,iMode),6,21));
saveas(gcf,sprintf('%d.jpg',iMode));
saveas(gcf,sprintf('%d.fig',iMode));
xlabel('X');