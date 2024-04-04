function plotVoroSphere(N,V,C,col)

for i=1:N
    V1=V(C{i},1);V2=V(C{i},2);V3=V(C{i},3);
    fill3(V1,V2,V3,col,'FaceColor','none','EdgeColor','w','FaceAlpha',1);
    if i==1
        hold on
    end
end
%scatter3(x,y,z,'k.','LineWidth',15);

hold off

axis equal
axis off
% grid on
% grid minor
view(180,0)

E_min=0.85;
E_max=1.7;
c=colorbar;
c.Location='southoutside';
c.Label.String='Second Moments $E_i$';
c.Label.Interpreter='latex';
c.Label.FontSize=15;
c.Limits=[E_min E_max];
caxis([E_min E_max])
colormap(gca,'jet')

end