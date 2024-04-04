tmax = size(U_t,3);

interval = 1:tmax;

clip = [0 0 0 0];
for t=interval
    DT = DT_t{t};
    X = DT.Points;
    clip(1) = min(clip(1),min(X(:,1)));
    clip(2) = max(clip(2),max(X(:,1)));
    clip(3) = min(clip(3),min(X(:,2)));
    clip(4) = max(clip(4),max(X(:,2)));
end


for t=interval
    DT = DT_t{t};
    X = DT.Points;
    com = mean(X);
    rmed = median(vecnorm(X-com,2,2));
    U = U_t(:,:,t);

    fig = figure(1);
    quiver(X(:,1),X(:,2),U(:,1),U(:,2));
%     quiver(X(1:N/2,1),X(1:N/2,2),U(1:N/2,1),U(1:N/2,2),'r');
%     hold on
%     quiver(X(N/2+1:end,1),X(N/2+1:end,2),U(N/2+1:end,1),U(N/2+1:end,2),'b');
    hold on
    scatter(0,0,15,'r','filled');
%     scatter([-5 5]',[0 0]',15,'r','filled');
    rectangle('Position',[com-rmed 2*rmed 2*rmed],'Curvature',[1 1],'EdgeColor','g');
    scatter(com(1),com(2),25,'k','+');   % CoM crosshairs
    hold off
    axis equal
    axis(clip);
    axis([-10 10 -10 10]);
    title(['$t=',num2str(t),'$'], 'Interpreter', 'Latex');
end