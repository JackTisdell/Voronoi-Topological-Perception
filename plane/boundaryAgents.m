com = mean(X);
[nbhd, nearest, d] = neighborhoods(DT);

CH = convexHull(DT);
CH = CH(1:end-1);
CH2nearest = X(nearest(CH),:) - X(CH,:);
CH2com = com - X(CH,:);
radialDeviation = acos(dot(normalize(CH2com,2,'norm'),normalize(CH2nearest,2,'norm'),2));

quiver(X(:,1),X(:,2),U(:,1),U(:,2));
hold on
% text(X(:,1),X(:,2),num2cell(1:N));
quiver(X(CH,1),X(CH,2),CH2nearest(:,1),CH2nearest(:,2),0,'r');

axis equal

hold off
