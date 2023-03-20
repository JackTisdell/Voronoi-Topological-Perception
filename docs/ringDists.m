function [rmed,rra] = ringDists(DT,c) 

arguments
    DT delaunayTriangulation
    c (1,2) double = [NaN NaN];
end

X = DT.Points;
if isnan(c(1))
    c = mean(X);
end

rmed = median(vecnorm(X-c,2,2));

tri = pointLocation(DT,c);
verts = DT.ConnectivityList(tri,:);
coords = X(verts,:);
rra = sqrt(poly_area(coords(:,1),coords(:,2)));

end
