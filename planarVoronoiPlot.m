%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the planar Voronoi diagram for a given Delaunay triangulation
% including (initial segments of) its unbounded edges.
% 
% Input:
%   DT  :   delaunayTriangulation object
%   E   :   OPTIONAL argmuent, Nx3 double where N is the number of rows in
%           DT.Points. E specifies the unbouded edges as follows. The i-th 
%           row of E represents the unbounded edge e which comes in 
%           counterclockwise from infinity of the Voronoi cell whose 
%           generator is the i-th element of convexHull(DT). Specifically, 
%           E(i,1) and E(i,2) are the x and y coordinates of the (finite) 
%           vertex of e and E(i,3) is the angle with respect to the
%           positive x axis at which e radiates to infinity.
%
% Notes:
%   - This routine assumes the values in E are correct but will run without
%     issue for any numerical values of E. The possibility to specificy E
%     as an input argument is intended only to dovetail planarVoronoiPlot
%     with voronoiEdgeIntersection, which computes E out of necessity. When
%     used together, providing E as an input here saves its redundant
%     computation.

function fig = planarVoronoiPlot(DT,E)

X = DT.Points;
[V,C] = voronoiDiagram(DT);

if nargin < 2
    [~,X_bd] = freeBoundary(DT);
    d = circshift(X_bd,-1,1) - X_bd;
    unbd_edge_arg = mod(atan2(d(:,2),d(:,1))-pi/2, 2*pi);
    
    CH = convexHull(DT);
    m = size(CH,1)-1;
    E = zeros(m,3);

    for l = 1:m
        c = V(C{CH(l)},:);  
        E(l,:) = [ c(2,:) unbd_edge_arg(l) ];
    end
end


c = 'blue';


% [minX minY maxX maxY] among all generators and (bounded) Voronoi vertices
ax_lim = [min([V; X],[],1) max([V(2:end,:); X],[],1)];

margin=.1;

draw_dist = 2*(ax_lim(3)-ax_lim(1) + ax_lim(4)-ax_lim(2));

fig = figure(1);

hold on

for i = 1:size(X,1)
    patch(V(C{i},1), V(C{i},2), c, 'FaceAlpha', 0);   
end
quiver(E(:,1),E(:,2), draw_dist*cos(E(:,3)), draw_dist*sin(E(:,3)), 'off', 'k', 'ShowArrowHead', 0);

hold off

axis equal
axis([ax_lim(1)-margin ax_lim(3)+margin ax_lim(2)-margin ax_lim(4)+margin]);

end