function [A,A0,T] = inwardTotalArea(DT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a Delaunay Triangulation object DT, this routine
% computes for each vertex p, the area of the (unique) neigboring 
% triangle in the triangulation whose projection from p contains the 
% center of mass of the vertices. 
%
% Output: 
%   A   : total area of the triangles described above
%   A0  : area of the (unique) triangle containing the center of mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = DT.Points;
N = size(P, 1);
cx = mean(P(:,1));
cy = mean(P(:,2));

com_triangle = pointLocation(DT,cx,cy); % triangle containg center of mass
com_triangle_verts = DT.ConnectivityList(com_triangle,:);
com_triangle_coords = DT.Points(com_triangle_verts',:);
A0 = 1/2*abs(det( [com_triangle_coords ones(3,1)] ));

angles2com = atan2(cy-P(:,2), cx-P(:,1));

VA = vertexAttachments(DT);
T = zeros(3,2,N);
A = 0;
for i = 1:N
    tri_idx = VA{i}';           % get IDs of m triangles containing i
    tris = DT.ConnectivityList(tri_idx,:);  % get vertices of those triangles (m x 3 array)
    verts = unique(tris);       % list those vertices (without repeats)
    verts = verts(verts~=i);    % exclude i itself
    pts = P(verts,:);   % get the coords of those vertices
    angs = atan2(pts(:,2)-P(i,2),pts(:,1)-P(i,1));  % get their angle from i
    [~,I] = sort( [angles2com(i); angs] );    % sort these verts AND the COM by angle
    pts_com = [ cx cy ; pts ];
    pts_com = pts_com(I,:);
    k = find(I==1);             % get position of COM after sorting
    pts_com = circshift(pts_com, 2-k, 1);

    T(:,:,i) = [ P(i,:) ; pts_com(1,:) ; pts_com(3,:) ];  
    A = A + 1/2*abs(det( [ P(i,:) 1 ; pts_com(1,:) 1 ; pts_com(3,:) 1 ] ));
end

A = A - 2*A0;

% hold on
% triplot(DT);
% for i = 1:N
%     fill(T(:,1,i), T(:,2,i), 'b', 'FaceAlpha', '0.2');
% end
% scatter([cx], [cy],'red');
% text(x,y,num2cell(1:N));
end