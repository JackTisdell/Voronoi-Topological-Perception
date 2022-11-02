function [A, verts] = voronoiForwardArea(DT,U)
arguments
    DT delaunayTriangulation
    U (:,2) double {mustBeCompatible(U,DT)}
end

N = size(U,1);
verts = cell(N,1);
A = Inf(N,1);
Ul = ([0 -1; 1 0]*U')';
Ur = ([0 1; -1 0]*U')';
Ql = voronoiProjectToBoundary(DT,Ul);
Qr = voronoiProjectToBoundary(DT,Ur);
infArea = (Ql(:,1)==Inf)|(Qr(:,1)==Inf);

X = DT.Points;
[V,C] = voronoiDiagram(DT);
al = atan2(Ql(:,2)-X(:,2), Ql(:,1)-X(:,1));
ar = atan2(Qr(:,2)-X(:,2), Qr(:,1)-X(:,1));

CH = convexHull(DT);
IN = setdiff(1:N, CH);
[~,X_bd] = freeBoundary(DT);
d = circshift(X_bd,-1,1) - X_bd;
unbd_edge_arg = atan2(d(:,2),d(:,1))-pi/2;
unbd_edge_arg = unbd_edge_arg + 2*pi*(unbd_edge_arg < -pi);

for l = 1:size(IN,2)
    i = IN(l);
    v = V(C{i},:);
    numverts = size(v,1);
    a = atan2( v(:,2)-X(i,2), v(:,1)-X(i,1) );
    [~,S] = sort([a' al(i) ar(i)]);
    S = circshift(S, 1-find(S==numverts+2,1));
    k = find(S==numverts+1,1);
    P = [Qr(i,:); v(S(2:k-1)',:); Ql(i,:)];
    verts{i} = P;
    A(i) = poly_area(P(:,1), P(:,2));
end

m = size(CH,1)-1;
for l=1:m
    i = CH(l);                      % idx of generator on convex hull
    if infArea(i) == true
        continue;
    end
    v = V(C{i},:);                  % vertices of i's voronoi cell
    numverts = size(v,1)-1;         % number of (finite) vertices of v
    a = zeros(numverts+2,1);        % initialize arrary of angles in [-pi,pi] w.r.t. generator
    a(1) = unbd_edge_arg(l);        % angle of clockwise-out unbouded edge
    a(2:end-1) = atan2( v(2:end,2)-X(i,2), v(2:end,1)-X(i,1) ); % angles to bounded vertices in counterclockwise order
    a(end) = unbd_edge_arg(l-1 + m*(l==1));      % angle of counterclockwise-out unbounded edge
    [~,S] = sort([a' al(i) ar(i)]); % append 9 and 3 o'clock angles w.r.t. U(i)'s 12 o'clock and sort
    S = circshift(S, 1-find(S==numverts+4,1));  % permute cyclically to put 3 o'clock in 1st position
    k = find(S == numverts+3, 1);   % locate sorted position of 9 o'clock
    kin = find(S==1, 1);            % -----------//------------ clockwise-out unbounded edge direction
    kout = find(S==numverts+2,1);   % -----------//------------ counterclockwise-out unbounded edge direction
    if kin < k || kout < k          % if either unbounded direction is forward-visible, continue to next iteration
        continue;
    end
    P = [Qr(i,:); v(S(2:k-1)',:); Ql(i,:)];     % else, find polygon
    verts{i} = P;
    A(i) = poly_area(P(:,1), P(:,2));
end

end

% validation
function mustBeCompatible(U,DT)
    if size(U,1) ~= size(DT.Points,1)
        eid = 'Size:notEqual';
        msg = 'Number of rows of second input must match the number of points in the triangulation.';
        throwAsCaller(MException(eid,msg));
    end
end