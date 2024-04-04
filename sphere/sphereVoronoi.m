function [Vertices, K,T] = sphereVoronoi(xyz)
%
% Compute the voronoi's diagram of points on the spheres S^2
%
% INPUT:
%   xyz is (3 x n) array, coordinates of n points in R^3
%   Requirement: they all must be be normalized to 1 (i.e., belong to the
%   2-sphere) and distincts
%
% OUTPUTS:
%   - Vertices is (3 x m) array, coordinates of the vertices of voronoi diagram
%   - K is (n x 1) cell, each K{j} contains the indices of the voronoi cell
%       vertices correspond to xyz(:,j).
%       Vertices are counter-clockwise oriented when looking from outside.
%   - T is the triangulation of the convex hull ---> which is the Delaunay
%   triangulation of the surface of the sphere. T can be used for example
%   to compute the topological distances between Voronoi generato
%
%   
% All credits to the original Author: Bruno Luong <brunoluong@yahoo.com>
% Script Obtained from the MATLAB File Exchange Forum 
% link -->  https://www.mathworks.com/matlabcentral/fileexchange/40989-voronoi-sphere
%
% Date of creation by Author: 28/March/2013 
%
% Last known update by Author: 17/July/2018:
%                
%
% Modified by Ivan Gonzalez @McGill University
% Last update: April 14th, 2019
%%
npnts = size(xyz,2);

T = convhull(xyz'); %triangulation of the convexhull surface
nt = size(T,1); %number of triangles in T

E = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])]; %Concatenated idices of vertices that form each edge of thetriangulation
E = sort(E,2); %the left endpoint of each edge needs to have a smaller index

% Now we check that each triangulation edge is not repeated (if it is there is a topological problem)
[~, ~, J] = unique(E, 'rows');
if ~all(accumarray(J,1)==2) %we count the nuumber of times that each edge is represented...we need to have accumarray(J,1)==2 for each edge to be OK
    error('Topology issue due to numerical precision')
end

% Which 2 seeds the edge vertices correspond?
allids = repmat((1:nt).',[3 1]);
k = accumarray(J, (1:3*nt).', [], @(k) {k});
k = [k{:}];
vid = allids(k.');

% each row is 2 cell ids of the edge
cellofedge = E(k(1,:),:); % ne x 2
ne = size(cellofedge,1);
edges = repmat((1:ne).',[2 1]);
edgeofcell = accumarray(cellofedge(:),edges, [], @(e) {e});

% Center of the circumscribed Delaunay triangles
Vertices = Center(xyz, T);

% Build the contour of the voronoi cells
K = cell(size(edgeofcell));
for k = 1:npnts
    % ordering and orientation of the edges
    v = cycling_edge(edgeofcell{k}, vid);
    v = oriented_edge(v, Vertices, xyz(:,k));
    K{k} = v(:,1);
end
Vertices=Vertices';
end % voronoisphere


%%
function P = Center(xyz, T)
% Center of the circumscribed Delaunay triangles
XYZ = reshape(xyz(:,T),[3 size(T)]);
A = XYZ(:,:,1);
B = XYZ(:,:,2);
C = XYZ(:,:,3);
A = A-C;
B = B-C;
A2B = bsxfun(@times, sum(A.^2,1), B);
B2A = bsxfun(@times, sum(B.^2,1), A);
AxB = cross(A,B,1);
P = cross(A2B - B2A, AxB, 1);
P = C + bsxfun(@times,P,1./(2*sum(AxB.^2,1)));
nP = sqrt(sum(P.^2,1));
P = bsxfun(@times, P, 1./nP);
s = dot(AxB,C);
P = bsxfun(@times, P, sign(s));
end % Center

%%
function v = cycling_edge(edges, vertexes)
% Chain the edges in cycle
u = vertexes(edges,:).';
n = size(u, 2);
[~, ~, I] = unique(u);
I = reshape(I,[2 n]);
J = repmat(1:n,[2 1]);
if ~all(accumarray(I(:), 1) == 2)
    error('Topology issue due to numerical precision')
end
K = accumarray(I(:), J(:), [], @(x) {x});
K = [K{:}];
v = zeros([n 2]);
p = 0;
q = 1;
% chain the edges
for j = 1:n
    i = K(:,q);
    if i(1) == p
        p = i(2);
    else
        p = i(1);
    end    
    i = I(:,p);
    if i(1) == q
        v(j,:) = u([1 2],p);
        q = i(2);
    else
        v(j,:) = u([2 1],p);
        q = i(1);
    end
end % for-loop
end % cycling_edge

%%
function v = oriented_edge(v, P, xyz)
% Orient the edges counter-clockwise
Q = null(xyz.');
E = P(:,v([1:end 1],1));
xy = Q'*E;
a = (xy(1,1:end-1)-xy(1,2:end))*(xy(2,1:end-1)+xy(2,2:end))';
if xor(a < 0, det([xyz Q]) < 0) % Combine orientation and directness of Q
    v = rot90(v,2);
end
end % oriented_edge
