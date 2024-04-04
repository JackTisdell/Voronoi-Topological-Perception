%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEIGHBORHOODS
% Get graph neighborhoods of vertices of Delaunay Triangulation
%
% Syntax:
%   [nbhd,nearest,dist] = neighborhoods(DT)
%   [nbhd,nearest,dist] = neighborhoods(DT,r)
% Input:
%   DT  :   delaunayTriangulation object
%   r   :   positive integer or Inf. topological radius of each
%           neighborhood. Default 1.
% Output:
%   nbhd:   Nx1 cell array where N is the number of points in DT. nbhd{i}
%           is a row vector containing the indices of points in DT
%           reachable from point i via a walk in DT of length at most r
%           excluding i itself.
%   nearest: Nx1 double. nearest(i) is the index of the vertex in DT
%           distinct from i which is closest to i.
%   dist:   Nx1 double. dist(i) is the distance between vertex i and its
%           nearest neighbor in DT.


function [nbhd,nearest,dist] = neighborhoods(DT,r)

arguments
    DT delaunayTriangulation
    r  (1,1) double {mustBePositiveIntegerOrInf} = 1
end

X = DT.Points;
N = size(X,1);
if r==Inf, r=N; end
E = edges(DT);
EE = [ E ; E(:,[2 1]) ];
v = ones(size(EE,1),1);
adj = sparse(EE(:,1),EE(:,2),v);
nbhd = cell(N,r);
nearest = zeros(N,1);
dist = zeros(N,1);
for i=1:N
    neighbors = find(adj(i,:)==1);
    [d,l] = min(arrayfun(@(j) norm(X(i,:)-X(j,:)), neighbors));
    dist(i) = d;
    nearest(i) = neighbors(l);
    nbhd{i,1} = neighbors;
end

neye = ones(N,N)-eye(N);  % false diagonal, true off-diagonal

if r > 1
    powAdj = eye(N);
    cumAdj = zeros(N,N);
    for k=1:r
        powAdj = powAdj*adj;
        cumAdj = cumAdj+powAdj;
        cumAdj = cumAdj & neye;
        for i=1:N
            nbhd{i,k} = find(cumAdj(i,:)~=0);
        end
        cumAdj = cumAdj | eye(N);
        if all(cumAdj,'all')
            break;
        end
    end
end

end

% validation
function mustBePositiveIntegerOrInf(a)
    if ~((a==floor(a) && a>0) || a==Inf)
        eid = 'neighborhoods:inputError';
        msg = 'Radius must be a positive integer or Inf';
        throwAsCaller(MException(eid,msg));
    end
end