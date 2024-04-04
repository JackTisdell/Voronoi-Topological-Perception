function [DT0,DT,V,C,idx_copy,idx_bdry] = periodicVoronoi(x,y,Omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the topology of a bounded Voronoi Tesselation in
% D0. For this we reflect the boundary generators on the four adjacent
% domains to D0
%
% INPUT: - x,y are the coordinates of the vertices of generators
%        - Omega contains the specifications of the rectangular domain.
%    
% OUTPUT: - DT0 is the Delaunay of the original N generators
%         - DT is the Delaunay of the N + 4 ghosts at "infinity" + ghost copies
%         - V is the list of vertices of all Voronoi cells
%           (including the copies)
%         - C is cell array specifying which generator contains which
%         vertices of V
%
% Written by Ivan Gonzalez @McGill University
% Created: May 23rd of 2018
% Last update: August 24rd of 2020 (quick adaptation to include it in the
%                                   crowd dynamics project)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx=Omega(2,1); Ly=Omega(3,2);
N=length(x);

%% Delaunay
DT0=delaunayTriangulation(x,y);
DT=DT0;

% Using Lmax "guarantees" that if the aspect ratio of the domain is two big
% then the Voronoi regions of the 4 ghost points at "infinity" do not
% intersect Omega.
Lmax=max([Lx Ly]);
DT.Points(end+(1:4),:)=[[Lmax/2; 6*Lmax; Lmax/2; -5*Lmax] [-5*Lmax; Lmax/2; 6*Lmax; Lmax/2]];

conList=DT.ConnectivityList; % we could also have used DT(:,:)

%% We find the points that will be labeled as "boundary" points.
CC=circumcenter(DT); %compute circumcenters
out=~inpolygon(CC(:,1),CC(:,2),Omega(:,1),Omega(:,2)); % tell if the CC (circumcenter) is outside the primary domain D0
idx_bdry=conList(out,:); %give the index of the generators who have the CC outside D0 (the output is a list of triangles from DT)
idx_bdry=unique(idx_bdry(:)); % concatenate the indices in a column vector, remove duplicates and sort in increasing order
idx_bdry=idx_bdry(1:end-4); % remove the four extra points we used to have bounded Voronoi regions (these are the last four since idx_bdry is in increasing order).

% Next we find the immediate neighboors of the boundary points
% This will give all the points we need to copy in each domain (D1,D2,...D8)
tri=vertexAttachments(DT,idx_bdry);
neigh=[];
for i=1:length(idx_bdry)
    NEIGH=DT(tri{i},:);
    neigh=[neigh;NEIGH(:)];
end
idx_copy=unique(neigh);
idx_copy(end-3:end)=[]; %the final indices of points to copy
                        % NOTE: idx_copy includes idx_bdry

%% We copy the ghost points in all the domains(D1,D2,...D8)
x3=x(idx_copy);y3=y(idx_copy);

x_copy=[x3;x3+Lx;x3+Lx;x3+Lx;x3;x3-Lx;x3-Lx;x3-Lx];
y_copy=[y3-Ly;y3-Ly;y3;y3+Ly;y3+Ly;y3+Ly;y3;y3-Ly];

p=8*length(x3); % number of copies that were maid (p==length(x_copy)=length(y_copy) )
%% We Incrementally update the Delaunay triangulation and Voronoi Tesselation
DT.Points(end+(1:p),:)=[x_copy y_copy];
[V,C]=voronoiDiagram(DT); % [V,C] gives the topology of the Voronoi
%C=C(1:N);
% diagram of the extended set of generators

end

