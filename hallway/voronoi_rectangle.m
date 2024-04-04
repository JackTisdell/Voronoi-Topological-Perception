function [DT0,DT,V,C,idx_bdry]=voronoi_rectangle(x,y,Omega)

% This routine computes the Delaunay triangulation and the Voronoi topology
% of a set of generators on the rectangle Omega.
% In order to obtain a bouded Voronoi diagram we insert ghost copies of
% the generators. The generators that are copied are those that have a
% Voronoi vertex outside Omega...equivalently those that are attached to a
% Delaunay triangle whose circumcenter is outside Omega.


Lx=Omega(2,1); Ly=Omega(3,2);

DT0=delaunayTriangulation(x,y);

DT=DT0;

% Using Lmax "guarantees" that if the aspect ratio of the domain is two big
% then the Voronoi regions of the 4 ghost points at "infinity" do not
% intersect Omega.
Lmax=max([Lx Ly]);
DT.Points(end+(1:4),:)=[[Lmax/2; 6*Lmax; Lmax/2; -5*Lmax] [-5*Lmax; Lmax/2; 6*Lmax; Lmax/2]];

conList=DT.ConnectivityList; % we could also have used DT(:,:)

% We find the points that will be labeled as "boundary" points.
CC=circumcenter(DT); %compute circumcenters
out=~inpolygon(CC(:,1),CC(:,2),Omega(:,1),Omega(:,2)); % tell if the CC (circumcenter) is outside the primary domain Omega
idx_bdry=conList(out,:); %give the index of the generators who have the CC outside D0 (the output is a list of triangles from DT)
idx_bdry=unique(idx_bdry(:)); % concatenate the indices in a column vector, remove duplicates and sort in increasing order
idx_bdry=idx_bdry(1:end-4); % remove the four extra points we used to have bounded Voronoi regions (these are the last four since idx_bdry is in increasing order).

x3=x(idx_bdry);y3=y(idx_bdry); % the points to reflect and copy

x_ref_copy=[x3;2*Lx-x3;x3;-x3]; %the reflected copies on each of the four edge-adjacent domains.
y_ref_copy=[-y3;y3;2*Ly-y3;y3];

pp=4*length(idx_bdry); % number of reflected copies that were maid (p==length(x_ref_copy))

% We Incrementally update the Delaunay triangulation and Voronoi Tesselation
DT.Points(end+(1:pp),:)=[x_ref_copy y_ref_copy];
[V,C]=voronoiDiagram(DT); % [V,C] gives the topology of the Voronoi
% diagram of the extended set of generators
end