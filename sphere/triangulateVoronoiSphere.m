function [vertices_tri,tri,inds_tri]=triangulateVoronoiSphere(vertices,faces,inds)
% triangulateVoronoiSphere triangulates a spherical Voronoi tessellation to
% obtain a smoother meshing. It is able to perform subtriangulations by
% calling the function several times, each time with the data from the
% previous call.
%
%   INPUTS: 
%   - "vertices" is a Kx3 array of the list of vertices making up the
%     polyhedron (K is the total number of unrepeated vertices)
%   - "faces" is EITHER cell array or a NUMERIC array such that faces{i}
%     or faces(i,:) contain the index of vertices forming face i.
%     NOTE 1: "faces" should be a cell array when faces have different number
%     of vertices. "faces" should be a numeric array when all faces have
%     the same number of vertices (e.g. in a triangulation).
%     NOTE 2: The faces NEED to be CONVEX.
%   - "inds" contains the index identification number of each "face" w.r.t the
%    parent faces, e.g. if inds(i)=7 this means that face i is obtain by 
%    subdivision of the original face 7. This is useful to plot all faces 
%    obtained by subdivision with the same color as the original face from 
%    which they were computed.
%
%   OUTPUTS
%   - "vertices_tri" is the list of vertices of the subdived triangulation.
%     it contains "vertices" as well as the mean/average point of the 
%     vertices of each "face" (projected onto the sphere).
%   - "tri" is the connectivity list of the new triangulation.
%   - "inds_tri" contains the identification number of the newly subdivided
%    triangle w.r.t its parent face in "inds".
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inspired by the function TriangulateFaces.m of the "geom3d" toolbox from 
% https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d
% Author: David Legland
% Created: 2008-09-08,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified from its original version by Ivan Gonzalez @McGill University
% Date: April 20th 2021
% Last update: April 20th 2021
%
%
%%
nv=size(vertices,1);
% number of original faces 
nf=length(faces);

% number of vertices in the new triangulation (original number of vertices
% plus the mean vertex of each face)
nv_tri=nv+nf;

% compute number of triangles at output
if isnumeric(faces) %for 2nd or subsequent layer of triangulation each face is already a triangle and it will be split into THREE triangles
    ni=3*ones(nf,1);
else
    ni=cellfun(@length,faces); % number of triangles per face for the first later (i.e. for not necessarely triangular faces).
end
nt = sum(ni); % total number of triangles in the new triangulation;

% To Compute and include new vertices made of the mean of each face's vertices
vertices_tri=zeros(nv_tri,3);
vertices_tri(1:nv,:)=vertices;

% allocate memory for triangle array
tri = zeros(nt, 3);
inds_tri = zeros(nt, 1);

% Compute mean vertices and triangulate faces.
t = 1;
for i = 1:nf
    if isnumeric(faces)
        face=faces(i,:)'; %face in the second and subsequent layers where the faces are in NUMERIC form
    else
        face=faces{i}; % face in the first layer where the faces are in CELL array form
    end
    
    meanvert=mean([vertices(face,1) vertices(face,2) vertices(face,3)]); % mean value of face's vertices
    meanvert=meanvert/norm(meanvert,2); %poject the mean point onto the surface of the sphere.
    vertices_tri(nv+i,:)=meanvert; % include new vertes in the list of vertices of the new triangulation.
    
    face_shift=[face(2:end); face(1)]; %shifted list of face vertices
    tri(t:(t+ni(i)-1),:)=[(nv+i)*ones(ni(i),1) face face_shift]; % insert new triangulation of parent face i. The shared vertex is the mean point of the parent vertices.
                                                                 % We insert each of the new ni(i) triangles of parent face i at the same time
    inds_tri(t:(t+ni(i)-1))=inds(i);
    t=t+ni(i);
end

end

