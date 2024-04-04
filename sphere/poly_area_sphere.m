function A=poly_area_sphere(Vertices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compute area with Girard's formula based on internal angles and centroid using a closed form formula obtained by
% applying Stokes Theorem to tansform the surface inegrals over the polygon
% into line integrals over the boundaries of the polygons.
%
%
% Inputs: - Vertices is the list of vertex coordinates of a SINGLE
%           spherical polygon...NOTE: WE DO NOT COMPUTE THE AREA OF EACH
%           ELEMENT IN A TESSELLATION
%
% Ouptuts: - A of the polygonal region.
%          
%
% EXTRACTED FROM THE FUNCTION poly_centroid_area_scaledEnergy_sphere.m AND
% ADAPTED TO COMPUE AREA OF A SINGLE POLYGON. 
%
% Writen by Ivan Gonzalez @McGill University
% Writen on April 16th, 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
v1=Vertices(:,1);v2=Vertices(:,2);v3=Vertices(:,3); % extract coordinates of vertices
Nneigh=length(v1); % Number of edges (same as the number of internal angles of the polygon)
vv1=[v1(end);v1;v1(1)];vv2=[v2(end);v2;v2(1)];vv3=[v3(end);v3;v3(1)]; % Doubly closed (circulatory) List of Vertices so that we can access both indices j-1 and j+1 in the loop over j down below

%% closed formula written in terms of cross product (use this one rather than the dot product one)
% It appears that this closed form formula is more accurate than the
% one expanded in terms of dot produc of vertices. We note that its
% accuracy w.r.t 4pi (the surface of the sphere) is comparable to the
% method consisting of triangulating the poilygon and using L'Huillier's
% formula. We conclude that this one is the one to go with.
alph=zeros(Nneigh,1); % allocation for internal angles of the polygon
for j=2:Nneigh+1
    %% auxilliary quantities
    vj=[vv1(j);vv2(j);vv3(j)]; % current vertex
    a=cross(vj,[vv1(j-1);vv2(j-1);vv3(j-1)]); %normal vector to the plane passing by the origin and vertices j and j-1
    b=cross(vj,[vv1(j+1);vv2(j+1);vv3(j+1)]); %normal vector to the plane passing by the origin and vertices j and j+1
    na=sqrt(a'*a); %norm of the normal vector a
    nb=sqrt(b'*b); %norm of the normal vector b
    %% Computation of the internal angles of the  polygon (needed for Girard's formula of area)
    % NOTE WE USE THE arc_cosine.m function defined below to avoid complex
    % numbers produced by numerical inaccuracies.
    alph(j-1)=arc_cosine((a'*b)/(na*nb)); % 
    
end
A=sum(alph)+(2-Nneigh)*pi; % Girard's formula: spherical excess formula in terms of the internal angles and the number of sides/angles of the polygon (https://en.wikipedia.org/wiki/Spherical_trigonometry)

end


function value = arc_cosine ( cc )

%*******************************************************************************
%
%% ARC_COSINE computes the arc cosine function, with argument truncation--> to avoid complex numbers.

%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%    Author: John Burkardt
%    Modified: 28 January 2005
%
%  Parameters:
%
%    Input, real C, the argument.
%    Output, real VALUE, an angle whose cosine is C.
%
  cc2 = cc;
  cc2 = max ( cc2, -1.0 );
  cc2 = min ( cc2, +1.0 );

  value = acos ( cc2 );

  return
end