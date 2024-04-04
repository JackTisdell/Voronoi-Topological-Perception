function [CMx,CMy,CMz,A,Energy,E]=poly_centroid_area_scaledEnergy_sphere(x,y,z,V,C,groundState)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the centroid, the area and energy of Voronoi Diagram on the Sphere S2.
% We compute area with Girard's formula based on internal angles and centroid using a closed form formula obtained by
% applying Stokes Theorem to tansform the surface inegrals over the Voronoi
% regions V_i into line integrals over the boundaries of the polygons.
%
%
% Inputs: - The topology of the Voronoi diagram is specified by the list of
%         vertices V and the Conectivity list C. The vertices should be in
%         counterclock wise order when seen from outside the sphere (a
%         small block of code performs that verification but has been 
%         commented since the function sphereVoronoi.m gives vertices in
%         good order.
%
% Ouptuts: - A is the area of each cell computed with the closed formula
%            you can see its accuracy by computing sum(A)-4*pi (order ~1e-12 for N=1000)
%          
%          - CMx,CMy,CMz are the coordinates of the centroid of the
%            Voronoi cell (on the sphere).
%
%          - E is an Nx1 vector containing the idividual scaled Voronoi cell
%          energies. E(i) is calculated using the dot product of the
%          x_i with the Euclian centroid (before projection onto the S2)
%
%          - Energy is the total scaled energy of the configuration
%
% Writen by Ivan Gonzalez @McGill University
% Writen on April 17th, 2020
% Last update: April 17th, 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Memory Allocation for variables 
N=length(C);
A=zeros(N,1); %Area of V_i obtained by our closed form formula (no need to use the triangulation of V_i)
CMx=A;CMy=A;CMz=A; % the centroid obtained with out closed form formula
%%
for i=1:N
   
    v1=V(C{i},1);v2=V(C{i},2);v3=V(C{i},3); % list of vertices of V_i
    Nneigh=length(v1); % Number of voronoi neighbors (same as the number of internal angles of the polygon)
    vv1=[v1(end);v1;v1(1)];vv2=[v2(end);v2;v2(1)];vv3=[v3(end);v3;v3(1)]; % Doubly Closed (circulatory) List of Voronoi Vertices so that we can access both indices j-1 and j+1 in the loop over j down below
    
    %% closed formula written in terms of cross product (use this one rather than the dot product one)
    % It appears that this closed form formula is more accurate than the
    % one expanded in terms of dot produc of vertices. We note that its
    % accuracy w.r.t 4pi (the surface of the sphere) is comparable to the
    % method consisting of triangulating every V_i and using L'Huillier's
    % formula. We conclude that this one is the one to go with.
    alph=zeros(Nneigh,1); % allocation for internal angles of the Voronoi polygon
    c1=alph; % x-coordinate of contribution to the euclidan centroid by each edge joining v_j to v_j+1 
    c2=alph; % y-coordinate of contribution to the euclidan centroid by each edge joining v_j to v_j+1 
    c3=alph; % z-coordinate of contribution to the euclidan centroid by each edge joining v_j to v_j+1 
    
    
    for j=2:Nneigh+1
        %% auxilliary quantities
        vj=[vv1(j);vv2(j);vv3(j)]; % current Voronoi vertex, we declare this variable because it will be used a lot
        a=cross(vj,[vv1(j-1);vv2(j-1);vv3(j-1)]); %normal vector to the plane passing by the origin and vertices j and j-1
        b=cross(vj,[vv1(j+1);vv2(j+1);vv3(j+1)]); %normal vector to the plane passing by the origin and vertices j and j+1
        na=sqrt(a'*a); %norm of the normal vector a
        nb=sqrt(b'*b); %norm of the normal vector b
        c=cross(b,vj); % (we called this vector u_j in our log book) it's the triple cross poduct [(v_j)x(v_j+1)]x(v_j) giving orthogonal vector to vj that lies on the plane generated by the origin, v_j and v_j+1
        % note that nc=nb because we are on the sphere
        c=c/nb; %normalize vector c, (in the log book we called it w_j, which is the normalized version of u_j)
        
        d=cross(vj,c); % fourth cross produc (v_j)x[(v_j)x(v_j+1)]x(v_j) appearing in the line integrals for the centroid
        ee=vj'*[vv1(j+1);vv2(j+1);vv3(j+1)]; % dot product between two consecutive vertices v_j and v_j+1
        
        %% Test the counterclockwise orientation of the vertcies when looked from above;
        % the paralellepiped generated by vj, v_j+1 and x_i needs to be
        % generated by "un tri�dre positif"
        % NOTE: So far all the configurations passed the test
        
%         if ([x(i) y(i) z(i)]*b)<0
%             warning('Polygon with Clockwise or undefined orientation detected')
%         end

        %% Computation of the internal angles of the Voronoi polygon (needed for Girard's formula of area)
        alph(j-1)=arc_cosine((a'*b)/(na*nb));
        
        %% Centroid calculations
        
        % for k=1 component (x-component)
        c1(j-1)=0.5*( arc_cosine(ee)*d(1) + nb^2*(c(2)*c(3)-vj(2)*vj(3)) +ee*nb*(vj(2)*c(3)+vj(3)*c(2)) );
        
        % for k=2 component (y-component)
        c2(j-1)=0.5*( arc_cosine(ee)*d(2) + nb^2*(c(1)*c(3)-vj(1)*vj(3)) +ee*nb*(vj(1)*c(3)+vj(3)*c(1)) );
        
        % for k=3 component (z-component)
        c3(j-1)=0.5*( arc_cosine(ee)*d(3) + nb^2*(c(1)*c(2)-vj(1)*vj(2)) +ee*nb*(vj(1)*c(2)+vj(2)*c(1)) );
        
    end
    A(i)=sum(alph)+(2-Nneigh)*pi; % Girard's formula spherical excess formula in terms of the internal angles and the number of sides/angles of the polygon (https://en.wikipedia.org/wiki/Spherical_trigonometry)
    CMx(i)=sum(c1)/A(i);CMy(i)=sum(c2)/A(i);CMz(i)=sum(c3)/A(i); %This is the euclidian centroid of the spherical sruface polygon that is V_i

end

%% Energy computed using the euclidian centroid
CM=[CMx CMy CMz]; %intermediary
X=[x y z];

E=2*A.*(1-sum(X.*CM,2)); %individual cell energies (NOT scaled)
Energy=sum(E); % total energy of the tessellation (NOT scaled)

E=(N*E)/groundState; % individual energies scaled by the average energy per cell in the soccer ball ground state
Energy=Energy/groundState; %total scaled energy of the configuration.
%% Normalization of the euclidian centroid to project onto the sphere
nor=sqrt(sum(CM.^2,2)); % euclidian norm of the centroids
CMx=CMx./nor;CMy=CMy./nor;CMz=CMz./nor; %constrained centroids projected onto the sphere from the euclidian centroids. 

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