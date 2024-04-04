%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function aims to compute the theoretical Energy ground state of CVTs
% on the unit Sphere
% ---> The optimal configuration is conjectured to be made out
% of 12 congruent and regular spherical hexagons and N-12 congruent and
% regular spherical hexagons (because by Euler's characteristics there
% is always 12 more pentagons than heptagons).
%
% In order to compute the ground energy we impose two constraints to be met
% by the configuration:
% 1) The total area of the hex/pen configuration needs to be equal to 4*pi (measure of the sphere)
% 2) For the configuration to tile the sphere we need that both hex and
% pent have the same length of sides.
%
% NOTE: for more detail on the formulas used refer to my red log book.
%
%
% Writen by Ivan Gonzalez @McGill University
% Writen on April 2019
%
% Important update: April 15th, 2020 (completely re-written on this day,
% in particular the expressions for energy were changed from 
% numerical quadrature formulas to an explicit closed form expression
% obtained using Stokes Theorem)
%
% Last update: April 20th, 2020 (addition of MeanA2ref, MeanIPRref_flat and
% MeanIPRref_sphere to get the normalization factor for all the regularity
% measures).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GroundState]=GroundState_Sphere(N)

tol_sym=1e-60; %tolerance for symbolic accuracy
tol_double=1e-15; %tolerance for double accuracy

digitsOld = digits; %store the default precision built in the system
digits(64) % set the precision of numerical values

%% 1) Find the dimensions of the triangles Th and Tp
syms theta
Pi=sym(pi); %sympolic pi
gamma=sqrt((3+sqrt(sym(5)))/6); %symbolic gamma;

f=theta*(N-12)+10*asin(gamma*sin(theta))-(Pi/3)*(N-2); %non linear equation to solve for half the internal angle of the hexagons

theta_star=vpasolve(f==0,theta,Pi/3); % the solution to the nonlinear equation for theta using vpasolve (thus we have the solution with 64 digits accuracy)
                                      % NOTE: even though f is symbolic we
                                      % can not use the symbolic solver
                                      % solve(.)...if we do it will
                                      % automatically switch to vpasolve(.)
                                      % since there is no symbolic
                                      % solution to such a transcendental
                                      % equation.

%check that theta_star is close to pi/3 (as shown in my log book)
%error_theta=theta_star-vpa(pp/3); %theta_star is always bigger that pi/3 (60 degrees) and converges to that value as N->Inf
%error_theta_deg=error_theta*180/vpa(pp);

phi_star=( N*(vpa(Pi/3)-theta_star)+12*theta_star-2*vpa(Pi/3) )/10; % the half internal angle of the pentagons

rh_star=arc_cosine(cot(theta_star)*vpa(sqrt(sym(3))));
rp_star=arc_cosine(cot(phi_star)*vpa(sqrt(1+2/sqrt(sym(5)))));

%% The length L of the "outmost side" of the triangles (i.e the one not touching the north pole)
% Further precision: L is the side length of the hexagon as well as the
% side length of the pentagon (they have the same value by contruction).
% The perimeter of the hexagon is 6L while the one of the pentagon is 5L

% Value of L using theta_star and the hexagon
L=arc_cosine(0.5+1.5*cot(theta_star)^2);

% Value of L using phi_star and the pentagon (one can check, they give the
% same value).
%L=arc_cosine(cos(2*Pi/5)+(1-cos(2*Pi/5))*cot(phi_star)^2*(1+2/sqrt(sym(5))));

%% For Th (i.e one of the 6 triangles of the regular hexagon)

% A) area with Girard's formula (since we already have the internal angles,
% pi/3, theta_star and theta_star again).
A_Th=vpa(2*theta_star-2*Pi/3);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A-bis) We could use  L'Huillier's formula instead but Girard's is clearly
% simpler
% each vector here represents the vector position of one vertex (for the closed form formula below v1 is the list of x-coordinates, v2 the list of y-coordinates etc...)
% v1=[0;0;1]; %the north pole
% v2=[sin(rh_star);0;cos(rh_star)]; % second vertex
% v3=[cos(vpa(Pi/3))*sin(rh_star);sin(vpa(Pi/3))*sin(rh_star);cos(rh_star)]; % third vertex
% a=rh_star;b=rh_star;c=arc_cosine(v2'*v3);
% S=(a+b+c)/2;
% A_Th=vpa(4*atan(sqrt(tan(0.5*S)*tan(0.5*(S-a))*tan(0.5*(S-b))*tan(0.5*(S-c)))));% area of the triangle Th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% B) Energy of Th using closed form integration over the boundaries

% NOTE: we only need to compute z-coordinate of the euclidian centroid
% since the energy formula implies taking the dot product of the euclidian
% centroid of Th with the generator located at the north pole (so x,y
% components are redundant)

v1=[0;sin(rh_star);cos(vpa(Pi/3))*sin(rh_star)]; % list of x-cordinates of the three vertices
v2=[0; 0; sin(vpa(Pi/3))*sin(rh_star)];% list of y-cordinates of the three vertices
v3=[1; cos(rh_star); cos(rh_star)]; %list of z-cordinates of the three vertices

vv1=[v1;v1(1)];vv2=[v2;v2(1)];vv3=[v3;v3(1)]; % Closed (circulatory) List of triangle vertices so that we can access index j+1 in the loop over j down below

c3=zeros(3,1); % z-coordinate of contribution to the euclidan centroid by each edge joining v_j to v_j+1 

for j=1:3
    % auxilliary quantities
    vj=[vv1(j);vv2(j);vv3(j)]; % current vertex, we declare this variable because it will be used a lot
    %a=cross(vj,[vv1(j-1);vv2(j-1);vv3(j-1)]); %normal vector to the plane passing by the origin and vertices j and j-1
    b=cross(vj,[vv1(j+1);vv2(j+1);vv3(j+1)]); %normal vector to the plane passing by the origin and vertices j and j+1
    %na=sqrt(a'*a); %norm of the normal vector a
    nb=sqrt(b'*b); %norm of the normal vector b
    c=cross(b,vj); % (we called this vector u_j in our log book) it's the triple cross poduct [(v_j)x(v_j+1)]x(v_j) giving orthogonal vector to vj that lies on the plane generated by the origin, v_j and v_j+1
    % note that nc=nb because we are on the sphere
    c=c/nb; %normalize vector c, (in the log book we called it w_j, which is the normalized version of u_j)
    
    d=cross(vj,c); % fourth cross produc (v_j)x[(v_j)x(v_j+1)]x(v_j) appearing in the line integrals for the centroid
    ee=vj'*[vv1(j+1);vv2(j+1);vv3(j+1)]; % dot product between two consecutive vertices v_j and v_j+1
    
    % Centroid calculations for k=3 component (z-component):
    c3(j)=0.5*( arc_cosine(ee)*d(3) + nb^2*(c(1)*c(2)-vj(1)*vj(2)) +ee*nb*(vj(1)*c(2)+vj(2)*c(1)) );  
end

CMz=sum(c3)/A_Th; % this is the z-coordinate of the euclidian centroid of the triangle Th

E_Th=2*A_Th*(1-CMz); %energy of the triangle Th


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B-bis) poor precision with numerical integration to recover the energy
% (obsolete method).
% each vector here represents the vector position of one vertex (for the closed form formula below v1 is the list of x-coordinates, v2 the list of y-coordinates etc...)
% v1=[0;0;1]; %the north pole
% v2=[sin(rh_star);0;cos(rh_star)]; % second vertex
% v3=[cos(vpa(Pi/3))*sin(rh_star);sin(vpa(Pi/3))*sin(rh_star);cos(rh_star)]; % third vertex
% 
% p=@(s,t) (s*v1+t*v2+(1-t-s)*v3);
% g_e=@(s,t) ((v1'*p(s,t))/(p(s,t)'*p(s,t))^2);
% 
% ener=abs(v3'*cross(v1,v2))*((1/40)*(g_e(0,0)+g_e(1,0)+g_e(0,1))+(1/15)*(g_e(0,0.5)+g_e(0.5,0)+g_e(0.5,0.5))+(9/40)*g_e(1/3,1/3));
% Eh=2*(A_Th-ener); %Energy of one triangle of the hexagon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% for Tp (i.e one of the 5 triangles of the regular pentagon)

% A) area with Girard's formula (since we already have the internal angles,
% 2*pi/5, phi_star and phi_star again).
A_Tp=vpa(2*phi_star-3*Pi/5);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A-bis) We could use  L'Huillier's formula instead but Girard's is clearly
% simpler
% each vector here represents the vector position of one vertex (for the closed form formula below v1 is the list of x-coordinates, v2 the list of y-coordinates etc...)
% v1=vpa([0;0;1]); %the north pole
% v2=[sin(rp_star);0;cos(rp_star)]; %second vertex
% v3=[cos(vpa(2*Pi/5))*sin(rp_star);sin(vpa(2*Pi/5))*sin(rp_star);cos(rp_star)]; %thrid vertex
% a=rp_star;b=rp_star;c=arc_cosine(v2'*v3);
% S=(a+b+c)/2;
% A_Tp=vpa(4*atan(sqrt(tan(0.5*S)*tan(0.5*(S-a))*tan(0.5*(S-b))*tan(0.5*(S-c)))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
% B) Energy of Tp using closed form integration over the boundaries
% NOTE: we only need to compute z-coordinate of the euclidian centroid
% since the energy formula implies taking the dot product of the euclidian
% centroid of Tp with the generator located at the north pole (so x,y
% components are redundant)

v1=[0;sin(rp_star);cos(2*Pi/5)*sin(rp_star)]; % list of x-cordinates of the three vertices
v2=[0;0;sin(2*Pi/5)*sin(rp_star)]; % list of y-cordinates of the three vertices
v3=[1;cos(rp_star);cos(rp_star)]; % list of z-cordinates of the three vertices

vv1=[v1;v1(1)];vv2=[v2;v2(1)];vv3=[v3;v3(1)]; % Closed (circulatory) List of triangle vertices so that we can access index j+1 in the loop over j down below

c3=zeros(3,1); % z-coordinate of contribution to the euclidan centroid by each edge joining v_j to v_j+1 

for j=1:3
    % auxilliary quantities
    vj=[vv1(j);vv2(j);vv3(j)]; % current vertex, we declare this variable because it will be used a lot
    %a=cross(vj,[vv1(j-1);vv2(j-1);vv3(j-1)]); %normal vector to the plane passing by the origin and vertices j and j-1
    b=cross(vj,[vv1(j+1);vv2(j+1);vv3(j+1)]); %normal vector to the plane passing by the origin and vertices j and j+1
    %na=sqrt(a'*a); %norm of the normal vector a
    nb=sqrt(b'*b); %norm of the normal vector b
    c=cross(b,vj); % (we called this vector u_j in our log book) it's the triple cross poduct [(v_j)x(v_j+1)]x(v_j) giving orthogonal vector to vj that lies on the plane generated by the origin, v_j and v_j+1
    % note that nc=nb because we are on the sphere
    c=c/nb; %normalize vector c, (in the log book we called it w_j, which is the normalized version of u_j)
    
    d=cross(vj,c); % fourth cross produc (v_j)x[(v_j)x(v_j+1)]x(v_j) appearing in the line integrals for the centroid
    ee=vj'*[vv1(j+1);vv2(j+1);vv3(j+1)]; % dot product between two consecutive vertices v_j and v_j+1
    
    % Centroid calculations for k=3 component (z-component):
    c3(j)=0.5*( arc_cosine(ee)*d(3) + nb^2*(c(1)*c(2)-vj(1)*vj(2)) +ee*nb*(vj(1)*c(2)+vj(2)*c(1)) );  
end

CMz=sum(c3)/A_Tp; % this is the z-coordinate of the euclidian centroid of the triangle Th

E_Tp=2*A_Tp*(1-CMz); %energy of the triangle Tp

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B-bis) poor precision with numerical integration to recover the energy
% (obsolete method).
% each vector here represents the vector position of one vertex (for the closed form formula below v1 is the list of x-coordinates, v2 the list of y-coordinates etc...)
% v1=vpa([0;0;1]); %the north pole
% v2=[sin(rp_star);0;cos(rp_star)]; %second vertex
% v3=[cos(vpa(2*Pi/5))*sin(rp_star);sin(vpa(2*Pi/5))*sin(rp_star);cos(rp_star)]; %thrid vertex
% p=@(s,t) (s*v1+t*v2+(1-t-s)*v3);
% g_e=@(s,t) ((v1'*p(s,t))/(p(s,t)'*p(s,t))^2);
% 
% ener=abs(v3'*cross(v1,v2))*((1/40)*(g_e(0,0)+g_e(1,0)+g_e(0,1))+(1/15)*(g_e(0,0.5)+g_e(0.5,0)+g_e(0.5,0.5))+(9/40)*g_e(1/3,1/3));
% Ep=2*(A_Tp-ener); %Energy of one triangle of the pentagon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tests on floating point precisions
flag1=abs((N-12)*(6*A_Th)+12*5*A_Tp-4*vpa(Pi)); %showing that the area is equals the surface of the sphere
flag2=abs(theta_star*(N-12)+10*phi_star-vpa(Pi/3)*(N-2)); %the first equation of the system relating phi and theta must be satified
flag3=abs(sin(phi_star)-vpa(gamma)*sin(theta_star)); % the second equation must be satisfied. 

if flag1>tol_sym || flag2>tol_sym || flag3>tol_sym
    warning('Problem with the SYMBOLIC value of the Theoretical Energy Ground State')
end


%% Ground state total regularity measures;
Energyref=(N-12)*6*E_Th+12*5*E_Tp; %Total energetic contribution from 6 triangles per hexagon and 5 triangles per pentagon (in Matlab sybmbolic class)
Energyref_double=double(Energyref);%Total Energy (in Matlab double class)


MeanA2ref=(N-12)*(6*A_Th)^2+12*(5*A_Tp)^2;
MeanIPRref_flat=(N-12)*6*(L^2/A_Th)+12*5*(L^2/A_Tp);
MeanIPRref_sphere=(N-12)*(6*L)^2/(6*A_Th*(4*vpa(Pi)-6*A_Th))+12*(5*L)^2/(5*A_Tp*(4*vpa(Pi)-5*A_Tp));

GroundState=[Energyref MeanA2ref MeanIPRref_flat MeanIPRref_sphere];

if abs(Energyref-Energyref_double)>tol_double
    warning('Problem with the DOUBLE value of the Theoretical Energy Ground State')
end

digits(digitsOld) % DON'T FORGET TO CHANGE THE PRECISION TO THE ORIGINAL VALUE,
                  % OTHERWISE EVERY CALCULATION MADE OUTSIDE THIS FUNCTION
                  % WILL STILL CARRY THOSE DIGITS AND WILL RENDER THINGS
                  % VERY SLOW.
                
end

function value = arc_cosine (cc)

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