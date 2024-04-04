%% For plotting/graphics
disp_plot=0; %display plot frames of the animation (1/0=Y/N)
disp_Earth=0; %display Earth (1/0=Y/N);

%% Paramaters of the scenario
Nr=150; % number of agents moving to the RIGHT
Nl=150; %number of agents moving to the LEFT
N=Nr+Nl; %total number of agents;

%% Parameters of the model
%L=3.46*2*sqrt(pi)/sqrt(N); %Euclidian Lengthscale
L=3.46*acos(1-2*pi/N); %Spherical Lengthscale
nu=1; % weighting parameter for the alignment component.
q=Inf; % topological fields of vision for the targets...if q~=Inf then we need p=q.
p=1; % topological fields of vision for averaging velocities....if q~=Inf then we need p=q.

%% Parameters of the implementation
dt=sqrt(pi)/10; %time scale comparable to our work on the torus
kMax=200; %maximal number of iterations

%% Allocating memory for data to save

% Voro/Delaunay
DT_t=cell(kMax,1);
V_t=cell(kMax,1);
C_t=cell(kMax,1);

%Position/Velocity
XYZ_t=zeros(N,3,kMax);
Vel_old_t=zeros(N,3,kMax);
Vel_t=zeros(N,3,kMax);

% Statistics and order parameters
Area_t=zeros(kMax,N);
Loc_alignment_t = zeros(kMax,N);
AngleVel_t=zeros(kMax,N);
Polarization_t=zeros(kMax,1);
Energy_t=zeros(kMax,1);
AvgVel_t=zeros(kMax,1);
StdVel_t=zeros(kMax,1);


%% Other stuff
if q~=Inf && p~=q
    error('Sizes "p" and "q" of the topological fields of vision are not compatible!!')
end
if L>pi
    error('length scale of interaction is larger than the the half circumference of the sphere')
end

if q==Inf, chara='\infty'; else, chara=num2str(q); end %for the title of the figure.

H_p=3*p*(p+1); % reference number of agents in topological field of vision of order p (in the hexagonal lattice)
ee=exp(1); % basis of the Neperian Logarithm. 
groundState=GroundState_Sphere(N); groundState=double(groundState); % Computation of the theoretical soccer ball gound state of the Voronoi energy. Obtained in MATLAB double class.

fig=figure(2); fig.Units='centimeters';  fig.Position=[0 2.328 28.045 26.035]; 
clf(fig);

addpath(genpath('PlotEarth')) %add path to the PlotEarth folder have the Earth as a graphic primitive
color1=[1 .42 0]; %orange 
color2=[0 .6 0]; % forest green
sz=20; % point size
%% Prescribed Positions and Velocities at time t=-1*dt

% Marsaglia uniform random points on the sphere
rng(3)
n1=floor(sqrt(N));n2=floor(sqrt(N)); n3=n1*n2; kk=3;
xx=2*rand(kk*n3,2)-1; x1=xx(:,1);x2=xx(:,2);
i=find(x1.^2+x2.^2<1);
if length(i)<N
    warning(['Need to produce more points for seeds# ',num2str(s)])
end
x11=x1(i);x22=x2(i);x33=1-2*(x11.^2+x22.^2);
x=2*x11(1:N).*sqrt(1-x11(1:N).^2-x22(1:N).^2);y=2*x22(1:N).*sqrt(1-x11(1:N).^2-x22(1:N).^2);z=x33(1:N);

% Objective: create randomly oriented Frenet basis of tangent plane at each x_i
XYZ=[x y z];
XYZ_old=XYZ; % for the sphere we need to know the agent's position at the current as well as at the previous time step (Note that we don't need this on flat spaces)

% 0) We find arbitrary POINTS belonging to the tangent planes of each
% agent's location. In order to avoid a division by 0 when using the
% algebraic form of these tangent planes, we extract the index of the
% coordinates (i.e. 1,2 or 3) of each vector in XYZ that achieve the
% highest Infinity norm of, i.e. the largest absolute component (that is
% automatically nonzero as each point in XYZ is on the unit sphere). The
% other components of the points are set to be all ONE. Then we simply use
% the algebraic expression of the tangent planes to express the remaining
% component of the point
[~,idxInfNorm]=max(abs(XYZ_old),[],2); idxInfNorm_lin=sub2ind([N 3],(1:N)',idxInfNorm);
temp1=ones(N,3); temp1(idxInfNorm_lin)=zeros(N,1); %setting to ONE the other components that did not achieve the Infinity norm.
temp2=(1-dot(XYZ_old',temp1'))'./XYZ_old(idxInfNorm_lin); % using algebraic equation of the plane to recover the constained remaining component (the one achieving the Infinity norm)...here is where we needed to guarantee a nonzero divisor.

Frenet1_old=temp1; Frenet1_old(idxInfNorm_lin)=temp2; % arbitrary POINTS belonging to the tangent planes of each x_i.

% 1) We obtain arbitrary VECTORS belonging to the tangent planes by
% substracting the agent's coordinates
Frenet1_old=Frenet1_old-XYZ_old; % aribitrary VECTORS belonging to the tangent planes of each x_i. These represent the 1st cardinal vector of the Frenet reference frames (the coordinates are given in terms of the fixed/inertial reference frame)
nor=sqrt(sum(Frenet1_old.^2,2)); Frenet1_old=Frenet1_old./nor;%normalization of the first cardinal vector of the Frenet reference frames. NOTE: we are divinding a Nx3 array by a Nx1 vector...each normalization factor is applied to a row of 3 elements.

% 2) Obtain the second Frenet basis vectors by cross product between the
% normal vectors and the first Frenet vectors.
Frenet2_old=cross(XYZ_old,Frenet1_old);

% NOTE: we can verify that we have indeed obtained an orthonormal Frenet basis
% for the tangent plane at each agent by computing the following
% dot(Frenet1_old',XYZ')' % should be all zeros up to machine precision
% dot(Frenet2_old',XYZ')' % should be all zeros up to machine precision
% sqrt(Frenet1_old(:,1).^2+Frenet1_old(:,2).^2+Frenet1_old(:,3).^2) % should be all ones up to machine precision
% sqrt(Frenet2_old(:,1).^2+Frenet2_old(:,2).^2+Frenet2_old(:,3).^2) % should be all ones up to machine precision

% 3) Obtain random unitary vectors belonging to the tangent planes by
% trigonometric linear combination of the Frenet basis vectors
angles=2*pi*rand(N,1);
Vel_old=Frenet1_old.*cos(angles)+Frenet2_old.*sin(angles);

% 4) Obtain the new XYZ postions after a rotation by and angle
% dt*||Vel_old||
XYZ=cos(dt*1).*XYZ_old+sin(dt*1).*Vel_old; % the 1s come from the fact that for the first update ww use Unitary velocities (this will not be true for subsequent updates)

%% Iterations over time

for k=1:kMax
    
%% Voronoi/Delaunay and p-topological field of view.
[V,C,DT]=sphereVoronoi(XYZ');
DT_t{k}=DT; V_t{k}=V; C_t{k}=C;
XYZ_t(:,:,k)=XYZ;

Neigh_Ui=cell(N,1); % the Union of Neigh with i itself, i.e in set notation {Neigh{i}}U{i}
idx_topoField=cell(N,1);

AdjMat=sparse(N,N); 
powerAdjMat=speye(N);
cumPowerAdjMat=sparse(N,N);

EdgeList= [DT(:,[1 2]); DT(:,[2 3]); DT(:,[3 1])]; %Concatenated idices of vertices that form each edge of the triangulation
EdgeList2= unique(EdgeList,'rows'); % contains 2 columns: for the starting-vertex and the end-vertex of each delaunay edge.
for i=1:N
    neigh=EdgeList2(EdgeList2(:,1)==i,2); %for each generator i, the list neigh contains the indices of Voronoi neighbors (including the ghost/copies)
    Neigh_Ui{i}=unique([i; neigh],'stable');    
    AdjMat(i,Neigh_Ui{i})=ones(1,length(Neigh_Ui{i}));
end

% artifact to make AdjMat symmetric in case the ordering from EdgeList
% skrew things up and we have an (i,j) dependency but not (j,i) or vice
% versa.
symmDiff=AdjMat-AdjMat';
if nnz(symmDiff)~=0
    %warning('AdjMat was forcibly made symmetric')
    idx_symmDiff=find(symmDiff~=0); AdjMat(idx_symmDiff)=1;
end

if p==1
    idx_topoField=Neigh_Ui;
else
    for i=1:p
        powerAdjMat=powerAdjMat*AdjMat;
        cumPowerAdjMat=cumPowerAdjMat+powerAdjMat;
    end
    
    for i=1:N
        idx_topoField{i}=find(cumPowerAdjMat(i,:)~=0);
        %NOTE: idx_topoField{i} contains i itself since AdjMat^2 has +1 in
        %the diagonals and thus cumPowerAdjMat has always nonzero diagonals.
    end
end

%% Repulsion term u
x=XYZ(:,1); y=XYZ(:,2); z=XYZ(:,3);
idx_u=zeros(N,1); % for index w.r.t x/y of the agent from which we move away
distEucl_u=zeros(N,1);% Euclidian distance on the sphere to the closest neighbour (bounded above by 2)
for i=1:N
    neigh=Neigh_Ui{i}; neigh(neigh==i)=[];
    [distEucl_u(i),idx_NearNeigh]=min(sqrt((x(i)-x(neigh)).^2+(y(i)-y(neigh)).^2+(z(i)-z(neigh)).^2)); %index of closest neighbour among the candidates (NOT w.r.t the indexing of xyz)
    idx_u(i)=neigh(idx_NearNeigh); %index of closest neighbor for each agent (w.r.t the indexing of xyz) 
end
dist_u=acos(1-(distEucl_u.^2)/2); % Geodesic distance between each agent and it closest neighbor, equivalent to acos(dot(XYZ',XYZ(idx_u,:)'))'

Uxyz=XYZ-XYZ(idx_u,:); %Euclidian vector pointing AWAY from the closest neighbor
Uxyz=Uxyz-dot(Uxyz',XYZ')'.*XYZ; %Projection of the Euclidian vectors onto the tangent planes.
                                 % Already verified that dot(Uxyz',XYZ')'=zeros(N,1) up to machine
                                 % precision, thus the vectors Uxyz belong
                                 % indeed to the tangent planes after the
                                 % last update.

nor=sqrt(sum(Uxyz.^2,2)); %we normalize the repulsion vectors
idx_nor_small=find(nor<=1e-10);
nor(idx_nor_small)=ones(size(idx_nor_small));
Uxyz(idx_nor_small,:)=zeros(length(idx_nor_small),3);
Uxyz=Uxyz./nor;

r_L=transitionFun(dist_u,L,3); %the last argument is the flag specifying the decaying function we want to use.

%% Velocity alignment
norm_a=zeros(N,1);
Axyz_tilda=zeros(N,3);
for i=1:N
    %% Extracting all agents in the p-topological field of view.
    % NOTE: In the flat case we excluded agent i itself from it's own field
    % of view. For the sphere we do not include it either, however some
    % computations needed between x_i and x_i_old are the same as the ones
    % needed between x_i and x_j_old (j in the topoField of view). Thus it is
    % practical to include x_i itlsef in its own field of view, simply to
    % ease computations.
    idx_topoField2=idx_topoField{i}; %we create a second idx_topoField index set but this time we remove agent i itself (for the displacement \tilda{a} the agent i will NOT average over its own velocity).
    idx_i=find(idx_topoField2==i); % idx of agent i itself w.r.t the list idx_topoField2 now that we are including agent i.
    ntopoField2=length(idx_topoField2)-1; % number of agents in p-topoField (here we do not count agent i).
    
    %% Computing the aligent for each agent.
    
    % 0) Defining useful quantities
    XYZ_old_TopoField=XYZ_old(idx_topoField2,:); % old location of current topoField agents at time t-dt will be used throughout the section so we declare it as a new variable to speed things up.
    Vel_old_TopoField=Vel_old(idx_topoField2,:); % old velocities of current topoField agents at time t-dt. We also declare it as a variable for easier use.
    vec_Eucl=XYZ(i,:)-XYZ_old_TopoField; % Euclidian vectors separating x_i from every x_j_old according to the current topological field of view (including x_i_old as explained earlier)
    RepXYZi=repmat(XYZ(i,:),ntopoField2+1,1); %repeated xyz component of x_i a total number of ntopoField2+1 times.
      
    
    % 1) Compute several different Frenet frames at x_i, each one having
    % the 1st axis pointing AWAY from a different x_j_old.
    Frenet1=vec_Eucl-dot(vec_Eucl',RepXYZi')'.*RepXYZi;
    nor=sqrt(sum(Frenet1.^2,2)); Frenet1=Frenet1./nor; %normalization of the 1st Frenet axis for each of the ntopoField2 frames (each axis points AWAY from one of the topoField agents).
    Frenet2=cross(RepXYZi,Frenet1); % 2nd Frenet axis for each agent in the topoField. Already unitary vectors.
    
    % 2) Compute the Frenet frames at x_j_old for each agent within the topoField of
    % view of x_i; the 1st axis needs to point TOWARDS x_i.
    % NOTE: because we aim to use the old velocities (time t-dt) we need to
    % compute these Frenet frames w.r.t to the location of the topofield
    % agents at time t-dt as well
    
    Frenet1_TopoField=vec_Eucl-dot(vec_Eucl',XYZ_old_TopoField')'.*XYZ_old_TopoField; % 1st Frenet axis of each agent in the topoField of x_i (belonging to the local tangent plane). It points TOWARDS x_i.
    nor=sqrt(sum(Frenet1_TopoField.^2,2)); Frenet1_TopoField=Frenet1_TopoField./nor; %normalization of the 1st Frenet axis for each agent in the topoField. 
    Frenet2_TopoField=cross(XYZ_old_TopoField,Frenet1_TopoField); % 2nd Frenet axis  for each agent in the topoField. Already unitary vectors.
    
    % 3) We decompose the old velocity of each topoField agent according to
    % the local Frenet frames we just computed; We will later use the
    % corresponding coordinates. NOTE: V_old needs to contain UNIT
    % VELOCITY vectors from the time t-dt.
    coord1=dot(Vel_old_TopoField',Frenet1_TopoField')';  % coordinate of Vel_old_TopoField w.r.t Frenet1_neigh (every vector is unitary so the coordinate is given by cos(theta)=dot(Vel_old_TopoField',Frenet1_neigh')'
    coord2=dot(Vel_old_TopoField',Frenet2_TopoField')';  % coordinate of Vel_old_TopoField w.r.t Frenet2_neigh
    
    % 4) So far we have created ntopoField2 different Frenet frames for x_i, each
    % having an equivalent Frenet frame at each topoField agent x_j_old. By equivalent we mean  
    % that they are obtained from one another by rotation along the piece
    % of great circle joining x_i to x_j_old (for each j in the current topoField).
    % Because of these equivalences, the components of the vel_j_old
    % are the same for the Frenet frame of x_j_old as they are for the Frenet
    % frame of x_i.
    Vel_old_TopoField=Frenet1.*coord1+Frenet2.*coord2; % unitary velocity vectors of each topoField agent but now belonging to the tagent plane at x_i.
    
    % 5) Now that every velocity vector has been brought to the tangent
    % plane of x_i we can compute angles and add them up for the weighted
    % average. First we extract the old velocity vector of agent x_i and
    % erase it from the list of velocity vectors from topoField agents
    % (Remember we only included i itself in the list because it was
    % computationally easier)
    Vel_old_i=Vel_old_TopoField(idx_i,:);
    Vel_old_TopoField(idx_i,:)=[];
    
    RepVel_i=repmat(Vel_old_i,ntopoField2,1); %repeated velocity component of x_i_old that were translated to the tangent plane of x_i (we repeat a total number of ntopoField2 times).
    theta_ij=acos(dot(Vel_old_TopoField',RepVel_i'))'; % angles in [0,pi] between v_i_old and every v_j_old with j in the p-topological field of vision.    
    g_ij=transitionFun(theta_ij,pi,3); %weights for the average of the velocities 

    a_til=sum(g_ij.*Vel_old_TopoField)/ntopoField2; %2x1 vector \tilda{a} obtained by the weighted average of velocities in the p-topological field
    
    % norm of a_i without the scalinfg factor \phi_i
    norm_a(i)=sqrt(sum(a_til.^2));
    
    % 6) Final alignment component of displacement 
    phi=ntopoField2/H_p;
    Axyz_tilda(i,:)=phi*a_til;

end
Loc_alignment_t(k,:)=norm_a(:);

%% Compute direction of velocity
Gamma_dot=(r_L.*Uxyz+nu*Axyz_tilda)./(r_L+nu);

% Update the old unitary velocities.
Vel_old=Gamma_dot; %update
nor=sqrt(sum(Vel_old.^2,2)); %we normalizization
idx_nor_small=find(nor<=1e-10);
nor(idx_nor_small)=ones(size(idx_nor_small));
Vel_old(idx_nor_small,:)=zeros(length(idx_nor_small),3);
Vel_old=Vel_old./nor;

% Save the data over time
Vel_old_t(:,:,k)=Vel_old;


%% Computation of the densities \overline{W} and the individual speeds proper to each agent.
W_tilda=Inf(N,1);

[~,~,~,Area_t(k,:),Energy_t(k),~]=poly_centroid_area_scaledEnergy_sphere(x,y,z,V,C,groundState(1));

for i=1:N
    
    if ~isempty(intersect(i,idx_nor_small)) %There are two cases, when vi=0 and when it's not.
        W_tilda(i)=Area_t(k,i)/2;
    else 
        % Our goal here is to consider the "forward" part of the Voronoi
        % cell in the direction of the current velocity vector v_i.
        % For this we consider the plane passing through x_i and orthogonal
        % to v_i. The Voronoi vetices are either ABOVE or BELOW the plane.
        % We will need to insert 2 new Voronoi vertices where there is a
        % crossing between above and below. Finally we'll remove all
        % vertices that are below the plane. The remaining spherical
        % polygon will be the "forward part" we are looking for.
        
        Vxyz=[V(C{i},1) V(C{i},2) V(C{i},3)]; %stacked coordinates of the vertices of the Voronoi cell of Vi. Easier for computations. 
        Vel_old_i=Vel_old(i,:); % declare the velocity, we'll use it several times. NOTE: because we updated V_old right in the previous section we have that Vel_old(i,:) is the CURRENT velocity and NOT the one from t-dt.
        Nneigh=length(Vxyz); %number of Voronoi vertices of the cell of x_i
        RepVel_i=repmat(Vel_old_i,Nneigh,1); % Repeating the unit velocity vector v_i 
                                                   
        
        % 0) We compute the dot product between the Voronoi vertices and the CURRENT unit velocity vector of x_i.   
        % The Voronoi vertices for which test is positive are the ones ABOVE the plane passing through x_i and orthogonal to v_i.
        test1=dot(RepVel_i',Vxyz')';  %if positive then the vertex is above and if negative then the vertex if below.
        test2=test1.*[test1(2:end);test1(1)]; %if the product is negative then there is a crossing 
        idx_crossing=find(test2<0); %there is a crossing of the plane between vertex idx_crossing and vertex idx_crossing+1
        
        if length(idx_crossing)~=2 %verifying that there is only two edges that cross
            warning('Wrong NUMBER of crossing edges')
        end
        
        % 1) We indentify which of the two crossings is below--> above
        % then we shift the order of the vertices such that 
        % the crossing below--> above happens between the 1st and 2nd vertex
        % REMEMBER that the function sphereVoronoi.m gives the Voronoi
        % vertices in ccw orientation when looking from outside the sphere.
         
        if test1(idx_crossing(1))<0 
            shift_idx=1; %For this case the crossing below-->above happens between vertices idx_crossing(1) and idx_crossing(1)+1
        else
            shift_idx=2; %For this case the crossing below-->above happens between vertices idx_crossing(2) and idx_crossing(2)+1
        end
        
        % We shift the list of vertices and their classification
        % "above/below" according to the shift index.
        Vxyz=circshift(Vxyz,-idx_crossing(shift_idx)+1);
        test2=circshift(test1,-idx_crossing(shift_idx)+1);
        
        if ~(test2(1)<0 && test2(2)>0) %check that there is indeed a crossing below-->above between the 1st and 2nd vertices in the new ordering.
            warning('Wrong ORDERING for crossing vertices')
        end
        
        idx_above=find(test2>0); %update the index of vertices above;
        
        % 2) We compute the FIRST crossing: as already
        %  mentioned, it happens between the first and second vertices in
        %  the new ordering.
        V_start=Vxyz(1,:); % starting vertex of the edge that crosses below-->above (remember the ccw ordering)
        V_end=Vxyz(2,:);  %  vertex where the edge ends
        
        t_star=(Vel_old_i*V_start')/(Vel_old_i*(V_start-V_end)'); % value of the parameter t for which the edge V_start-->V_end crosses the separating plane.
        
        if t_star<0 || t_star >1
            warning('Problem with PARAMETRIZATION of 1st crossing edge')
        end
        
        V_crosses1=V_start+t_star*(V_end-V_start); % First new vertex obtained by pluging t_star into the parametric form of the edge V_start-->V_end
        V_crosses1=V_crosses1/sqrt(sum(V_crosses1.^2)); %project it onto the sphere.
        
        % 3) We compute the SECOND crossing:
        % We obtain the starting vertex of this crossing edge as the LAST
        % index that is above, i.e. idx_above(end).
        % NOTE: in the event that there is ONLY ONE vertex below
        % (located automatically in first place in the new ordering) we need to
        % loop around the list of vertices so that the ending vertex is the
        % is the first one (vertex below the plane)
        
        V_start=Vxyz(idx_above(end),:);
        if idx_above(end)==Nneigh %case when there is only one vertex below (
            V_end=Vxyz(1,:);
        else
            V_end=Vxyz(idx_above(end)+1,:);
        end
        
        t_star=(Vel_old_i*V_start')/(Vel_old_i*(V_start-V_end)'); % value of the parameter t for which the edge V_start-->V_end crosses the separating plane.
        
        if t_star<0 || t_star >1
            warning('Problem with PARAMETRIZATION of 2nd crossing edge')
            
        end
        
        V_crosses2=V_start+t_star*(V_end-V_start); % First new vertex obtained by pluging t_star into the parametric form of the edge V_start-->V_end
        V_crosses2=V_crosses2/sqrt(sum(V_crosses2.^2)); %project it onto the sphere.
        
        % 3) We The vertices of V_tilda and compute its area  
        V_tilda=[V_crosses1; Vxyz(idx_above,:); V_crosses2]; % Complete list of vertices of the halfed Voronoi region.
        W_tilda(i)=poly_area_sphere(V_tilda); % surface area of the spherical polygon V_tilda.
        
        if W_tilda(i)>Area_t(k,i)
            warning('Area of halfed Voronoi cell has a problem')
        end
    end
end

speed=tanh(ee*W_tilda/(pi*(1-cos(L)))); % scaling w.r.t half the area of a spherical cap of geodesic radius L.

%% Plot/Diaplay Graphics
if disp_plot==1
    if k==1
        hold on
        if disp_Earth==1
            Earth=plotearth('MapType','bluemarble');
            Earth.View=[151.4 23.8000];
        else
            [x_sph,y_sph,z_sph]=sphere(100);
            plot_sph=surf(0.980*x_sph,0.980*y_sph,0.980*z_sph,'FaceLighting','gouraud','EdgeColor','none','FaceColor','black','FaceAlpha',.6);
            lightangle(60,20)
        end
        
        quiv=quiver3(XYZ(:,1),XYZ(:,2),XYZ(:,3),Vel_old(:,1),Vel_old(:,2),Vel_old(:,3),0.15);
        quiv.ShowArrowHead='off'; quiv.AutoScale='on'; quiv.Color=color1; quiv.LineWidth=1.8;
        %quiv.Marker='o'; quiv.MarkerSize=4;quiv.MarkerFaceColor=color1;
        
        hold off
        
        axis equal
        axis manual
        axis off
        titl=title({['Sphere with $N=',num2str(N),': L=',num2str(L),', p=',num2str(p),', q=',chara,', \nu=',num2str(nu),'$'];['$\#$iter $=',num2str(k-1),'$']},'Interpreter','latex');
    end
    
    quiv.XData=XYZ(:,1); quiv.YData=XYZ(:,2); quiv.ZData=XYZ(:,3); %update quiver's arrow starting point
    quiv.UData=Vel_old(:,1); quiv.VData=Vel_old(:,2); quiv.WData=Vel_old(:,3); %update quiver's arrows directions
    titl.String{2}=['$\#$iter $=',num2str(k-1),'$']; %update second line of the title (to update iteration number
    drawnow limitrate % display updates
end

%% Computing time series of global metrics
% Current velocity (not normalized)
Vel=speed.*Gamma_dot;
nor_Vel=sqrt(sum(Vel.^2,2));

% Save current velocity.
Vel_t(:,:,k)=Vel;

%a) Polarization in terms of Angular Velocity
% NOTE: since we are on the unit sphere, the norm of Ang_Vel is the same as
% the norm of the Vel vectors
Ang_Vel=cross(XYZ,Vel); %vectors of Angular velocities (The expression simplifies a bit because XYZ are on the unit sphere and Vel are tangent to the sphere).
Polarization_t(k)=sqrt(sum(sum(Ang_Vel).^2))/sum(nor_Vel); % Polarization= norm of the sum of ANGULAR VELOCITIES vectors divided by the sum of the norm of the ANGULAR VELOCITIES vectors

%b) Average and standard deviation of the norm of the velocity vectors
AvgVel_t(k)=mean(nor_Vel);
StdVel_t(k)=std(nor_Vel);

%c) deviation from average angular velocity;
Avg_Ang_Vel=mean(Ang_Vel); Avg_Ang_Vel=Avg_Ang_Vel/sqrt(sum(Avg_Ang_Vel.^2));
Rep_Avg_Ang_Vel=repmat(Avg_Ang_Vel,N,1);
Ang_Vel=Ang_Vel./nor_Vel;
AngleVel_t(k,:)=acos(dot(Rep_Avg_Ang_Vel',Ang_Vel'));


% % c) Angles of velocity vectors. We need to bring all vectors in Vel to a
% % common tangent plane, WLOG let's take that tangent plane to be at the
% % NORTH POLE=[0 0 1]. It simplifies calculations a lot.
% Rep_NP=repmat([0 0 1],N,1); % repeated coordinates of the North Pole
% vec_Eucl=[0 0 1]-XYZ; % Euclidian vectors separating NorthPole from every x_i
% Frenet1_NP=[vec_Eucl(:,1:2) zeros(N,1)]; %First Frenet axis at the North Pole pointing towards each of of the N agents, because of NORTH POLE=[0 0 1] the projection of vec_Eucl onto the Tangent planes simplfify drastically.
% nor=sqrt(sum(Frenet1_NP.^2,2)); Frenet1_NP=Frenet1_NP./nor; %normalization of the 1st Frenet axis at the North Pole
% Frenet2_NP=cross(Rep_NP,Frenet1_NP); % 2nd Frenet axis at the North Pole
% 
% Frenet1_Agent=vec_Eucl-dot(vec_Eucl',XYZ')'.*XYZ; % 1st Frenet axis of each agent belonging to the tangent plane at each agent. They all points TOWARDS North Pole.
% nor=sqrt(sum(Frenet1_Agent.^2,2)); Frenet1_Agent=Frenet1_Agent./nor; %normalization of the 1st Frenet axis for each agent
% Frenet2_Agent=cross(XYZ,Frenet1_Agent); % 2nd Frenet axis for each agent. Already unitary vectors.
% 
% coord1=dot(Vel_old',Frenet1_Agent')';  % 1st coordinate of unitary velocity w.r.t Frenet1_Agent (every vector is unitary so the coordinate is given by cos(theta)=dot(Vel_old',Frenet1_Agent')'
% coord2=dot(Vel_old',Frenet2_Agent')';  % coordinate of Vel_old_TopoField w.r.t Frenet2_neigh
% Vel_old_NP=Frenet1_NP.*coord1+Frenet2_NP.*coord2; % unitary velocity vectors of each agent but now belonging to the tagent plane at the North Pole.

% AngleVel_t(k,:)=atan2(Vel_old_NP(:,2),Vel_old_NP(:,1))';

%% Update position of agents and enforce boudary/door conditions.

XYZ_old=XYZ; %Update the old position vectors
XYZ=cos(dt*nor_Vel).*XYZ+sin(dt*nor_Vel).*(Vel./nor_Vel); % Udpate the current position vector

disp(num2str(k))

end

save("data.mat",'dt','Nr','Nl','L','nu','q','p','DT_t','V_t','C_t','XYZ_t','Vel_t','Vel_old_t','Loc_alignment_t','Polarization_t','Area_t','Energy_t','AngleVel_t','AvgVel_t','StdVel_t');
