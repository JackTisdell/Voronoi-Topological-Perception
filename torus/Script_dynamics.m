%% For the plotting/graphics

disp_plot=1; %display plot frames of the animation

fig=figure(2); fig.Units='centimeters';  fig.Position=[0 5.9972 28.3633 22.3661];
%% Paramaters of the scenario

Nr=150; % number of agents moving to the RIGHT
Nl=150; %number of agents moving to the LEFT
Lx=10; %length of the hallway Omega
Ly=10; % width of the hallway Omega

Omega=[[0;Lx;Lx;0],[0;0;Ly;Ly]]; % hall way domain
N=Nr+Nl; %total number of agents;

sz=15; % plotting size for quiver object
%% Parameters of the model
L=1; % length scale of interaction and personal space
nu=2.5; % weighting parameter for the alignment component.
q=Inf; % topological fields of vision for the targets...if q~=Inf then we need p=q.
p=1; % topological fields of vision for averaging velocities....if q~=Inf then we need p=q.
if q~=Inf && p~=q
    error('Sizes "p" and "q" of the topological fields of vision are not compatible!!')
end

H_p=3*p*(p+1); % number of agents in topological field of vision of order p in the hexagonal lattice

%% Parameters of the implementation
dt=.5; %time scale
kMax=20; %maximal number of iterations

ee=exp(1);
rE=N^2*18*sqrt(3)/(5*(Lx*Ly)^2); %normalization factor for Voronoi energy
sectorx_basis=[0;0;2*Lx;2*Lx]; % To help define \overline(Wi)
sectory_basis=[Lx;-Lx;-Lx;Lx];
if q==Inf, chara='\infty'; else, chara=num2str(q); end %for the title of the figure.
c_del=1.5; % delta-vecinity factor for points to be consider by ID_copy (for coloring purposes), i.e we will use the color-delta neighborhood of the domain Omega
%% Initial configurations

rng(2)
xy=rand(N,2); 
x=Lx*xy(:,1); y=Ly*xy(:,2);

%% Prescribed velocities at time t=-1*dt


% 2) the unitary velocities have uniform random angle on [0,2pi)
angl=2*pi*rand(N,1);
vx_old=cos(angl); vy_old=sin(angl);

velx_old=vx_old; vely_old=vy_old;

%% Initializing data to save
DT_t = cell(kMax,1);
idx_ghost_t=cell(kMax,1);
idx_ghost2original_t = cell(kMax,1);
idx_copy_t=cell(kMax,1);
loc_alignment_t = zeros(kMax,N);

Vx_t=zeros(kMax,N);
Vy_t=zeros(kMax,N);
Vx_old_t=zeros(kMax,N);
Vy_old_t=zeros(kMax,N);
Area_t=zeros(kMax,N);
Angle_vi_t=zeros(kMax,N);
Polarization_t=zeros(kMax,1);
Energy_t=zeros(kMax,1);
AvgVel_t=zeros(kMax,1);
StdVel_t=zeros(kMax,1);


%% Loop over time

for k=1:kMax

%% Delaunay and Voronoi
[~,DT,V,C,idx_copy,idx_bdry] = periodicVoronoi(x,y,Omega);
x4=DT.Points(:,1); y4=DT.Points(:,2); % full set of generators 
N_copy=length(idx_copy);

DT_t{k} = DT;
idx_copy_t{k}=idx_copy;

idx_ghost_col=(N+5:length(x4))'; % indices w.r.t x4/y4 of all ghost copies created (NOT including the 4 ghosts at "infinity");
idx_ghost_t{k}=idx_ghost_col;

idx_ghost=reshape(idx_ghost_col,[N_copy 8]); % reshape of the indices such that each column represents one of the 8 copies that were created to obtain the periodic boundary conditions
idx_ghost2copy=mod(idx_ghost-N-4,N_copy); idx_ghost2copy(idx_ghost2copy==0)=N_copy; %matrix that maps every element of idx_ghost to the index copy (result is index w.r.t idx_copy)
idx_ghost2original=idx_copy(idx_ghost2copy); %transforms the indices of idx_ghost2copy into indices w.r.t x and y.
idx_ghost2original_col=idx_ghost2original(:);
idx_ghost2original_t{k}=idx_ghost2original_col;
%%
Neigh_Ui=cell(N,1); % the Union of Neigh with i itself, i.e in set notation {Neigh{i}}U{i}
Neigh_Ui_ghosts=cell(N,1);
idx_topoField=cell(N,1);

AdjMat=sparse(N,N); 
powerAdjMat=speye(N);
cumPowerAdjMat=sparse(N,N);

EdgeList= [DT(:,[1 2]); DT(:,[2 3]); DT(:,[3 1])]; %Concatenated idices of vertices that form each edge of the triangulation
EdgeList2= unique(EdgeList,'rows'); % contains 2 columns: for the starting-vertex and the end-vertex of each delaunay edge.

for i=1:N
    neigh_ghost=EdgeList2(EdgeList2(:,1)==i,2); %for each generator i, the list neigh contains the indices of Voronoi neighbors (including the ghost/copies)
    idx=find(neigh_ghost>N);
    idx_bool=ismember(idx_ghost_col,neigh_ghost(idx));
    
    
    neigh=neigh_ghost;
    neigh(idx)=idx_ghost2original_col(idx_bool);
    Neigh_Ui_ghosts{i}=[i; neigh_ghost];
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

%% velocity alignment component \tilda{a}

ax_tilda=zeros(N,1);
ay_tilda=zeros(N,1);
norm_a=zeros(N,1);
for i=1:N
    %%%%%% NOTE: vi_old and vj_olf need to be already normalized%%%%%
    temp=idx_topoField{i}; %temporary variable 
    idx_topoField2=temp(temp~=i); %we create a second idx_topoField index set but this time we remove agent i itself (for the displacement \tilda{a} the agent i will NOT average over its own velocity).
    ntopoField2=length(idx_topoField2); % number of agents in p-topoField (again...we don't count agent i itself).
    
    Vi_old=repmat([vx_old(i);vy_old(i)],[1 ntopoField2]); % old unit velocity vector (i.e at step t=(n-1)*dt) of agent i that is repeated/concatenated ntopoField times.
    Vj_old=[vx_old(idx_topoField2)'; vy_old(idx_topoField2)']; %old unit velocity vectors of all agents within p-topoField of agent i (once more, without considerin i itlsef).
    dotproduct=dot(Vi_old,Vj_old); %dot product of unit velocity vectors at time t=(n-1)*dt.
    
    theta_ij=acos(dotproduct); % angles in [0,pi] of the vi with every vj that lands in its p-topological field of vision.    
    g_ij=transitionFun(theta_ij,pi,3); %weights for the average of the velocities 
    
    phi=ntopoField2/H_p;
    a_tilda=sum(Vj_old.*g_ij,2)/ntopoField2; %2x1 vector \tilda{a} obtained by the weighted average of velocities in the p-topological field
    
    % norm of a_i without the scalinfg factor \phi_i for colouring of
    % the Voronoi regions
    norm_a(i)=sqrt(a_tilda(1).^2+a_tilda(2).^2);
    
    a_tilda=phi*a_tilda;
    ax_tilda(i)=a_tilda(1); %extract x and y coordinates separately.
    ay_tilda(i)=a_tilda(2);
end
loc_alignment_t(k,:) = norm_a(:);

%% Repulsion term u
idx_u_ghost=zeros(N,1); %for index w.r.t x4/y4 of agent from which we move away
%idx_u=zeros(N,1); % for index w.r.t x/y of the agent from which we move away
dist_u=zeros(N,1); %geodesic distance on the torus to the closest neighbour

for i=1:N
    neigh_ghosts=Neigh_Ui_ghosts{i}; neigh_ghosts(neigh_ghosts==i)=[];
    [dist_u(i),idx_NearNeigh]=min(sqrt((x(i)-x4(neigh_ghosts)).^2+(y(i)-y4(neigh_ghosts)).^2)); %closest neighbour among the candidates.
    idx_u_ghost(i)=neigh_ghosts(idx_NearNeigh); %
    %idx_u(i)=Neigh_Ui{i}(Neigh_Ui_ghosts{i}==idx_u_ghost(i)); % to use this way for extracting idx_u we need to have the option "stable"
                                                               % with the command "unique" when defining Neigh_Ui so that both
                                                               % Neigh_Ui and Neigh_Ui_ghost have the same corresponding ordering.
end

%We compute unitary direction for u with discontinuity at 0 (here we
% consider 1e-10 (roughly the order of precision allowed buy the Voronoi diagram)
ux=x-x4(idx_u_ghost); uy=y-y4(idx_u_ghost);
nor=sqrt(ux.^2+uy.^2); %Nx 1 vector of the norms of the vectors u
idx_nor_notsmall=find(nor>=1e-10); idx_nor_small=(1:N)'; idx_nor_small(idx_nor_notsmall)=[]; % finding indices for which we need to normalize and those for which we will not (we actually enforce ux=uy=0 for the latter case)
ux(idx_nor_notsmall)=ux(idx_nor_notsmall)./nor(idx_nor_notsmall); uy(idx_nor_notsmall)=uy(idx_nor_notsmall)./nor(idx_nor_notsmall); 
ux(idx_nor_small)=0; uy(idx_nor_small)=0;

r_L=transitionFun(dist_u,L,3); %the last argument is the flag specifying the decaying function we want to use.

%% Compute direction of velocity
gammax_dot=(r_L.*ux+nu*ax_tilda)./(r_L+nu);
gammay_dot=(r_L.*uy+nu*ay_tilda)./(r_L+nu); 

% set this as the old velocity for next iteration, plus we normalize it.
vx_old=gammax_dot; vy_old=gammay_dot;
nor=sqrt(vx_old.^2+vy_old.^2); %Nx 1 vector of the norms of the vectors u
idx_nor_notsmall=find(nor>=1e-10); idx_nor_small=(1:N)'; idx_nor_small(idx_nor_notsmall)=[]; % finding indices for which we need to normalize and those for which we will not (we actually enforce ux=uy=0 for the latter case)
vx_old(idx_nor_notsmall)=vx_old(idx_nor_notsmall)./nor(idx_nor_notsmall); vy_old(idx_nor_notsmall)=vy_old(idx_nor_notsmall)./nor(idx_nor_notsmall); 
vx_old(idx_nor_small)=0; vy_old(idx_nor_small)=0;

Vx_old_t(k,:)=vx_old';
Vy_old_t(k,:)=vy_old';

%% Computation of the densities \overline{W} and the individual speeds proper to each agent.
W_tilda=Inf(N,1);
V_tilda=cell(N,1);
angle_vi=zeros(N,1); %for angles w.r.t the horizontal of the vector v_i

Ener=zeros(N,1);

for i=1:N
    
    if ~isempty(intersect(i,idx_nor_small)) %There are two cases, when vi=0 and when it's not.
        [Ener(i),Area_t(k,i)]=poly_area_energy(V(C{i},1),V(C{i},2));
        W_tilda(i)=Area_t(k,i)/2;
    else 
        angle_vi(i)=atan2(vy_old(i),vx_old(i));
        sectorx=cos(angle_vi(i))*sectorx_basis-sin(angle_vi(i))*sectory_basis;
        sectory=sin(angle_vi(i))*sectorx_basis+cos(angle_vi(i))*sectory_basis;
        
        V_tilda{i}=sutherlandHodgman([V(C{i},1) V(C{i},2)],[sectorx+x(i) sectory+y(i)]);
        W_tilda(i)=poly_area(V_tilda{i}(:,1),V_tilda{i}(:,2));
        
        [Ener(i),Area_t(k,i)]=poly_area_energy(V(C{i},1),V(C{i},2));
    end
end
Energy_t(k)=mean(Ener)*rE;
Angle_vi_t(k,:)=angle_vi';

speed=tanh(2*ee*W_tilda/(pi*L^2));

%% Computing time series of global metrics
% Current velocity (not normalized)
velx=speed.*gammax_dot; vely=speed.*gammay_dot; %x and y coordinates of the N velocity vectors.
nor_vel=sqrt(velx.^2+vely.^2); %vector with the 2-norms of the N velocity vectors.

Vx_t(k,:)=velx';
Vy_t(k,:)=vely';

%a) Polarization
Polarization_t(k)=sqrt(sum(velx)^2+sum(vely)^2)/sum(nor_vel); % Polarization= norm of the sum of velocity vectors divided by the sum of the norm of the velocity vectors

%b) Average and standard deviation of the norm of the velocity vectors
AvgVel_t(k)=mean(nor_vel);
StdVel_t(k)=std(nor_vel);
%% PLOT alignement (New version using HSV colors where the Hue is the angle w.r.t the horizontal and the saturation is the area of the Voro cell)

if disp_plot==1
    
    %1) Domain
    scatter(0,0,'.w','MarkerEdgeColor','white')
    hold on
    rectangle('Position',[0 0 Lx Ly])
   
    
    % 2) Voronoi Diagram of the cylinder
    %     vor=voronoi(DT,'k-');
    %     vor(1).MarkerEdgeColor='none';
    %     vor(1).MarkerFaceColor='none';

    
    % 4) Unitary velocity vector field (simply to showcast heading of agents)
    quiv=quiver(x(1:Nr),y(1:Nr),vx_old(1:Nr)*Ly/100,vy_old(1:Nr)*Ly/100,0);
    quiv.ShowArrowHead='off';
    quiv.AutoScale='off';
    quiv.Color='k';
    quiv.LineWidth=1;
    
    quiv=quiver(x(Nr+1:end),y(Nr+1:end),vx_old(Nr+1:end)*Ly/60,vy_old(Nr+1:end)*Ly/60,0);
    quiv.ShowArrowHead='off';
    quiv.AutoScale='off';
    quiv.Color='k';
    quiv.LineWidth=1;
    
    % 5) plot all agents
    scatter(x(1:Nr),y(1:Nr),sz,'k.')
    scatter(x(Nr+1:N),y(Nr+1:N),sz,'k.')
    
    axis equal
    axis([0 Lx 0 Ly ])
    hold off
    yticks([0 L])
    yticklabels({'' 'L'})
    xticks([])
    title({['Torus with $N=',num2str(N),': L=',num2str(L),', p=',num2str(p),', q=',chara,', \nu=',num2str(nu),'$'];['$\#$iter $=',num2str(k-1),'$']},'Interpreter','latex')
    
    pause(0.000001)

    
    pic_count=pic_count+1;
    
end
%% Update position of agents and enforce boudary/door conditions.

x_new=x+dt*speed.*gammax_dot;
y_new=y+dt*speed.*gammay_dot;

x=mod(x_new,Lx);
y=mod(y_new,Ly);

disp(num2str(k))

end
%%

save("data.mat",'Nr','Nl','L','nu','q','p','Lx','Ly','DT_t','idx_copy_t','idx_ghost2original_t','idx_ghost_t','loc_alignment_t','Vx_t','Vy_t','Vx_old_t','Vy_old_t','Polarization_t','Area_t','Angle_vi_t','Energy_t','AvgVel_t','StdVel_t');
