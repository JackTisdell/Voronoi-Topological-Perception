%% Implementation of bidrectional flow on a the straight hallway using version 1.5 of the VTP model
%
% Script_dynamics.m runs the time evolution of the system, it calculates
% all important dynamical quantities and metrics over time and stores it all in the "hallwayData.mat" file.
%
% To obtain a complete frame-by-frame animation, the user should refer to the accompanying
% Script_Graphics.m file that loads the "hallwayData.mat" file to produce each frame.
% However, for convenience and easier appreciation of the dynamics, the user can choose between displaying the state of the system
% at each iteration or to fully ommit any visualization (by letting setting disp_plot true or false).
%
% INPUT:
%   The fundamental degrees of freedom and input variables the user can play with:
%   - L, the repulsive falloff distance (preferred inter-personal distance)
%   - nu, the alignment coefficient
%   - Lc, the preferred inter-personal distance for agents about to enter the
%     hallway. Lc effectively captures the same degree of freedom as N
%   - p, the size of the topological field of view for alignment; i.e.
%     the alignment component is obtained by averaging the velocities
%     of all neighbors within a Delaunay topological distance of p from
%     agent x_i. The DEFAULT value is p=1.
%   - q, the size q of the topological field of view for awareness of targets;
%     i.e., agent x_i will have a nonzero homing component only if its target set
%     is at most at a Delaunay topological distance of q. The DEFAULT value
%     is q=INF so that every agent is capable to "see" its target set from
%     anywhere in the hallway.
%   - dt, the time scale built into the gouverning equations. REFERENCE value
%     is dt=1.
%
%    NOTE: changing p and q from their default values p=1 and q=Inf can substantially
%    increase computation time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Ivan Gonzalez and Jack Tidsell at McGill University
%
% Created June 2021
% Last Update March 2023:
% -added the computation of queueing structures and metrics
% -stored all dynamical data and metrics in MATLAB structures rather than
% cell arrays for easier later access 
% -made the plotting of the hallway more efficient by updating graphic
% handles rather than deleting them and creating new ones at every
% iteration
% -general clean up and commenting of the code
%
% should any bug or error arise please contact
% ivan.gonzalez@mail.mcgill.ca and jack.tisdell@mail.mcgill.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO:
% - change method for normalization of vectors using the
% temp=[x y]; issmall=sum(abs(v),2)<tol;  v=normalize(v,"norm"); v(issmall)=0;
% - change calculations of distances to use @hypot 
% - create separate function for the entrance stochastic process .
% - simplifying the calculation of repulsion component of displacement
% -change name for Neigh_Ui
% -vectorize the calculations for the homing component when q~=Inf using a
%   cellfun 
% - vectorize calculations of area and Voronoi energy using cellfun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the optional display of graphics (logical Yes/No options)

disp_plot=true; %display plot frames of the animation: "true" to display, "false" to ommit
save_frames=false; % "true" if you want to save the frames, false for otherwise. (NOTE: if disp_plot=false then save_frames is futile)
computeQueueing=false; % Do you wish to compute the queueing structures and associated metrics? These are time consuming to compute. 

%% Declaring parameter values

nu=2; % weighting parameter for the alignment component.
L=1; % falloff distance (preferred radius of personal space)
Lc=3; % lenght-scale for the maximal congestion allowed near each entrance.

p=1; % topological fields of vision for averaging velocities....if q~=Inf then we need p=q for compatibility. DEFAULT VALUE is p=1.
q=Inf; % topological fields of vision for the targets...if q~=Inf then we need p=q. DEFAULT VALUE is q=Inf.

dt=.5; % time scale. REFERENCE VALUE is dt=1.
randSeed=1; %seed of pseudo-random generator used to generate agents,

%% Declaring parameters of the scenario
Nr_total=2000; % maximal total number of agents moving to the right that can be inserted on the left boundary;
Nl_total=2000; % maximal total number of agents moving to the left that can be inserted on the right boudary;

kMax=20; %maximal number of iterations. The dynamics will run from time step k=1 to k=kMax or until only one agent is left in the halllway

Lx=60; % length (horizontal) of the hallway
Ly=12; % width (vertical) of the hallway (labelled "l" in the original paper)

rng(randSeed); %fixing the random seed generator used in the entry process for reproductibility. For different results change the argument in rng(.)
%% Other useful quantities to declare

N_total=Nr_total+Nl_total; % maximal total number of agents to be inserted over time (At both left and right entrances).

Omega=[[0;Lx;Lx;0],[0;0;Ly;Ly]]; % stacked coordinates of the rectangular hallway's corners
sectorx_basis=[0;0;2*Lx;2*Lx]; sectory_basis=[Lx;-Lx;-Lx;Lx]; % to help computing Area_forward

tol_exit=1e-5; % numerical tolerance for entrance/exit (any value between 1e-3 and 1-7 will do)
H_p=3*p*(p+1); % number of agents in topological field of vision of order p in the hexagonal lattice (quantity needed in the scaling of the alignment component)
rE_factor=18*sqrt(3)/(5*(Lx*Ly)^2); % constant factor used to scale Voronoi energy (clustering metric)

% Discretized left and right doors to use as the domain of the probability density functions (PDF) that govern the agent's entrance distribution
res_PDF=150; %resolution of the discretized domain of the PDF (range of the random variable)
Y=linspace(0,Ly,res_PDF)'; Y(end)=[]; Y=Y+Y(2)/2; % discretized domain of the PDF.

%% Catching possible error and question dialogue with user about graphics

if q~=Inf && p~=q %error message in case p and q parameters are incompatible
    error('Sizes "p" and "q" of the topological fields of vision are not compatible! Please set q=Inf or choose INTEGER values s.t. p=q')
end

if disp_plot==true
    answerQ= questdlg('The state of the system will be plotted at each time step, are you sure you wish to plot?','Watch Out','Yes','No','No');
    switch answerQ
        case 'Yes'
            if save_frames==true
                answerQ2 = questdlg('Are you sure you wish to save the plot frame at every time step?','Watch Out','Yes','No','No');
                if strcmp(answerQ2,'No')
                    save_frames=false;
                else
                    if ~isfolder('saved_frames') %if there is no subfolder called 'saved_frames', create one
                        mkdir saved_frames
                    end
                    save_dirctry=strcat(pwd,'/saved_frames'); % directory for the 
                end
            end
            
            % Create RGB color triplets using the rgb.m function (to plot agents and their velocity vectors).
            addpath(genpath('rgbColors')) %add path to the rgb function in the rgbColors folder
            quivColorNr=rgb('OrangeRed'); quivColorNl=rgb('DarkGreen'); %colors for plotting agents and their velocity vectors 
            forestColorNr=rgb('Orange'); forestColorNl=rgb('MediumSeaGreen'); %colors for plotting the queueing structures (forests).
            arrwSize=Ly/40;% size of arrows in quiver object to represent velocity 

            % Create graphic object handles
            fig=figure(1); clf; %create figure object
            fig.Units='centimeters'; fig.Position=[4.5 13.5 37.5 11.5]; % fix some of its properties
            axeFig=axes(fig); %creating an axes object such that axeFig.Parent=fig
            axeFig.Units='normalized'; axeFig.Position=[0.01 0.01 0.980 0.99]; %  normalized units represent a proportion of the Parent object fig
            
            % Writing down a title to the figure.
            if q==Inf, chara='\infty'; else, chara=num2str(q); end  %for the title of the figure.
            titleFig=title({['Bi-directional flow with parameter values $\nu=',num2str(nu),', L=',num2str(L),', L_c=',num2str(Lc),' , p=',num2str(p),' , q=',chara,'$']},'Interpreter','latex'); %we leave the third line of the title as an empty character that will be updated at every iteration with important information of the current iteration
            
            rectngl=rectangle(axeFig,'Position',[0 0 Lx Ly]); % draws the boundary of the hallway using a rectangle graphics object
            axeFig.NextPlot='add'; % equivalent to "hold on": in order to add new graphic objects as children of axeFig. To see the difference type axeFig.Children when toggling between axeFig.NextPlot='add' and axeFig.NextPlot='replace'.
           
            %create dummy quiver objects for position and velocity of out agents: will be updated at every iteration
            quivNr=quiver(axeFig,-1,-1,0,0); set(quivNr,'Color',quivColorNr,'Marker','.','MarkerSize',Ly,'ShowArrowHead','off','AutoScale','off','LineWidth',1); % for population moving RIGHT
            quivNl=quiver(axeFig,-1,+1,0,0); set(quivNl,'Color',quivColorNl,'Marker','.','MarkerSize',Ly,'ShowArrowHead','off','AutoScale','off','LineWidth',1); % for population moving LEFT
            
            % dummy line object to plot the queueing structures.
            % NOTE: we chose to use a line graphics object rather that a
            % GraphPlot object as the latter cannot be fully updated at every
            % iteration in order to plot the changes.
            queuePltNr=line(axeFig,-1,-1,'Color',forestColorNr,'LineWidth',1); 
            queuePltNl=line(axeFig,-1,-1,'Color',forestColorNl,'LineWidth',1);
            
            % dummy line object to plot the bi-direction Voronoi interface
            lineInterface=line(axeFig,-1,-1,'Color','magenta','LineWidth',1.5);
            
            % dummy line objects to plot the Delaunay triangulation and Voronoi
            % tessellation of the whole population
            delaunayPlt=line(axeFig,-1,-1,'Color','blue','LineStyle',':');
            voronoiPlt=line(axeFig,-1,-1,'Color','black','LineStyle',':');
           
            % text objects for the labels of agents
            labelHandleNr=text(axeFig,0,0,''); %dummy text handles to label agents 
            labelHandleNl=text(axeFig,0,0,''); 
            
            % Set some basic axes properties for axeFig.
            % NOTE: the following axeFig.Porperties need to be set after declaring all graphic objects children to axeFig since most graphic functions (such as "plot" or "rectangle") set
            % new values of axe.Properties and override any previous attempt to fix them...MATLAB is strange.
            axis equal % specifies the x:y scale of the axeFig axes; it's equivalent to setting axeFig.DataAspectRatio=[1 1 1] along with other more technical properties of axeFig
            axis([0 Lx 0 Ly]) % adjusts the graphic limits; acts the same as setting the XLim, YLim and ZLim properties of the axeFig axes.
            axis off %making all gridlines invisible: acts the same as setting axeFig.Visible='off'
            
        case 'No'
            disp_plot=false;
    end
end

%% Allocating memory for metrics and dynamical quantities to store.

% 0) We summarize the state of every agent in a matrix of size N_total x 4 :
% --> column #1 is NaN if agent NEVER entered the hallway, is 1 if agent is in the
% hallway by the time the simualtion ends and has value 0 if it entered
% and eventually exited via the rightful door.
% --> column #2 indicates the time iteration k at which the agent
% entered the hallway (NaN if the agent never entered)
% --> column #3 total amount of time iterations the agent spent in the
% hallway (0 if the agent never entered)
% --> column #4 is +1 if agent is part of the population moving RIGHT
% and is -1 if the agent is part of the population moving LEFT (0 if
% the agent never entered).
% NOTE: at the end of the script, columns #1 and #4 will be transformed into categorical arrays before saving all variables in the "hallwayData.mat" file.
agent_state=[NaN(N_total,1) NaN(N_total,1) zeros(N_total,1) NaN(N_total,1)];

% 1) Combinatorial/connectivity structures. i.e. graphs (need to be stored in cell arrays)
DT_t =cell(kMax,1); % Delaunay triangulation of the overall population (including the ghost copies used to efficiently enforce the boudness of Voronoi regions)
cuttedForestNr_t=cell(kMax,1); cuttedForestNl_t=cell(kMax,1); % Queueing structures of both populations (i.e. the weighted graphs that we obtain by our recursive reduction method of the minimal weight spanning forest (MWSF)).
QueuesLengthNr_t=cell(kMax,1); QueuesLengthNl_t=cell(kMax,1); 

% 2) Vector quantities (need to be stored in cell arrays)
Area_t=cell(kMax,1); % area of Voronoi cells
Vx_t=cell(kMax,1); Vy_t=cell(kMax,1);% x and y components of current velocities (i.e. at time t=k)
Vx_old_t=cell(kMax,1); Vy_old_t=cell(kMax,1); % x and y components of UNITARY old velocities (i.e. at time t=(k-1)dt )

% 3) Sacalar quantities (can be stored in numeric arrays)
AvgVel_t=zeros(kMax,1); % average and standard deviation of norms of current velocities
StdVel_t=zeros(kMax,1);

Nr_t=NaN(kMax,1); Nl_t=Nr_t;  %number of agents moving RIGHT and LEFT at each time step.

LengthInterface_t=Nr_t; %length of Voronoi interface between the two populations.

Polarization_t=Nr_t; % Polarization observables
PolarizationNr_t=Nr_t;
PolarizationNl_t=Nr_t;

Energy_t=Nr_t; %Voronoi energy, a.k.a clustering observables/metric
EnergyNr_t=Nr_t;
EnergyNl_t=Nr_t;

MomentAngAbs_t=Nr_t; % Angular and Absolute Angular momentum observables
MomentAng_t=Nr_t;

% Metrics pertinent to the observed queues and that are calculated over
% the cutted forests (i.e. the queueing formations).
QueueingMetricNr_t=Nr_t; QueueingMetricNl_t=Nr_t; % Our global queueuing metric defined over the cutted forests
LeafRatioMWSF_Nr_t=Nr_t;LeafRatioMWSF_Nl_t=Nr_t; %measure of ramification of the MWSF (the more ramified the less queueing is observed)  
LeafRatioCF_Nr_t=Nr_t;LeafRatioCF_Nl_t=Nr_t; %measure of ramification of the cutted forest, i.e. of the final queuing structures
MaxCuttingWeightNr_t=Nr_t; MaxCuttingWeightNl_t=Nr_t; % largest weight of the final cutted forest.


%% Time iterations
for k=1:kMax
    
    %% I) Stochastic process for entering the hallway
    % I.0) Initial insertion of agents in the hallway
    if k==1 % By default we insert one agent at random at each door at the beginning of the simulation.
        x=[tol_exit;Lx-tol_exit]; y=Ly*[rand(1);rand(1)];
        idxAgent_right=1; idxAgent_left=2;
        Nr_enter=1; Nl_enter=1;
        agentID_Nr=1; agentID_Nl=2; %we label the ID of the two agents first entering the hallway
        
        %Random angles of incidence: the one moving to the right has an incidence
        % angle in (-pi/2,+pi/2), while the one moving to the left has an angle in (+pi/2,3*pi/2)
        angl=pi*rand(2,1);
        vx_old=[cos(angl(1)-pi/2);cos(angl(2)+pi/2)]; vy_old=[sin(angl(1)-pi/2);sin(angl(2)+pi/2)];
    end
    
    %Update of labels (indices) of each moving group as well as the number of agents in it.
    idx_agent=[idxAgent_right; idxAgent_left]; %index of agents that are still inside the hallway.
    Nr=length(idxAgent_right); Nl=length(idxAgent_left); N=Nr+Nl; %counting the number of agents still in the hallway that are moving right/left
    
    if N<=1, break, end % break the time loop if only one agent remains in the hallway
    
    % To keep bringing agents inside the hallway at time t=k*dt via our stochastic entering process we need to find
    % the closest agent to every point y belonging to the entering set;
    % i.e. we need to solve argmin_i ||y-x_i|| as a function of y
    % for all points y of the entering set (i.e. the doors of the
    % hallway).
    % To make this search tractable and more time efficient:
    % i) we discretize the doors using the points in Y.
    % ii) for each point y of the discretized doors Y, we restrict our
    % search to candidates x_i whose Voronoi cell "can see" the
    % point y (i.e. y is within the Voronoi cell V_i).
        
    if Nr_enter+Nl_enter<N_total % We keep inserting agents provided there are some still in "reserve"
        
        % I.1) Preliminary quantities need to build the probability
        % density functions (PDF) of the entering process at each door.
        [~,DT,~,~,idx_bdry]=voronoi_rectangle(x,y,Omega);
        N_bdry=length(idx_bdry); %number of agents having a boundary agents
        idx_lateralcopies=N+4+[(N_bdry+1:2*N_bdry)' (3*N_bdry+1:4*N_bdry)']; % index of the lateral copies in terms of the ordering of x4 and y4 (in order to find the Voro cells touching the doors); see voronoi_rectangle.m file
        
        idx_touches_rightdoor=[]; idx_touches_leftdoor=[]; %initialize indices of agents whose Voronoi region touches the right and left door respectively.
        tempList= [DT(:,[1 2]); DT(:,[2 3]); DT(:,[3 1])]; %Concatenated idices of vertices that form each edge of the triangulation
        EdgeList=unique(tempList,'rows'); % contains 2 columns: for the starting-vertex and the end-vertex of each delaunay edge.
        
        % we find the agents whose Voronoi cell intersects the door:
        % i) we can restrict our search to boundary agents i \in idx_bdry.
        % ii) agent i touches the right door if it is a
        %    direct delaunay neighbour to its reflected copy to the
        %    right, the index of the copy is idx_lateralcopies(...,1)).
        % iii) similarly, agent i touches the left door iff it is in direct
        %   contact with its own reflected copy to the left. the index of the copy is idx_lateralcopies(...,2)).
        for i=1:N
            neigh=EdgeList(EdgeList(:,1)==i,2); %for each generator i, the list neigh contains the indices of Voronoi neighbors (including the ghost/copies created in voronoi_rectangle.m)
            idx=find(idx_bdry==i); % index of agent i w.r.t the ordering of idx_bdry (it's empty if agent i is not on the boundary)
            if ~isempty(idx)
                %right door
                intersection=intersect(neigh,idx_lateralcopies(idx,1));
                if ~isempty(intersection)
                    idx_touches_rightdoor=[idx_touches_rightdoor; i]; %#ok<*AGROW>
                end
                %left door
                intersection=intersect(neigh,idx_lateralcopies(idx,2));
                if ~isempty(intersection)
                    idx_touches_leftdoor=[idx_touches_leftdoor; i];
                end
            end
        end
        
        % Extracting info of agents moving right/left: these two lines are a technicallity to allow adding agents over time without afecting the
        % current implementation of the dynamics
        xr=x(1:Nr); yr=y(1:Nr); xl=x(Nr+1:end); yl=y(Nr+1:end);
        vxr_old=vx_old(1:Nr); vyr_old=vy_old(1:Nr); vxl_old=vx_old(Nr+1:end); vyl_old=vy_old(Nr+1:end);
        
        % I.2) Insertion of angents on the LEFT door according to the prescribed PDF
        idx_candidates_left_entry=idx_touches_leftdoor;
        while Nr_enter<Nr_total
            PDF=zeros(size(Y)); % initializing the discretized Probability Density Function (PDF) of the entering random variable.
            distY=sqrt((zeros(res_PDF-1,1)-x(idx_candidates_left_entry)').^2 +(Y-y(idx_candidates_left_entry)').^2); % matrix of size (res_PDF)x(idx_candidates_left_entry) containing the distance between each discretized door point and every agent in the list of candidates.
            [mdistY,idx_mdistY]=min(distY,[],2); % idx_mdistY is an index set w.r.t idx_candidates_left_entry
            
            idxY_supercritical=find(mdistY>=Lc); %indices w.r.t the vector Y (discretized left door) that are in the supercritical case, i.e points of Y whose minimal distance to an agent is bigger than Lc
            if isempty(idxY_supercritical)  % when idxY_supercritical is an empty array we are in the case where the boundary has been saturated and thus we move on to the rest of the time iteration.
                break
            end
            
            % As opposed to the old version v1.4, in v1.5 we have removed
            % the PDF's dependency on the velocity's angle of the incoming agents i
            % already in the hallway whose V_i touch the left door.
            PDF(idxY_supercritical)=1-transitionFun(mdistY(idxY_supercritical),L,3); % PDF supported on the discretized left door
            
            xr(end+1)=tol_exit; x(end+1)=xr(end); %#ok<*SAGROW> %the x-coordinate of the agent who just entered is always on the LEFT door (we use tol_exit to obtain numerical stability w.r.t the Voronoi diagram of the cloud of points)
            yr(end+1)=randsample(Y,1,true,PDF); y(end+1)=yr(end); % generate y-coordinate of new agent being inserted on the LEFT door according the PDF we just computed
            
            % As opposed to the old version v1.4, in v1.5 agents now enter
            % with a UNIT velocity vector having random angle.
            angl=pi*(rand()-0.5); % random angle of incidence for people moving to the RIGHT is in (-pi/2,+pi/2);
            vxr_old(end+1)=cos(angl); vx_old(end+1)=vxr_old(end);
            vyr_old(end+1)=sin(angl); vy_old(end+1)= vyr_old(end);
            
            idx_candidates_left_entry(end+1)=N+1; %since the new agent we inserted is at the LEFT door, it automatically becomes a candidate to explore in the construction of the PDF for the next person to enter.
            
            agentID_Nr(end+1)=Nr_enter+Nl_enter+1; %attributing an ID to the agent we just inserted.
            
            Nr_enter=Nr_enter+1; %update on the number of agents inserted thus far
            N=N+1;    
        end
        
        % I.3) Insertion of angents on the RIGHT door according to the prescribed PDF
        idx_candidates_right_entry=idx_touches_rightdoor;
        while Nl_enter<Nl_total
            PDF=zeros(size(Y)); %discretized Probability Density Function (PDF) of the random variable.
            distY=sqrt((Lx*ones(res_PDF-1,1)-x(idx_candidates_right_entry)').^2 +(Y-y(idx_candidates_right_entry)').^2); % matrix of size (res_PDF)x(idx_candidates_right_entry) containing the distance between each discretized door point and every agent in the list of candidates.
            [mdistY,idx_mdistY]=min(distY,[],2); % idx_mdistY is and index set w.r.t idx_candidates_right_entry
            
            idxY_supercritical=find(mdistY>=Lc); %indices w.r.t the vector Y (discretized right door) that are in the supercritical case, i.e points of Y whose minimal distance to an agent is bigger than Lc
            
            if isempty(idxY_supercritical)  % when idxY_supercritical is an empty array we are in the case where the boundary has been saturated and thus we move on to the rest of the time iteration.
                break
            end
            
            % As opposed to the old version v1.4, in v1.5 we have removed
            % the PDF's dependency on the velocity's angle of the incoming agents i
            % already in the hallway whose V_i touch the right door.
            PDF(idxY_supercritical)=1-transitionFun(mdistY(idxY_supercritical),L,3); % PDF supported on the discretized right door
            
            xl(end+1)=Lx-tol_exit; x(end+1)=xl(end);
            yl(end+1)=randsample(Y,1,true,PDF); y(end+1)=yl(end); % generate y-coordinate of new agent moving RIGHT being inserted according the PDF we just computed
            
            % As opposed to the old version v1.4, in v1.5 agents now enter
            % with a UNIT velocity vector having random angle.
            angl=pi*(rand()+0.5); % random angle of incidence for people moving to the LEFT is in (+pi/2,+3*pi/2);
            vxl_old(end+1)=cos(angl); vx_old(end+1)=vxl_old(end);
            vyl_old(end+1)=sin(angl); vy_old(end+1)= vyl_old(end);
            
            idx_candidates_right_entry(end+1)=N+1; %since the new agent we inserted is at the RIGHT door, it automatically becomes a candidate to explore in the construction of the PDF for the next person to enter.
            
            agentID_Nl(end+1)=Nr_enter+Nl_enter+1; %attributing an ID to the agent we just inserted.
            
            Nl_enter=Nl_enter+1;  %update on the number of agents inserted thus far
            N=N+1;
        end
        % I.4) Technical issue: if xr, xl, yr and yl are row vectors we transform them
        % into column vectors...same goes for the old velocities
        xr=xr(:); xl=xl(:); yr=yr(:); yl=yl(:); vxr_old=vxr_old(:); vxl_old=vxl_old(:); vyr_old=vyr_old(:); vyl_old=vyl_old(:);
        
        % Technical problem: if xr, xl, yr and yl are row vectors we transform them
        % into column vectors...same goes for the old velocities
        xr=xr(:); xl=xl(:); yr=yr(:); yl=yl(:);
        vxr_old=vxr_old(:); vxl_old=vxl_old(:); vyr_old=vyr_old(:); vyl_old=vyl_old(:);
        
        % We reorganize the order of our positions and old velocities so that all
        % agents moving to the right come first in the ordering and the agents
        % moving left come in last.
        x=[xr; xl]; y=[yr; yl];
        vx_old=[vxr_old; vxl_old]; vy_old=[vyr_old; vyl_old];  
        agentID=[agentID_Nr(:); agentID_Nl(:)];
        
        % Recompute useful quantities.
        Nr=length(xr); Nl=length(xl); N=Nr+Nl;
        idxAgent_right=(1:Nr)'; idxAgent_left=(Nr+1:N)';
        %idx_agent=[idx_agent_right; idx_agent_left]; %index of agents that are still inside the hallways
        
    end
    
    % I.5) Write down agent_state information after all new entries into the
    %hallway at time k.
    agentIDEnter=agentID(isnan(agent_state(agentID,2)));%among all agent actives, which are the ones that just entered this iteration? answer: the ones whose second column does not yet bear the time stamp of their entrance.
    agent_state(agentIDEnter,1)=1; %indicate that agent has entered the hallway: "NaN" if it never enters, "1" when currently in the hallway, "0" if it entered and already exited.
    agent_state(agentIDEnter,2)=k; %indicate at what time iteration the agent has entered the hallway
    agent_state(agentID,3)=agent_state(agentID,3)+1; %add 1 to the counter of time iterations each agent spends inside the hallway.
    agent_state(agentID_Nr,4)=+1; agent_state(agentID_Nl,4)=-1; %indicate if agent is moving right "+1" OR left "-1".
    
    Nr_t(k)=Nr; Nl_t(k)=Nl; %store the number of agents in each subpopulation at time k.
    Vx_old_t{k}=vx_old; Vy_old_t{k}=vx_old;%store the old velocity vectors (i.e. velocity vectors at time t=k-1 that will be used throughout the iteration t=k).
    
    % IMPORTANT: at this stage all newly created agents at iteration k have
    % been inserted and:
    % i) their position coordinates are included in the [x,y] variables
    % ii) their (unitary) OLD velocities have been included into the [vx_old,vx_old] variables
    % iii) their corresponding entry in the "agent_state" matrix has been updated.

    %% II) Delaunay and Voronoi diagram
    [~,DT,V,C,idx_bdry]=voronoi_rectangle(x,y,Omega); %Compute Delaunay topology and Voronoi diagram of all generators present in the hallway at iteration k.
    x4=DT.Points(:,1); y4=DT.Points(:,2); % full set of generators with ghost copies (see voronoi_rectangle.m file)
    N_bdry=length(idx_bdry); %number of boundary agents (i.e. agents in the hallway whose Voronoi region intersects the at least one of thhe hallway boundaries)
    DT_t{k}=DT; %store the Delaunay topology of all agents at time t=k*dt
    
    %% III) Topological field of view
    Neigh_Ui=cell(N,1); % the Union of Neigh with i itself, i.e in set notation {Neigh{i}}U{i}
    idx_topoField=cell(N,1); % idx_topoField{i} will be filled with all the agents that are in the p-Voronoi neighborhood of agent i (including i itself)
    
    AdjMat=sparse(N,N); % adjacency matrix of the Delaunay graph
    powerAdjMat=speye(N); % initializing the p power of the adjacency matrix, i.e. powerAdjMat=AdjMat^p
    cumPowerAdjMat=sparse(N,N); %initializing the cumulative sum of powers of the adjacency matrix, i.e. cumPowerAdjMat=AdjMat + AdjMat^2 + AdjMat^3 +...+ AdjMat^p
    
    idx_lateralcopies=N+4+[(N_bdry+1:2*N_bdry)' (3*N_bdry+1:4*N_bdry)']; % index of the lateral copies in terms of the ordering of x4 and y4 (in order to find the Voronoi cells touching the doors)
    idx_touches_rightdoor=[];idx_touches_leftdoor=[]; %initialize array of agents whose Voronoi region touches each door respectively
    
    tempList= [DT(:,[1 2]); DT(:,[2 3]); DT(:,[3 1])]; %Concatenated idices of vertices that form each edge of the triangulation
    EdgeList= unique(tempList,'rows'); % contains 2 columns: for the starting-vertex and the end-vertex of each delaunay edge.
    
    for i=1:N
        neigh=EdgeList(EdgeList(:,1)==i,2); %for each generator i, the list neigh contains the indices of Voronoi neighbors (including the ghost/copies)
        
        % We find the agents whose Voro cell intersects the door.
        % a) we can restrict our search to boundary agents i \in idx_bdry.
        % b) agent i touches the right door if it is a
        %    direct delaunay neighbour with its reflected copy to the
        %    right, the index of the copy is idx_lateralcopies(...,1)).
        % c) similarly, agent i touches the left door iff it is in direct
        %   contact with its own reflected copy to the left.the index of the copy is idx_lateralcopies(...,2)).
        
        idx=find(idx_bdry==i); % index of agent i w.r.t the ordering of idx_bdry (it's empty if agent i is not on the boundary)
        if ~isempty(idx)
            %right door
            intersection=intersect(neigh,idx_lateralcopies(idx,1));
            if ~isempty(intersection)
                idx_touches_rightdoor=[idx_touches_rightdoor; i];
            end
            %left door
            intersection=intersect(neigh,idx_lateralcopies(idx,2));
            if ~isempty(intersection)
                idx_touches_leftdoor=[idx_touches_leftdoor; i];
            end
        end
        
        neigh=neigh(neigh<=N); %we only care about the neighbours that are among the original set of N agents. (at this point agent i is not included in its own set of neighbors)
        Neigh_Ui{i}=[i; neigh]; %storing neighbors WITH agent i itself.
        AdjMat(i,neigh)=ones(1,length(neigh)); %#ok<SPRIX>
    end
    
    % Artifact to make AdjMat symmetric in case the ordering from EdgeList
    % skrew things up and we have an (i,j) dependency but not (j,i) or vice
    % versa.
    % Can only occur due to floating point error in the use of the built-in matlab library for Delaunay triangulations and never
    % arises in practice.
    symmDiff=AdjMat-AdjMat';
    if nnz(symmDiff)~=0
        warning('AdjMat was forcibly made symmetric')
        AdjMat(symmDiff~=0)=1;%#ok<SPRIX>
    end
    
    % Finally, we fill idx_topoField with the indices of agents in the topological fields
    if p==1
        idx_topoField=Neigh_Ui;
    else
        for i=1:p % if p>=2 then the topological field of view is obtained via with cumulative sum of powers of the adjacency matrix
            powerAdjMat=powerAdjMat*AdjMat;
            cumPowerAdjMat=cumPowerAdjMat+powerAdjMat;
        end
        
        for i=1:N
            idx_topoField{i}=find(cumPowerAdjMat(i,:)~=0);
            %NOTE: idx_topoField{i} contains i itself since AdjMat^2 has +1 in
            %the diagonals and thus cumPowerAdjMat has always nonzero diagonals.
        end
    end
    
    %% IV) Homming component h (i.e. targeting component)
    hy=zeros(N,1); %in the scenario of the simple hallway the y-component of h is always zero.
    if q==Inf % In case every agent is always able to "see" its target set
        hx=ones(N,1); hx(idxAgent_left)=-ones(Nl,1); % hx is 1 for agents moving right and -1 for agents moving left
        % we don't need the y-component of h since it's always 0 for this particular shape of hallway.
            
    else % In case agents can only see their target set (desired exit) if it falls whithin a topological distance of "q"
        hx=zeros(N,1); % initializing x-coordinate of homing component
        idx_hx_right=[]; idx_hx_left=[];% initializing array for the indices of agents effectively moving to the RIGHT and LEFT
        
        for i=1:Nr % the RIGHT exit is within a topo distance of q from agent i only if there is some agent j touching the RIGHT door that is in idx_topoField{i}
            if ~isempty(intersect(idx_touches_rightdoor,idx_topoField{i}))
                idx_hx_right=[idx_hx_right; i];
            end
        end
        hx(idx_hx_right)=+1; %updating x-coordinate of homing component for agents effectively moving RIGHT
      
        for i=Nr+1:N % the LEFT exit is within a topo distance of q from agent i only if there is some agent j touching the LEFT door that is in idx_topoField{i}
            if ~isempty(intersect(idx_touches_leftdoor,idx_topoField{i}))
                idx_hx_left=[idx_hx_left; i];
            end
        end
        hx(idx_hx_left)=-1; %updating x-coordinate of homing component for agents effectively moving LEFT
    end
    
    %% V) Alignment component
    ax=zeros(N,1); ay=zeros(N,1);% initializing arrays for final alignment component
    for i=1:N
        %%%%%% NOTE: vi_old and vj_old need to be already normalized%%%%%
        temp=idx_topoField{i}; %temporary variable
        idx_topoField2=temp(temp~=i); %we create a second idx_topoField index set but this time we remove agent i itself (for the displacement \tilda{a} the agent i will NOT average over its own velocity).
        ntopoField2=length(idx_topoField2); % number of agents in p-topoField of agent i (again...we don't count agent i itself).
        
        Vi_old=repmat([vx_old(i);vy_old(i)],[1 ntopoField2]); % old unit velocity vector (i.e at step t=(k-1)*dt) of agent i that is repeated/concatenated ntopoField times.
        Vj_old=[vx_old(idx_topoField2)'; vy_old(idx_topoField2)']; %old unit velocity vectors of all agents within the p-topological field of view of agent i (once more, without considering i itlsef).
        dotproduct=dot(Vi_old,Vj_old); %dot product of unit velocity vectors at time t=(n-1)*dt.
        
        theta_ij=acos(dotproduct); % angles in [0,pi] of the old vi with every old vj that lands in its p-topological field of vision.
        g_ij=transitionFun(theta_ij,pi,3); %weights for the average of the velocities (obtained via a selected transition function);
        
        alignment=Vj_old*g_ij'/H_p; % 2x1 vector of the SCALED alignment component of displacement of agent i.
        ax(i)=alignment(1); ay(i)=alignment(2);%extract x and y coordinates separately of the final alignment component of displacement.
    end
    
    %% VI) Repulsion component r
    % To include repulsion from the top and bottom boundaries of the hallway we
    % use the ghost copies produced by the voronoi_rectangle.m function
    
    idx_touches_rightdoor_right=idx_touches_rightdoor(idx_touches_rightdoor<=Nr); %indices of agents wanting to move to the RIGHT that touch the right door
    idx_touches_leftdoor_left=idx_touches_leftdoor(idx_touches_leftdoor>=(Nr+1)); %indices of agents wanting to move to the LEFT that touch the left door
    
    idx_repul=zeros(N,1);  % memory allocation; will contain the indices of the original agents or ghost copies that we want to move away from.
    delta=zeros(N,1); % memory allocation; will contain distance between agent i and the neighbor or ghost copy it will move away from.
    idx_inside=(1:N)'; idx_inside(idx_bdry)=[]; %indices of agents whose Voronoi region does not touch any boundary of the hallway.
    
    %VI.1) closest original agent
    [dist_NearNeigh,idx_NearNeigh] = min((x-x').^2+(y-y').^2+diag(inf(N,1)),[],2); % Surprisingly, for N under a couple thousands, it's less time consuming to brute-force distance calculations than to go over the already computed Delaunay topology
    dist_NearNeigh=sqrt(dist_NearNeigh); %distances to nearest neighbor
    
    % VI.2) for agents whose Voronoi region does not touch the boundary (i.e those
    % with index in idx_inside) it is immediate that they will move away from
    % their closest neighbor and not away from a ghost copy.
    idx_repul(idx_inside)=idx_NearNeigh(idx_inside);
    delta(idx_inside)=dist_NearNeigh(idx_inside);
    
    % VI.3) for agents whose Voronoi cell touches the boundary (idx_bdry) we
    % need to compare the distance to their closest neighbour (original agent) against the
    % half distance to their closest ghost neighbor...this will tell us if
    % we need to move away from an agent or away from a boundary piece of the hallway.
    
    % VI.3.0) At first we make each boundary agent move away from its closest ghost copy.
    idx_medialcopies=N+4+[(1:N_bdry)' (2*N_bdry+1:3*N_bdry)']; %index of the bottom and top copies we made of the boundary agents (in terms of the ordering of x4 and y4)
    idx_copies=[idx_medialcopies idx_lateralcopies]; % array of size N_bdry x 4: each row contains the indices (indexed w.r.t [x4 y4]) of the 4 copies we made of the boundary agents.
    
    X_bdry=repmat(x(idx_bdry),[1 4]); Y_bdry=repmat(y(idx_bdry),[1 4]); % columns are coordinates of original agents and are repeated% To illustrate this execute: [idx_bdry idx_copies]
    X_copies=x4(idx_copies); Y_copies=y4(idx_copies); % x and y coordinates of each of the 4 copies ghost copies of each boundary agent
    
    dist_copies=sqrt((X_bdry-X_copies).^2+(Y_bdry-Y_copies).^2); %brute force method to calculate the distance between boundary agents and all 4 of their respective ghost copies
    [sorted_dist_copies,idx_order]=sort(dist_copies,2); %increasingly sorted distance between a boundary agent and each of the 4 copies what what created of him.
    
    sorted_dist_copies=sorted_dist_copies/2; %We only care about the half distance (this is the distance to the boundary).
    
    idx_linear1=sub2ind(size(idx_copies),(1:N_bdry)',idx_order(:,1)); % linear indices to recover FIRST column of the sorting, e.g dist_copies(idx_linear1)==sorted_dist_copies(:,1) is a pure 1-boolean vector of size N_bdry x 1
    idx_linear2=sub2ind(size(idx_copies),(1:N_bdry)',idx_order(:,2)); % linear indices to recover SECOND column of the sorting, e.g dist_copies(idx_linear2)==sorted_dist_copies(:,2) is a pure 1-boolean vector of size N_bdry x 1
    
    idx_closecopies=idx_copies(idx_linear1); %indices of the CLOSEST copy/ghost w.r.t x4 and y4
    idx_close2copies=idx_copies(idx_linear2); %index of the SECOND CLOSEST copy/ghost w.r.t x4 and y4
    
    idx_repul(idx_bdry)=idx_closecopies;
    delta(idx_bdry)=sorted_dist_copies(:,1);
    
    % VI.3.1) For the agents moving right and touching the right door we change the
    % repusion to their SECOND closest neight ONLY WHEN the closest neighbour is a
    % RIGHT lateral copy.
    % Similarly for agents moving left and touching the left door we change the
    % repulsion to their SECOND closest neight ONLY WHEN the closest neighbour is a
    % LEFT lateral copy
    
    idx_rightcopy_closest=idx_bdry(idx_order(:,1)==3); % indices of the agents whose closest copy is also their RIGHT copy
    idx_leftcopy_closest=idx_bdry(idx_order(:,1)==4); % indices of the agents whose closest copy is also their LEFT copy
    
    idx_right_modify=intersect(idx_touches_rightdoor_right,idx_rightcopy_closest); %indices of the agents touching the right door and moving to the right for which we'll modify the repulsion
    idx_left_modify=intersect(idx_touches_leftdoor_left,idx_leftcopy_closest); %indices of the agents touching the left door and moving to the left for which we'll modify the repulsion
    
    % WATCH OUT BELOW: THE ORDER OF THE ARGUMENTS OF "ismember" IS VERY IMPORTANT
    bool_wrtbdry_right=ismember(idx_bdry,idx_right_modify); % boolean indices of values in idx_right_modify but in terms of idx_bdry, e.g idx_bdry(bool_wrtbdry_right)==idx_right_modify is all 1-boolean
    bool_wrtbdry_left=ismember(idx_bdry,idx_left_modify); % boolean indices of values in idx_left_modify but in terms of idx_bdry, e.g idx_bdry(bool_wrtbdry_left)==idx_left_modify is all 1-boolean
    
    idx_repul(idx_right_modify)=idx_close2copies(bool_wrtbdry_right); % we replace the closest copy by the second closest clopy for the appropriate agents touching their exiting doors.
    idx_repul(idx_left_modify)=idx_close2copies(bool_wrtbdry_left);
    
    delta(idx_right_modify)=sorted_dist_copies(bool_wrtbdry_right,2); %we replace the distance to the closest copy to be the distance to the second closest copy for the appropriate agents touching their exiting doors.
    delta(idx_left_modify)=sorted_dist_copies(bool_wrtbdry_left,2);
    
    % VI.3.3) Finally, for the agents in idx_bdry we compare the distances to the closest agent and the half
    % distances to the copy (closest or second closest depending on what was found in VI.3.1)
    % NOTE: because of our earlier calculation in VI.3.1), the agents in idx_bdry
    % are repelled by the boundary (by chost/copy agents). The last step is to
    % see if they should be repelled by their closer agent instead and in this
    % case change the repulsion.
    
    [delta(idx_bdry),flag]=min([dist_NearNeigh(idx_bdry)  delta(idx_bdry)],[],2); %compute minimal distance between original agents and ghosts agents...the factor 1/2 is to account for the symmetry at the boudary.
    
    idx1=find(flag==1); % indicate that the neighbour we need to move away from is an ORIGINAL agent and not a ghost copy.
    idx_repul(idx_bdry(idx1))=idx_NearNeigh(idx_bdry(idx1));
    
    % VI.3.2) We compute unitary direction for repulsion r with discontinuity at 0 (here we
    % consider 1e-10 (roughly the order of precision allowed buy the Voronoi diagram)
    repulx=x-x4(idx_repul); repuly=y-y4(idx_repul);
    nor=sqrt(repulx.^2+repuly.^2); % Nx1 vector of the norms of the vectors r
    idx_normlz=find(nor>=1e-10); % finding indices for which we need to normalize and those for which we will not (we actually enforce rx=ry=0 for the latter case)
    
    rx=zeros(N,1); ry=zeros(N,1); %initializing normalized repulsion vectors (r^\hat in the paper).
    rx(idx_normlz)=repulx(idx_normlz)./nor(idx_normlz); ry(idx_normlz)=repuly(idx_normlz)./nor(idx_normlz); %normalization of the repulsion component
    
    sigma=transitionFun(delta,L,3); %coefficient sigma for the repulsion component of displacement, see transitionFun.m file for details on the decaying cut-off function.
    
    %% VII) Compute overall displacement vector (weighted sum of homing, repulsion and alignment)
    
    displacementx=((1-sigma).*hx+sigma.*rx+nu*ax)/(1+nu); %total displacement components
    displacementy=((1-sigma).*hy+sigma.*ry+nu*ay)/(1+nu);
    
    % Next we compute the unit length vector of displacement (normalizing it using
    % Heavyside function). These unitary vectors will be used throughout the
    % rest of iteration t=k*dt.
    nor=sqrt(displacementx.^2+displacementy.^2); % Nx1 vector of the norms of the displacement vectors
    idx_normlz=find(nor>=1e-10);% finding indices for which we need to normalize
    
    normlzd_displx=zeros(N,1); normlzd_disply=zeros(N,1); %initializing normalized displacement vectors
    normlzd_displx(idx_normlz)=displacementx(idx_normlz)./nor(idx_normlz); normlzd_disply(idx_normlz)=displacementy(idx_normlz)./nor(idx_normlz);
    
    %% VIII) Computation the individual density-velocity response proper to each agent and final velocity.
    % We also take advantage and compute some metrics in this section (mainly Voronoi energy and Voronoi area)
    
    % VIII.1) We find the agents whose direct (straight) line of site intersects the
    % correct boundary without "leaving" the Voronoi cell. These agents agents will have maximal speed of 1.
    idx_unit_speed=[]; %initialization for the index of the agents whose Voro cell intersects the desired door and the velocity vector points towards the door "without leaving the Voro cell"
    
    % VIII.1.RIGHT) We look among the agents moving RIGHT and touching the RIGHT door.
    for ii=1:length(idx_touches_rightdoor_right)
        i=idx_touches_rightdoor_right(ii);
        tt=(1+1e-9)*(Lx-x(i))/normlzd_displx(i); %value of the parameter for which the line with unit director vector [normlzd_displx(i), normlzd_disply(i)] emanating from [x(i),y(i)] intersects the RIGHT door
        if tt>0 && tt~=Inf % the case t==Inf is when vi=0, we leave that case for later
            [x_intercept,~] = polyxpoly([x(i);Lx+1e-9],[y(i); tt*normlzd_disply(i)+y(i)],V(C{i},1),V(C{i},2)); %x-coordinate of the intercept of the line with the RIGHT door
            if abs(x_intercept-Lx)<=1e-7 %if that x-coordinate is "numerically close" to the RIGHT door then we add this agent to the list of unit speed.
                idx_unit_speed=[idx_unit_speed; i];
            end
        end
    end
    % VIII.1.LEFT) We look among the agents moving LEFT and touching the LEFT door.
    for ii=1:length(idx_touches_leftdoor_left)
        i=idx_touches_leftdoor_left(ii);
        tt=(1+1e-9)*(-x(i))/normlzd_displx(i); %value of the parameter for which the line with unit director vector [vx_old(i), vy_old(i)] emanating from [x(i),y(i)] intersects the LEFT door
        if tt>0  && tt~=Inf %the case t==Inf is when vi=0, we leave that case for later
            [x_intercept,~] = polyxpoly([x(i);-1e-9],[y(i); tt*normlzd_disply(i)+y(i)],V(C{i},1),V(C{i},2)); %x-coordinate of the intercept of the line with the LEFT door
            if abs(x_intercept)<=1e-7  %if that x-coordinate is "numerically close" to the LEFT door then we add this agent to the list of unit speed.
                idx_unit_speed=[idx_unit_speed; i];
            end
        end
    end
    
    % VIII.2) For agents NOT having have maximal speed of 1 we compute their individual speeds
    % according to the fundamental diagram of "density vs speed"
    idx_notunit_speed=(1:N)'; idx_notunit_speed(idx_unit_speed)=[]; %gives the index of agents for which we need to compute the "perceived" density 
    angleVel=atan2(normlzd_disply,normlzd_displx); % array Nx1 containing the angles w.r.t horizontal of the CURRENT velocities of the N agents.
    
    Area=zeros(N,1); % area of each whole Voronoi region
    Ener=zeros(N,1); % Voronoi energy (clustering metric) of each Voronoi region
    AreaForward=Inf(length(idx_notunit_speed),1); %initialize the forward Voronoi area in front of every agent
    
    % VIII.2.1) For the agents NOT having unit speed we compute several quantities
    % that are required for the dynamics as well as ordering parameters and
    % metrics.
    for ii=1:length(idx_notunit_speed)
        i=idx_notunit_speed(ii);
        %There are two cases, when displacement component calculated in step step VII) is 0 and when it's not
        if ~ismember(i,idx_normlz) %case where the displacement component is 0 (rarely arises in practice)
            [Ener(i),Area(i)]=poly_area_energy(V(C{i},1),V(C{i},2));
            AreaForward(ii)=Area(i)/2;
        else %case where the calculated displacement component is different than 0
            sectorx=cos(angleVel(i))*sectorx_basis-sin(angleVel(i))*sectory_basis; % "sector" represents a large rectangle oriented/rotated counterclockwise by angle_vi(i) so that it points in the direction of the current velocity of agent i.
            sectory=sin(angleVel(i))*sectorx_basis+cos(angleVel(i))*sectory_basis;
            
            V_tilda=sutherlandHodgman([V(C{i},1) V(C{i},2)],[sectorx+x(i) sectory+y(i)]); %to obtain the "forward Voronoi cell" V_tilda we crop the whole Voronoi polygon against the large rectangle discribed by "sector". 
            AreaForward(ii)=poly_area(V_tilda(:,1),V_tilda(:,2)); % calculating the area of the "forward Voronoi cell"
            [Ener(i),Area(i)]=poly_area_energy(V(C{i},1),V(C{i},2)); %calculating area and energy of the WHOLE Voronoi cell (for metric purposes, it does not affect the dynamics)
        end
    end
    
    % VIII.2.2) For the agents having unit speed we only need to compute Energy,
    % and Area.
    % NOTE: in the case of agents with unit speed these quantities are actually
    % not required for the dynamics (we only compute them to have these metrics at
    % hand).
    for ii=1:length(idx_unit_speed)
        i=idx_unit_speed(ii);
        [Ener(i),Area(i)]=poly_area_energy(V(C{i},1),V(C{i},2));
    end
    
    % VIII.3.1) Density-velocity response based on our locally perceived density
    rho=ones(N,1);
    rho(idx_notunit_speed)=tanh(2*exp(1)*AreaForward/(pi*L^2));
   
    % VIII.3.2) Current velocity vector 
    velx=rho.*displacementx; vely=rho.*displacementy; % x and y coordinates of the N velocity vectors.
    
    nor_vel=sqrt(velx.^2+vely.^2); %vector with the 2-norms of the N velocity vectors.
    Vx_t{k}=velx; Vy_t{k}=vely; % store current velocity vector
    
    %% IX) Computing global metrics
    % IX.0) Preliminary calculations for current velocity (not normalized)
    shiftx=x-Lx/2; shifty=y-Ly/2; % position vectors w.r.t the center of the hallway
    nor_shift=sqrt(shiftx.^2+shifty.^2);
    crossprod=shiftx.*vely-shifty.*velx;
    
    % IX.1) Scaled Voronoi energies (clustering metrics).
    rE=N^2*rE_factor; % Normalization factor for Voronoi energy. NOTE: DON'T MOVE THE PLACE OF rE as it depends on N
    Energy_t(k)=mean(Ener)*rE; %Total Voronoi Energy
    EnergyNr_t(k)=rE*sum(Ener(idxAgent_right))/N;  %Voronoi Energy of people moving RIGHT
    EnergyNl_t(k)=rE*sum(Ener(idxAgent_left))/N; %Voronoi Energy of people moving LEFT
    
    % IX.2) Area 
    Area_t{k}=Area;
    
    % IX.3) Polarization for the whole group of agents as well as for the subgroup
    % going right and the one going left.
    Polarization_t(k)=sqrt(sum(velx)^2+sum(vely)^2)/sum(nor_vel); % Polarization= norm of the sum of velocity vectors divided by the sum of the norm of the velocity vectors
    PolarizationNr_t(k)=sqrt(sum(velx(idxAgent_right))^2+sum(vely(idxAgent_right))^2)/sum(nor_vel(idxAgent_right));
    PolarizationNl_t(k)=sqrt(sum(velx(idxAgent_left))^2+sum(vely(idxAgent_left))^2)/sum(nor_vel(idxAgent_left));
    
    % IX.4) Angular Momentum and Absolute Angular Momentum (measured w.r.t the center of the hallway)
    MomentAng_t(k)=abs(sum(crossprod))/(nor_vel'*nor_shift); % Angular Momentum=norm of sum of cross products divided by the dot product of the norms of velocities and relative positions
    MomentAngAbs_t(k)=sum(abs(crossprod))/(nor_vel'*nor_shift); % Absolute angular momentum=sum of norm of cross product divided by the dot product of norms of r and velocities.
    
    % IX.5) Average and standard deviation of the norm of the velocity vectors
    AvgVel_t(k)=mean(nor_vel);
    StdVel_t(k)=std(nor_vel);
    
    % IX.6) Queuing formations
    if computeQueueing % only compute queueing structure and metric if the user specifies it by setting computeQueueing=true
        % IX.6.RIGHT) We extract the connectivity of the population going RIGHT
        graphNr=graph(AdjMat(1:Nr,1:Nr)); %#ok<*UNRCH> % subgraph of the Delaunay triangulation resticted to the population moving RIGHT (creates graph data structure)
        [cuttedForestNr_t{k},temp,QueuesLengthNr_t{k}] = queuingStructureAndMetric(graphNr,x(1:Nr),y(1:Nr),normlzd_displx(1:Nr),normlzd_disply(1:Nr),'right');
         QueueingMetricNr_t(k)=N*temp;

        % IX.6.LEFT) We extract the connectivity of the population going LEFT
        graphNl=graph(AdjMat(Nr+1:N,Nr+1:N)); %create graph data structure
        [cuttedForestNl_t{k},temp,QueuesLengthNl_t{k}] = queuingStructureAndMetric(graphNl,x(Nr+1:N),y(Nr+1:N),normlzd_displx(Nr+1:N),normlzd_disply(Nr+1:N),'left');
        QueueingMetricNl_t(k)=N*temp; 
    end
    
    % IX.7) Length of Voronoi interface between the two populations
    % We find the pairs of agents who share a Voronoi edge and in order to avoid counting
    % the same pair twice we filter it to the condition j<i
    % NOTE: Surprisingly, computing the interface in the following way is up to 20 times
    % faster than computing it directly inside a loop running over i=1:N
    [ii,jj]=find(AdjMat~=0); % extract connectivity list of the Delaunay
    idx=find(jj<ii); %to remove duplications we enforce j<i
    ii=ii(idx); jj=jj(idx);
    
    idx=[find(jj<=Nr & ii>=Nr+1);find(ii<=Nr & jj>=Nr+1)]; %we only care about agents going to the right having neighbours going to the left and vice-versa.
    idx_agents_interface=[ii(idx) jj(idx)]; % list of pairs of agents going in the opposite direction from each other that share a Voronoi edge.
    
    sz_interface=size(idx_agents_interface,1); % size of the interface (i.e. the number of Voronoi edges forming the interface)
    vx_interface=zeros(3,sz_interface); vy_interface=zeros(3,sz_interface); % list of x&y coordinates of all the Voronoi vertices forming the interface
    length_interface=zeros(sz_interface,1); %allocate memory of array to store the length of the Voronoi interface between the pair of agents in idx_agents_interface
    
    for l=1:sz_interface %loop over each pair of agents [i j] in idx_agents_interface
        idx_intersect=intersect(C{idx_agents_interface(l,1)},C{idx_agents_interface(l,2)}); %finding the index of the two common Voronoi vertices between agents i  and j (indexed w.r.t the Voronoi vertex list V)
        Vx=V(idx_intersect,1); Vy=V(idx_intersect,2);
        length_interface(l)=sqrt( (Vx(1)-Vx(2))^2+ (Vy(1)-Vy(2))^2 );
        vx_interface(:,l)=[Vx; NaN]; vy_interface(:,l)=[Vy; NaN];  
    end
    vx_interface=vx_interface(:); vy_interface=vy_interface(:); % list of x and y coordinates of the Voronoi vertices making the interface between
    LengthInterface_t(k)=sum(length_interface); %total length of the Voronoi interface between the two populations
  
    %% X) GRAPHICS
    if disp_plot && k<=kMax %only produce graphics if user indicated so
        
        % COMMENT OR UNCOMMENT THE UPDATE OF THE DIFFERENT GRAPHIC
        % HANDLES 1) THROUGH 7) BELOW IF YOU WISH TO DISPLAY THEM IN THE FRAME OR NOT
        
        %IMPORTANT: All graphics are obtained by updating at every iteration the dummy graphic
        %handles that were initialized in the preamble. This is significantly faster than
        %having to erase and create new graphic handles for every object
        %at every iteration.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1) Title of the figure:
        % Update the third line of the graphics title with information
        % about the current state of the hallway
        titleFig.String{2}=['$ n_{r}=',num2str(Nr),'\quad n_{l}=',num2str(Nl),'\quad \#$iter $=',num2str(k),'$'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2) Position and velocity:
        % Update the quiver graphics object to represent the current
        % position and velocity of agents
        set(quivNr,'XData',x(1:Nr),'YData',y(1:Nr),'UData',normlzd_displx(1:Nr)*arrwSize,'VData',normlzd_disply(1:Nr)*arrwSize)
        set(quivNl,'XData',x(Nr+1:end),'YData',y(Nr+1:end),'UData',normlzd_displx(Nr+1:end)*arrwSize,'VData',normlzd_disply(Nr+1:end)*arrwSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3) Bi-directional Voronoi interface:
        % Update the line graphical object displaying the Voronoi interface
        % between the two populations
%         set(lineInterface,'XData',vx_interface,'YData',vy_interface)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4) Voronoi tessellation of the whole population: 
%         updated_voronoiPlt=voronoi(DT); set(updated_voronoiPlt,'Visible','off')
%         set(voronoiPlt,'XData',updated_voronoiPlt(2).XData,'YData',updated_voronoiPlt(2).YData);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 5) Delaunay triangulation of the whole population:
        % Update the graphics of the Delaunay triangulation (of the overall
        % population).
        % IMPORTANT: we do not plot Delaunay edges connecting agents
        % to ghost copies outside the hallway (nor edges between
        % adjacent ghost copies)
%         updated_delaunayPlt=triplot(DT,'Visible','off'); %the triplot function outputs a line graphics object that we make invisible in our axeFig axes.
%         xData=updated_delaunayPlt.XData; yData=updated_delaunayPlt.YData; %we extract the x&y coordinates of all the endpoints of that line object
%         idxInside=((xData>0 & xData<Lx) & (yData>0 & yData<Ly)) | isnan(xData); % we clean it by only retaining the points whose x&y coordinates are inside the hallway. NOTE: we also need to keep the NaN values marking breakpoints for the line object
%         set(delaunayPlt,'XData',xData(idxInside),'YData',yData(idxInside)); %update the triplt line object declared in the preamble.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 6) Queueing formations:
        % Update the line graphics objects representing the queue
        % formations encoded in the connectivity of our Cutted Forests 
        if computeQueueing %ONLY UPDATE QUEUE LINE OBJECTS IF USER SPECIFIED TO COMPUTE THEM 
            endNodesCF_Nr=cuttedForestNr_t{k}.Edges.EndNodes;
            linEdge=reshape(endNodesCF_Nr.',[],1); lgth=length(linEdge); %we flatten the matrix of indices of endpoints of queuePltNr to obtain a stacked vector and compute its number of elements
            lineQueueNr=NaN((3*lgth-2)/2,2); idxNotNaN=(1:size(lineQueueNr,1))'; idxNotNaN(3*( 1:(lgth/2-1)))=[];
            lineQueueNr(idxNotNaN,:)=[x(linEdge) y(linEdge)]; %we create a list of [x,y] corrdinates of the edges of our queue structure and separate them by NaN values to indicate a break for the line graphical object
            set(queuePltNr,'XData',lineQueueNr(:,1),'YData',lineQueueNr(:,2)) %update the graphpltNr line graphical object in axeFig
            
            endNodesCF_Nl=cuttedForestNl_t{k}.Edges.EndNodes;
            linEdge=reshape(endNodesCF_Nl.'+Nr,[],1); lgth=length(linEdge); %we flatten the matrix of indices of endpoints of queuePltNl to obtain a stacked vector and compute its number of of elements
            lineQueueNl=NaN((3*lgth-2)/2,2); idxNotNaN=(1:size(lineQueueNl,1))'; idxNotNaN(3*( 1:(lgth/2-1)))=[];
            lineQueueNl(idxNotNaN,:)=[x(linEdge) y(linEdge)]; %we create a list of [x,y] corrdinates of the edges of our queue structure and separate them by NaN values to indicate a break for the line graphical object
            set(queuePltNl,'XData',lineQueueNl(:,1),'YData',lineQueueNl(:,2)) %update the graphpltNl line graphical object in axeFig
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 7) Agent labels:
        % Added to axeFig.Children as Text graphic handles
%         delete([labelHandleNr;labelHandleNl]); %We delete the current Text handles and then add the newly updated ones. There must be a better way to update the text objects in axeFig.Children but I couldn't come up with one.
%         labelHandleNr=text(x(1:Nr)+Lx/100,y(1:Nr),string(idxAgent_right),'Color',quivColorNr); % adding new Text handles with data from iteration k.
%         labelHandleNl=text(x(Nr+1:end)-Lx/60,y(Nr+1:end),string(idxAgent_left),'Color',quivColorNl);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        drawnow % callback to update all graphic objects in our axeFig axes
        
        if save_frames==true %only save the graphics if user indicated so
            tit=['000000' num2str(k)];  tit=tit(end-4:end); % number labeling of .png files with constant nnumber of characters
            print(tit,'-dpng') 
            movefile([tit '.png'],save_dirctry) % the .png image created is moved to the subfolder 'saved_frames'
        end 
    end
    
    %% XI) Compute new position of agents and enforce boudary conditions.
    
    x_new=x+dt*rho.*displacementx;
    y_new=y+dt*rho.*displacementy;
    
    vx_old=normlzd_displx; vy_old=normlzd_disply; % we update the old velocities: the unitary direction of displacement at time t=k*dt becomes v_old for t=(k+1)*dt.
    
    % Agents who are crossing the top or bottom boundaries retain their
    % current y-coordinate in the next step (while possibly changing their
    % x-coordinate)
    idx_stayput_top=find(y_new>=Ly-tol_exit); y_new(idx_stayput_top)=y(idx_stayput_top);
    idx_stayput_bottom=find(y_new<=tol_exit); y_new(idx_stayput_bottom)=y(idx_stayput_bottom);
    
    % Agents moving to the RIGHT cannot cross the LEFT boundary. Similarly,
    % agents moving LEFT cannot cross the RIGHT boundary
    idx_stayput_right=find(x_new(idxAgent_right)<=tol_exit); x_new(idxAgent_right(idx_stayput_right))=x(idxAgent_right(idx_stayput_right));
    idx_stayput_left=find(x_new(idxAgent_left)>=Lx-tol_exit); x_new(idxAgent_left(idx_stayput_left))=x(idxAgent_left(idx_stayput_left));

    % Enforce door conditions (exits)
    idx_exit_right=find(x_new(idxAgent_right)>=Lx-tol_exit); idx_exit_left=find(x_new(idxAgent_left)<=tol_exit); % indices w.r.t idx_agent_right and idx_agent_left that crossed the doors.
    idx_remove=[idxAgent_right(idx_exit_right);idxAgent_left(idx_exit_left)]; %indices w.r.t [x_new,y_new] of agents that crossed their desired door and are thus to be removed from the halllway.
    
    % Agents moving right that cross the right boundary of the hallway get
    % removed. Sames goes for agents moving left that already crossed the
    % left door.
    x_new(idx_remove)=[]; y_new(idx_remove)=[];
    vx_old(idx_remove)=[]; vy_old(idx_remove)=[];
    
    agent_state(agentID(idx_remove),1)=0; %once an agent has been absorbed we tag its state as "0" (i.e. its "Status" will correspond to "Exited" in DynamicalData.StateOfAgents)
    agentID(idx_remove)=[];
    
    % Update agents' positions.
    Nr=Nr-length(idx_exit_right); Nl=Nl-length(idx_exit_left); N=Nr+Nl;
    idxAgent_right=(1:Nr)'; idxAgent_left=(Nr+1:N)';
    agentID_Nr=agentID(1:Nr); agentID_Nl=agentID(1+Nr:N);
    
    x=x_new;
    y=y_new;
    
    if mod(k,5)==0 %we print the iteration number every few iterations. NOTE: printing every iteration number in the Command Winsow can increase computational time up to 2%
        disp(['#iter=',num2str(k)]) % display current iteration number in Command Window
    end
end %end loop over iterations

%% Save dynamical quantities and metrics
SimulationParameters=table(L,nu,Lc,kMax,p,q,Lx,Ly,dt,randSeed,'VariableNames',{'L','nu','Lc','IterMax','p','q','Lx','Ly','dt','randSeed'});

DynamicalData.NumberOfAgents=table(Nr_t,Nl_t,'VariableNames',{'Nr','Nl'});
temp=array2table(agent_state,'VariableNames',{'Status','EnteringTimeStamp','TimeInHallway','MovingToThe'});
temp.Status=categorical(temp.Status,[0 1 NaN],["exited","entered","hasNotEntered"]); %transform  1st column of agent_state into categorical array
temp.MovingToThe=categorical(temp.MovingToThe,[1 -1 NaN],["right","left","TBD"]); %transform  4th column of agent_state into categorical array
DynamicalData.StateOfAgents=temp;

DynamicalData.DelaunayTriangulations=DT_t; %store delaunay triangulation over iterations
DynamicalData.Velocities=table(Vx_t,Vy_t,Vx_old_t,Vy_old_t,'VariableNames',{'Velx','Vely','Velx_old','Vely_old'}); %store current and old velocities 
DynamicalData.VoronoiAreas=table(Area_t,'VariableNames',{'VoronoiAreas'}); % store the area of each of the N Voronoi cells at every iteration

Metrics.Polarizations=table(Polarization_t,PolarizationNr_t,PolarizationNl_t,'VariableNames',{'Polarization','PolarizationNr','PolarizationNl'}); %Polarizations
Metrics.AngularMomentums=table(MomentAng_t,MomentAngAbs_t,'VariableNames',{'AngularMomentum','AbsoluteAngularMomentum'}); % regular and absolute angular momentums
Metrics.Energies=table(Energy_t,EnergyNr_t,EnergyNl_t,'VariableNames',{'Energy','EnergyNr','EnergyNl'}); %Voronoi (clustering) energy
Metrics.VelocityDistribution=table(AvgVel_t,StdVel_t,'VariableNames',{'AverageVelocity','StdVelocity'}); % mean and std of the magnitude of velicities (after the 
Metrics.VoronoiInterfaceLength=table(LengthInterface_t,'VariableNames',{'InterfaceLength'});
Metrics.QueueStructures=table(cuttedForestNr_t,cuttedForestNl_t,QueuesLengthNr_t,QueuesLengthNl_t,'VariableNames',{'CuttedForestNr','CuttedForestNl','QueuesLengthNr','QueuesLengthNl'});
Metrics.QueueMetrics=table(QueueingMetricNr_t,QueueingMetricNl_t,'VariableNames',{'QueueingMetricNr','QueueingMetricNl'});
%%
save("hallwayData.mat",'SimulationParameters','DynamicalData','Metrics');

