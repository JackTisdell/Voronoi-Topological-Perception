%% Frame rendering of bidrectional flow on a the straight hallway using version 1.5 of the VTP model
%
% Script_Graphics.m loads the data of a SINGLE simulation contained in
% "hallwayData.mat" thought to be originally produced by the Script_dynamics.m script. 
% Some dynamical quantities found in "hallwayData.mat" are used here to produce
% frame-by-frame graphics of the configuration of agents in the hallway as time
% evolves. Some selected metrics are also displayed to showcase remarkable
% regimes, e.g. Voronoi interface length to measure "percolation"; our
% simple notion of queueing structures and queueing metrics to quantify line-up, etc.
%
% The "hallwayData.mat" file should contain the following Workspace
% variables. 
%
% - "SimulationParameters" which is a 1x10 table with the values of fixed parameters
% needed to uniquely distinguish a simulation produced by Script_Graphics.m:
%       * L= inter-personal distance preferred by agents
%       * Lc= inter-personal distance preferred by agents near the doors
%       (entrances and exits). Effectively replaces number of agents in the hallway as degree of freedom 
%       * nu= the alignment coefficient weighting alignment versus
%       repulsion and homing in the VTP model
%       * p and q are the sizes of topological fields of view for
%       alignment and homing respectively
%       * Lx and Ly are the dimensions of the hallway (Lx is length and Ly
%       is width)
%       * dt=the time step used in the simulation
%       * IterMax= number of time iterations the dynamics evolved over
%       * randSeed= pseudo-random number generator's seed specified for
%       reproductibility of the stochastic process used to insert agents in
%       the hallway.
%
% - "DynamicalData" is a 1x1 structure with the following fields:
%       * "NumberOfAgents" is table of size IterMax*2 with the number of agents
%       of each population at every time iteration
%       * "StateOfAgents" is a 4 column table recording the state of all
%       agents: if they entered the hallway, if they exited via their
%       respective door, time stamp of the iteration they entered and how long did
%       they spend inside the hallway, if their target is the right or left
%       door.
%       * "DelaunayTriangulations" contains all delaunayTriangulations
%       structures of the overall population present in the hallway at each time iteration.
%       NOTE: these delaunayTriangulations MATLAB structures contain the
%       "ghost" points used to enforce bounded Vornoi regions cropped against the
%       hallway; see the read-me description file for more details on the
%       this construction.
%       * "Velocities" contains the x&y components of the velocity vector
%       for each agent at each iteration. The vectors Velx_old and Vely_old
%       represent the x&y components of the UNITARY directions of the velocity
%       vectors at time t=(k-1)*dt that are used by the VTP model for the alignment component of displacement.
%       * "VoronoiAreas" contains the area of the Voronoi region of
%       every agent present in the hallway at iteration t=k*dt.
% 
% - "Metrics" is a 1x1 structure with fields:
%       * "Polarizations", measure the direction consensus of agents
%       * "AngularMomentums" contains both the original as well as the
%       absolute angular momentums used to measure rotational consensus
%       around an given point (here taken to be the center of the hallway).
%       * "Energies" contains values of Voronoi energy used to measure clustering of the hallway
%       * "VelocityDistribution" measures average and standard deviation
%       of the speed of all agents in the hallway at any given iteration
%       * "VoronoiInterfaceLength" is the total lenght of Voronoi
%       boundaries shared between the two populations. This is a measure of
%       percolation, i.e. how much a certain population is venturing into
%       the other.
%       * "QueueStructures" contains graph structures that correspond to
%       weighted graph forests, see the accompanying paper for details on their
%       construction.
%       % "QueueMetrics" contain several queueuing metrics aiming
%       to quantify how much queuing is present in a given configuration.
%       These metrics are calculated over the weighted forests in "QueueStructures"
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Ivan Gonzalez and Jack Tidsell at McGill University
%
% Last Update March 2023:
%
% should any bug or error arise please contact
% ivan.gonzalez@mail.mcgill.ca and jack.tisdell@mail.mcgill.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User specified settings for the frame-by-frame rendering

save_frames=true; % Do you want to save the frames? Set it to be logical true/false or numeric 1/0)
iterNumb="all"; % Either positive integer OR "all" OR "last". iterNumb sets the iteration number to be displayed, if you want all iterations 1:IterMax to be plotted set iterNumb='all'.
iterInterval=0; % this gives the radius of an interval centered around iterNumb for iterations to be displayed. DEFAULT value is 0

%% Load the data produced by Script_Dynamics and import simulation parameters
if exist('DynamicalData')==false || exist('Metrics')==false || exist('SimulationParameters')==false %#ok<EXIST> % If the necessary variables are not already in the Workspace we load them up
    disp("Data is being loaded into the Workspace, it may take a few seconds.")
    load("hallwayData.mat")
    disp("Data loaded successfully!")
end
simuParam=num2cell(SimulationParameters{1,:}); % MATLAB trick of using a cell array to distribute data to multiple variables without the need to access each field of the structure SimulationParameters separately.
[L,nu,Ls,kMax,p,q,Lx,Ly,~,~]=simuParam{:};

%% Validate user-defined inputs for iterations to display and if frames should be saved
[save_frames,save_dirctry]=doWeSaveFrames(save_frames); % Local function with question dialogues so that user corroborates if they want to save the frames produced
iterDisplay=validateDisplayIter(kMax,iterNumb,iterInterval); % Local function to validate the user-defined iterations to be displayed

%% Pre-process some of the imported data for easier plotting.

LengthInterface=Metrics.VoronoiInterfaceLength.InterfaceLength/Ly; %Normalized length of the bi-directional Voronoi interface (normalized w.r.t the width of the hallway)
[VelxNrmlzd,VelyNrmlzd]=myCellFunNormalization(DynamicalData.Velocities.Velx,DynamicalData.Velocities.Vely); %normalizing the velocity vectors (just to plot unit arrows in the quiver objects)

%% I) Prepare graphic objects and handles

% Create RGB color triplets using the rgb.m function (to plot agents and their velocity vectors).
addpath(genpath('rgbColors')) %add path to the rgb function in the rgbColors folder
quivColorNr=rgb('OrangeRed'); quivColorNl=rgb('DarkGreen'); %colors for plotting agents and their velocity vectors
forestColorNr=rgb('Orange'); forestColorNl=rgb('MediumSeaGreen'); %colors for plotting the queueing structures (forests).
arrwSize=Ly/40;% size of arrows in quiver object to represent velocity

% I.0)Create a figure and axes
fig=figure('Name','Bi-directional corridor: dynamics and metrics');clf;
fig.Units='centimeters'; fig.Position=[4.5 13.5 37.5 11.5]; % fix some of its properties
axeHallway=axes(fig); %creating an axes object such that axeFig.Parent=fig
axeHallway.Units='normalized'; axeHallway.Position=[0.01 0.01 0.980 0.99]; %  normalized units represent a proportion of the Parent object fig
axeHallway.NextPlot='add'; % equivalent to "hold on": in order to add new graphic objects as children of axeFig. To see the difference type axeFig.Children when toggling between axeFig.NextPlot='add' and axeFig.NextPlot='replace'.


% I.3) Initialize graphic handles of objects to plot
if q==Inf, chara='\infty'; else, chara=num2str(q); end %for the title of the figure.

titleFig=title(axeHallway,{'Regime (V)';['Bi-directional flow with parameter values $\nu=',num2str(nu),', L=',num2str(L/Ly),', L_s=',num2str(Ls/Ly),'$'];''},'FontSize',14,'Interpreter','latex'); %we leave the third line of the title as an empty character that will be updated at every iteration with important information of the current iteration
            
rectngl=rectangle('Position',[0 0 Lx Ly]); % draws the boundary of the hallway using a rectangle graphics object

%create dummy quiver objects for position and velocity of out agents: will be updated at every iteration
quivNr=quiver(axeHallway,-1,-1,0,0); set(quivNr,'Color',quivColorNr,'Marker','.','MarkerSize',Ly,'ShowArrowHead','off','AutoScale','off','LineWidth',1); % for population moving RIGHT
quivNl=quiver(axeHallway,-1,+1,0,0); set(quivNl,'Color',quivColorNl,'Marker','.','MarkerSize',Ly,'ShowArrowHead','off','AutoScale','off','LineWidth',1); % for population moving LEFT

% dummy line object to plot the queueing structures. NOTE: we chose to use a line graphics object rather that a GraphPlot object as the latter
% cannot have its edge-connectivity updated at every iteration in order to plot the changes.
queuePltNr=line(axeHallway,-1,-1,'Color',forestColorNr,'LineWidth',1);
queuePltNl=line(axeHallway,-1,-1,'Color',forestColorNl,'LineWidth',1);

% dummy line object to plot the bi-direction Voronoi interface
lineInterface=line(axeHallway,-1,-1,'Color','magenta','LineWidth',1.5);

% dummy line objects to plot the Delaunay triangulation and Voronoi
% tessellation of the whole population
delaunayPlt=line(axeHallway,-1,-1,'Color','blue','LineStyle',':');
voronoiPlt=line(axeHallway,-1,-1,'Color','black','LineStyle',':');

% text objects for the labels of agents
labelHandleNr=text(axeHallway,0,0,''); %dummy text handles to label agents
labelHandleNl=text(axeHallway,0,0,'');

% Set some basic axes properties for axeFig.
% NOTE: the following axe Porperties need to be set after declaring all graphic objects children to axeHallway since most graphic functions (such as "plot" or "rectangle") set
% new values of axe.Properties and override any previous attempt to fix them...MATLAB is strange.
axis equal % specifies the x:y scale of the axeHallway axes; it's equivalent to setting axeHallway.DataAspectRatio=[1 1 1] along with other more technical properties of axeHallway
axis([0 Lx 0 Ly]) % adjusts the graphic limits; acts the same as setting the XLim, YLim and ZLim properties of the axeHallway axes.
axis off %making all gridlines invisible: acts the same as setting axeHallway.Visible='off'

%% Time loop over the selected iteration numbers
for j=1:length(iterDisplay)
    iter=iterDisplay(j);
    iterArray=(1:iter); %do NOT change iterArray into a column vector, leave it as a ROW.
    
    %% III) Update the hallway configuration
    Nr=DynamicalData.NumberOfAgents{iter,1}; Nl=DynamicalData.NumberOfAgents{iter,2}; N=Nr+Nl;
    DT=DynamicalData.DelaunayTriangulations{iter}; % Delaunay triangulation of the current iteration
    x4=DT.Points(:,1); y4=DT.Points(:,2); %extracting position of agents from the Delaunay triangulation (including ghost copies of agents)
    dirx=VelxNrmlzd{iter}; diry=VelyNrmlzd{iter}; %direction of each agent (extracted out of the normalized velicity vectors)
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) Title of the figure:
    % Update the last line of the graphics title with information
    % about the state of the current iteration
    %--------------------------------------------------------------------------
    titleFig.String{3}=['$n_{r}=',num2str(Nr),'\quad n_{l}=',num2str(Nl),'\quad $iter\# $=',num2str(iter),'$'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) Position and velocity:
    % Update the quiver graphics object to represent the current
    % position and direction of agents
    %--------------------------------------------------------------------------
    set(quivNr,'XData',x4(1:Nr),'YData',y4(1:Nr),'UData',dirx(1:Nr)*arrwSize,'VData',diry(1:Nr)*arrwSize)
    set(quivNl,'XData',x4(Nr+1:N),'YData',y4(Nr+1:N),'UData',dirx(Nr+1:N)*arrwSize,'VData',diry(Nr+1:N)*arrwSize)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3) Bi-directional Voronoi interface:
    % Update the line graphical object displaying the Voronoi interface
    % between the two populations
    %--------------------------------------------------------------------------
%     [V,C]=voronoiDiagram(DT); %recover topology of Voronoi tessellation from DT
%     tempList=[DT(:,[1 2]); DT(:,[2 3]); DT(:,[3 1])]; EdgeList= unique(tempList,'rows'); %Concatenated idices of agents that form each edge of the triangulation
%     EdgeList=EdgeList(EdgeList(:,1)<=N,:); EdgeList=EdgeList(EdgeList(:,2)<=N,:); EdgeList=EdgeList(EdgeList(:,1)<EdgeList(:,2),:); % Clean up EdgeList: we remove Delaunay edges with at least one endpoint outside the hallway and remove duplicated edges.
%     idx=find( (EdgeList(:,2)<=Nr & EdgeList(:,1)>=Nr+1) | (EdgeList(:,1)<=Nr & EdgeList(:,2)>=Nr+1) ); % Further clean up: since we only care about agents whose Voronoi region are at the interface, we only keep edges whose endpoints are agents of the opposite population.
%     idx_agents_interface=EdgeList(idx,:); % index w.r.t x4&y4 of pairs of agents whose Voronoi boundary is included in the interface.
%     sz_interface=size(idx_agents_interface,1); vx_interface=NaN(3,sz_interface); vy_interface=NaN(3,sz_interface); %declaring NaN values is useful to declare breaks in the Line Graphics object
%     for l=1:sz_interface %loop over every Voronoi edge of the interface
%         idx_intersect=intersect(C{idx_agents_interface(l,1)},C{idx_agents_interface(l,2)}); %finding the index of the two common Voronoi vertices between agents i  and j (indexed w.r.t the Voronoi vertex list V)
%         vx_interface(1:2,l)=V(idx_intersect,1); vy_interface(1:2,l)=V(idx_intersect,2);  % 
%     end
%     set(lineInterface,'XData',vx_interface(:),'YData',vy_interface(:)) %Using NaN values is a clever trick to mark breaks in the Line Graphic Object
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4) Voronoi tessellation of the whole population:
    % NOTE: Sadly, it is ~20 faster to use the function voronoi(DT) to plot
    % the Voronoi tessellation than to work out a clever way to plot it
    % using the topology already extracted via the call [V,C]=voronoiDiagram(DT)
    %--------------------------------------------------------------------------
    updated_voronoiPlt=voronoi(DT); set(updated_voronoiPlt,'Visible','off') % function "voronoi" does not allow to set "Visible" propery directly from the call; needs to be set posteriorly
    set(voronoiPlt,'XData',updated_voronoiPlt(2).XData,'YData',updated_voronoiPlt(2).YData);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5) Delaunay triangulation of the whole population:
    % Update the graphics of the Delaunay triangulation (of the overall
    % population).
    % IMPORTANT: we do not plot Delaunay edges connecting agents
    % to ghost copies outside the hallway (nor edges between
    % adjacent ghost copies).
    %--------------------------------------------------------------------------
%     updated_delaunayPlt=triplot(DT,'Visible','off'); %the triplot function outputs a line graphics object that we make invisible in our axeFig axes.
%     xData=updated_delaunayPlt.XData; yData=updated_delaunayPlt.YData; %we extract the x&y coordinates of all the endpoints of that line object
%     idxInside=((xData>0 & xData<Lx) & (yData>0 & yData<Ly)) | isnan(xData); % we clean it by only retaining the points whose x&y coordinates are inside the hallway. NOTE: we also need to keep the NaN values marking breakpoints for the line object
%     set(delaunayPlt,'XData',xData(idxInside),'YData',yData(idxInside)); %update the triplt line object declared in the preamble.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6) Queueing formations:
    % Update the line graphics objects representing the queue
    % formations encoded in the connectivity of our Cutted Forests
    %--------------------------------------------------------------------------
%     endNodesCF_Nr=Metrics.QueueStructures.CuttedForestNr{iter}.Edges.EndNodes;
%     linEdge=reshape(endNodesCF_Nr.',[],1); lgth=length(linEdge); %we flatten the matrix of indices of endpoints of queuePltNr to obtain a stacked vector and compute its number of elements
%     lineQueueNr=NaN((3*lgth-2)/2,2); idxNotNaN=(1:size(lineQueueNr,1))'; idxNotNaN(3*( 1:(lgth/2-1)))=[];
%     lineQueueNr(idxNotNaN,:)=[x4(linEdge) y4(linEdge)]; %we create a list of [x,y] corrdinates of the edges of our queue structure and separate them by NaN values to indicate a break for the line graphical object
%     set(queuePltNr,'XData',lineQueueNr(:,1),'YData',lineQueueNr(:,2)) %update the graphpltNr line graphical object in axeFig
%     
%     endNodesCF_Nl=Metrics.QueueStructures.CuttedForestNl{iter}.Edges.EndNodes;
%     linEdge=reshape(endNodesCF_Nl.'+Nr,[],1); lgth=length(linEdge); %we flatten the matrix of indices of endpoints of queuePltNl to obtain a stacked vector and compute its number of of elements
%     lineQueueNl=NaN((3*lgth-2)/2,2); idxNotNaN=(1:size(lineQueueNl,1))'; idxNotNaN(3*( 1:(lgth/2-1)))=[];
%     lineQueueNl(idxNotNaN,:)=[x4(linEdge) y4(linEdge)]; %we create a list of [x,y] corrdinates of the edges of our queue structure and separate them by NaN values to indicate a break for the line graphical object
%     set(queuePltNl,'XData',lineQueueNl(:,1),'YData',lineQueueNl(:,2)) %update the graphpltNl line graphical object in axeFig
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 7) Agent labels added to axeFig.Children as Text graphic handles
    % WATCH OUT: Aents's labels are NOT the same as their ID. The label is
    % their index number referenced w.r.t the current list of N agents [x,y].
    % Whenever agents exit and/or enter, the label of the remaining agents in
    % the hallway changes. On the other hand, an agent's ID is a fixed number
    % that discribes it throughout the experiment.
    %--------------------------------------------------------------------------
%     delete([labelHandleNr;labelHandleNl]); %We delete the current Text handles and then add the newly updated ones. There must be a better way to update the text objects in axeFig.Children but I couldn't come up with one.
%     labelHandleNr=text(x4(1:Nr)+Lx/100,y4(1:Nr),string(1:Nr),'Color',quivColorNr); % adding new Text handles with data from iteration k.
%     labelHandleNl=text(x4(Nr+1:N)-Lx/60,y4(Nr+1:N),string(Nr+1:N),'Color',quivColorNl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    drawnow %update al graphic objects
    
    if save_frames==true %only save the graphics if user indicated so
        tit=['000000' num2str(iter)];  tit=tit(end-4:end);  % number labeling of .png files with constant nnumber of characters
        print(tit,'-dpng','-r120') % save to .png figure with specified resolution
        movefile([tit '.png'],save_dirctry) % the .png image created is moved to the subfolder 'saved_frames'
    end
    
end


%% Local functions
function varargout=myCellFunNormalization(u,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes two cell arrays u&v, each containing numeric vectors in their
% cells such that u{k} corresponds to the x-component of N vectors in 2D and
% similarly for v{k} for the y-component.
% IMPORTANT: N can vary with the index k (that's the whole point of the function)
% The function then normalizes all these N=N(k) vectors inside each of the
% k=1,...,K cells.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nrm=cellfun(@(var1,var2) hypot(var1,var2),u,v,'UniformOutput',false); %each cell of nrm contains numeric vectors with the norm of every 2D vector
    varargout{1}=normalizeVecCell(u,nrm); %entry-wise divide every element of u{k} and v{k} by nrm
    varargout{2}=normalizeVecCell(v,nrm);
    function nrmVec=normalizeVecCell(var,nor) %divide vector inside each cell by their respective norms
        nrmVec=cellfun(@(dummyVar,dummyNorm) dummyVar./dummyNorm,var,nor,'UniformOutput',false);
    end

end

function [save_frames,save_dirctry]=doWeSaveFrames(save_frames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function pops a dialogue box to make sure the user wants to save the
% rendered images (as .png files). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(islogical(save_frames) || ismember(save_frames,[0 1]))
    error('Please specify save_frames as a logical true/false or as a numeric 1/0 value')
end
save_dirctry=[]; %default directory
if save_frames==true
    answerQ2 = questdlg('Are you sure you wish to save the plotted frame? The running time can increase significantly?','Just making sure','Yes, save','No, run without saving','No, run without saving');
    if strcmp(answerQ2,'No, run without saving')
        save_frames=false;
    else
        if ~isfolder('saved_frames') %if there is no subfolder called 'saved_frames', create one
            mkdir saved_frames
        end
        save_dirctry=strcat(pwd,"/saved_frames"); % directory for saving the frames
    end
end
end

function iterDisplay=validateDisplayIter(kMax,iterNumb,iterInterval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function validates if the iterations to be displayed are compatible
% with the user defined parameters and with the data produced by Script_dynamics.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnumeric(iterInterval) && iterInterval>=0
    if mod(iterInterval,1)~=0
        iterInterval=ceil(iterInterval);%in case iterInterval is not an integer we round upwards
        disp(['The value you gave for iterInterval has been rounded up to iterInterval=',num2str(iterInterval)])
    end
else 
    disp('The value you gave for iterInterval is nonsensicle: the DEFAULT value iterInterval=0 will be used instead')
    iterInterval=0;
end

if isstring(iterNumb) || ischar(iterNumb)
    try % if iterNumb is given as text, we validate it against "all" OR "last".
        validStr=validatestring(iterNumb,["all" "last"]);
    catch
        error("The string or char you entered for iterNumb needs to be validated against ""all"" OR ""last"".")
    end
    if validStr=="all" 
        iterNumb=0;
        iterInterval=kMax;
    else %case where validStr="last"
        iterNumb=kMax;
    end
elseif isnumeric(iterNumb) && iterNumb>0 && mod(iterNumb,1)==0 %  if iterNumb is of numeric class, positive and is an integer
    if iterNumb>kMax
        disp('The value you gave for iterNumb is larger than kMax: the value iterNumb=kMax will be used instead')
        iterNumb=kMax;
    end
else
    error("Please set iterNumb to be a positive integer or a string/char that can be validated against ""all"" OR ""last"". ")
end

iterDisplay=(max(1,iterNumb-iterInterval):min(kMax,iterNumb+iterInterval)).'; % integer interval centered at iterNumb of radius iterInterval and cropped against [1:kMax];
end


