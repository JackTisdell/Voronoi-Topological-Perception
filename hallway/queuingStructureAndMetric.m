function [queueStructure,queueMetric,queueLengthDistri] = queuingStructureAndMetric(restrictDT,x,y,ux,uy,movingToThe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: This function constructs the queuing structure (graph)
% associated to a subpopulation in the bi-directional hallway and computes
% queuing data on it.
% Each connected component of the construction is interpreted as a distinc lane
% and the construction itslef satifsies the 4 postulates that a "reasonable"
% queuing structure should satisfy:
%  i) the structures should be a subgraph of the overall Delaunay
%  triangulation and only include edges linking two people from the same
%  subpopulation
% ii) all vertices (agents) of the structure have degree=1,2, i.e. its a
% non-ramified graph without singletons.
% iii) it's a forest, i.e. each connected component (lane) is acyclic.
% iv) to each edge we associate a weight measuring its contribution to
% "good" queing based on: the relative position of agents, the homing
% direction, the current velocity of both agents at the endpoints.
%
% NOTE: more details about the specific steps of the construction are found at
% the end of this file.
%--------------------------------------------------------------------------
% INPUT:
% 1) the delaunay triangulation RESTRICTED to one subpopulation. It needs
%   to be a MATLAB "graph" type variable.
% 2) the position [x,y] of the subpopulation in question.
    %IMPORTANT: ONLY input the [x,y] coordinates of the Nr/Nl agents of the
    %subpopulation. DO NOT input the [x,y] of all the N=Nr+Nl agents.
% 3) the UNITARY velocity vectors [ux,uy] of the the subpopulation in question
%   % IMPORTANT: ONLY input the [ux,uy] coordinates of the Nr/Nl agents of the
    % subpopulation.
% 4) movingToThe need to be set to 'right' for the subpopulation moving RIGHT 
%    and 'left'for the subpopulation moving to the LEFT

% OUTPUTS:
% 1) The queuing structure associated to the population. Each connected
% component of this graph is to be interpreted as a lane.
% 2) The queuing metric that measures que quality of a queuing formation
% acording to the four criteria:
%       - number of lanes (connected components of queueStructure)
%       - overall length of the lanes
%       - overall weight of each lane
%       - total number of agents within the subpopulation engaged in the
%         a lane.
% 3) Distribution of the topological length (i.e. number of edges) of lanes: 
%    e.g. if queueStructure has four lanes with corresponding lengths
%                lane 1 --> length=4
%                lane 2 --> length=7
%                lane 3 --> length=3
%                lane 4 --> length=7
%
%   then queueLengthDistri=[4 7 3 7];
%
%  NOTE: If there are NO lanes then queueLengthDistri=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Ivan Gonzalez on May 2023
% Last update: May 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Properties to fix depending on the supopulation.
switch movingToThe
    case 'right'
        sortOrder='ascend';
    case 'left'
        sortOrder='descend';
end

%% Calculate weights to associate to the restricted Delaunay Traingulation
idx_endpts=restrictDT.Edges.EndNodes; % pair of nodes at endpoints of UNDIRECTED edges indexed w.r.t [x,y].

if size(idx_endpts,1)~=1 %when there is only ONE edge in restrictDT, MATLAB messes up the size of the arrays and we need to treat that case separately
    [~,s_idx]=sort(x(idx_endpts),2,sortOrder); % sorting each row according to x-ccordinate of agents in the cc and extracting corresponding sorting index
    idx_endpts=idx_endpts(sub2ind(size(idx_endpts),repmat((1:size(idx_endpts,1))',1,2),s_idx)); % indices w.r.t [x,y] of starting and ending points of DIRECTED EDGES (first column is starting point and second column is ending point)
else
    [~,s_idx]=sort(x(idx_endpts),sortOrder);
    idx_endpts=idx_endpts(s_idx);
end
temp=[(x(idx_endpts(:,2))-x(idx_endpts(:,1)))'; (y(idx_endpts(:,2))-y(idx_endpts(:,1)))']; % vector coordinates of directed edges
nor=sqrt(temp(1,:).^2+temp(2,:).^2); idx_normlz_edge=find(nor>=1e-10); %norm of each Delaunay edge

dir_edge=zeros(2,size(idx_endpts,1)); %initialize vector coordinates of directed edges.
dir_edge(:,idx_normlz_edge)=temp(:,idx_normlz_edge)./nor(idx_normlz_edge); %we normalize each Delaunay directed edge if its numerically stable as we only need the unitary vectors.

dotprod1=dot([ux(idx_endpts(:,1))'; uy(idx_endpts(:,1))'],dir_edge); % dot product between directed edge and current velocity of agent at STARTING POINT of the directed edge
dotprod2=dot([ux(idx_endpts(:,2))'; uy(idx_endpts(:,2))'],dir_edge); % dot product between directed edge and current velocity of agent at ENDING POINT of the directed edge
weights=0.5*(acos(dotprod1)'+acos(dotprod2)'); %each edge weight is the average of the required angles measured at its endpoints.% NOTE: because of Matlab's weird indexing structure, applying the sorting index "si" to idx_endpts_xy is complicated

%IMPORTANT: MATLAB's acos(.) function is not numerically stable; it
%may yield complex weights. In that case we use the custom made
%arc_cosine.m function to properly calculate our graph weights
if ~isreal(weights)
    weights=0.5*(arc_cosine(dotprod1)'+arc_cosine(dotprod2)');
end

restrictDT.Edges.Weight=weights; % assign the new weights to edges.

%% Obtain the Minimum Weight Spanning Forest (MWSF) 

minWeighSpanForest=minspantree(restrictDT,'Type','forest'); %minimal weight spanning forest of restrictDT, i.e. the minimum spanning tree over each connected component of restrictDT

endnodesMWSF=minWeighSpanForest.Edges.EndNodes;
if ~isempty(endnodesMWSF) %there is a technical issue when endnodesMWSF is an empty array. In that case we need to manually make weightsMWSF an empty array as well.
    weightsMWSF=minWeighSpanForest.Edges.Weight;
else
    weightsMWSF=[];
end
[~,s_idx]=sort(weightsMWSF); % sorting all weights in the MWSF

%% Cutting the forest
% We start cutting the forest until in decreasing order of weights until only vertices with degree<=2 are left
cuttedForest=minWeighSpanForest; deg=degree(cuttedForest);
j=1;
while any(deg>=3)
    remEdge_idx=endnodesMWSF(s_idx(length(weightsMWSF)-j+1),:); % indices w.r.t [x,y] of both end points of the edge to remove (i.e. the edje with the jth largest weight)
    cuttedForest=rmedge(cuttedForest,remEdge_idx(1),remEdge_idx(2));
    deg=degree(cuttedForest);
    j=j+1;
end

%% Optional step: remove any connected component having only two vertices (i.e. only one edge)
[bins,binsizes]=conncomp(cuttedForest,'OutputForm','cell'); %computing the connected components (Trees) of the forest.
idx_filter_cnncomp=cell2mat(bins(binsizes==2)'); % we identify the connected components who have only two nodes (i.e. only one edge)
if ~isempty(idx_filter_cnncomp) %we remove any connected component that has only two nodes (one edge) since lines are considered to be made of 3 agents or more.
    cuttedForest=rmedge(cuttedForest,idx_filter_cnncomp(:,1),idx_filter_cnncomp(:,2));
end
%% Process data on the final cutted forest and compute metrics defined on it.

idx_lanesCC=find(binsizes>=3); % indices w.r.t bins of the connected components that have three vertices or more and thus are considered to be part of some lane (i.e. connected components having 2 edges or more).

%IMPORTANT: after performing the removal of connected components having only 2 vertices (length 1) it is possible to
%recover an empty queuing structure. In that case the queuing metric is set to be infinity
if isempty(idx_lanesCC)
    % In case the final cuttedForest is an empty graph we set
    % queueMetric=Inf, queueLengthDistri=0; and queueStructure=[];
    queueMetric=Inf;
    queueLengthDistri=0;
    queueStructure=[];
else
    %In case the cuttedForest is not empty we calculate our metric over
    %each of its connected components (lanes)
    weightsCF=cuttedForest.Edges.Weight;
    endNodesCF=cuttedForest.Edges.EndNodes;
    laneLengths=binsizes(idx_lanesCC).'-1; %array with the topological length of each lane (number of edges in each lane)
    
    nbAgents=sum(laneLengths+1); % Number of agents within the subpopulation that are in the queing structure.
    nbLanes=length(idx_lanesCC); % number of distinct lanes (connected components of the cutted forest).
    
    idx_agentsLanes=bins(idx_lanesCC);% indices w.r.t [x,y] of agents in each lane
    laneSumWeight=zeros(nbLanes,1); %preallocate for the sum of weights of each lane.
    for cc=1:length(idx_lanesCC)
        idx_edgesLanes=any(ismember(endNodesCF,idx_agentsLanes{cc}),2); % index w.r.t cuttedForest.Edges of the edges that belong to the connected component (lane) cc.
        laneSumWeight(cc)=sum(weightsCF(idx_edgesLanes)); % sum of the weights of all edges in the connected lane cc
    end
    queueMetric=((1./laneLengths.^2).'*laneSumWeight)/(nbAgents*nbLanes); % coefficients (ROW vector)
    queueLengthDistri=laneLengths;
    queueStructure=cuttedForest;
end




end

%%%%%%%%%%%%%%%%%%%%%%% Details on our queueuing structure construction %%%%%%%%%%%%%%%%%%%%%
%      First, we aim to compute the minimal weight spanning forest (MWSF) of each populations (moving RIGHT and LEFT respectively),
%      i.e. a subgraph of the Delaunay traingulation restricted to each population
%      respectively. In other words, the minimal weight spanning tree (MWST) of
%      each connected component.
%
%      Specifically, the weight of each EDGE in the original connectivity graph is calculated as follow:
%      a) we MOMENTARILY provide the edge with a direction pointing towards
%      the target (i.e. the dot product of the directed version of the edge with
%      the homing vector should always be non-negative).
%
%      b) For both endpoints of each edge we compute the angles made
%      between the CURRENT velocities and the directed version of the edge
%      from a).
%
%      c) We endow each edge with a weight resulting from the mean of the
%      two angles computed in b).
%
%      IMPORTANT: The resulting weighted graph that at is constructed at this point is UNDIRECTED.
%      The preferred direction for queues to be pointing towards the
%      target is already encoded in the directed version of the edge we used in step a).
%
%      Next, we run the minimum weight spanning tree algorithm (e.g.
%      the famous Prim's or Kruskal's algorithms) to obtain a minimum
%      weight spanning tree for each connected component of the weighted graphs
%      defined up to step c)...the resulting disconnected graph is called a minimum weight spanning
%      FOREST (MWSF).
%
%      After the MWSF are obtained
%      d) We sort all edges in the MWSF by decreasing order of weights and ITERATIVELY
%      remove edges from the sorted list until some stopping criteria is met. In this
%      particular implementation, the stopping criteria is when there are no nodes
%      of degree >=3 left in the cutted forest.
%
%      e) In general, we recover a large number of distinct connected components in the cutted forest
%      compared to the starting MWSF. In particular we follow the convention that a connected component
%      with only two agents (only one edge) is not representative of a lining-up of agents and thus these
%      are removed from the graph.
%      After this final removal, we are left with a (cutted) forest whose connected components represent each a
%      different line-up of 3 agents or more. This is our final queuing structure and it is a subgraph of
%      both the Delaunay triangulation of the overall population as well as the restricted Delaunay triangulation
%      to the specific population considered.
%
%      MOTIVATION: The MWSF will tend to choose edges that are as
%      "horizontal" as possible and whose endpoints (agents) have velocities poiting towards the target.
%      Conversely, the largest possible weight an edge can have is pi;
%      this happens whenever both agents move in the exact opposite
%      direction of the directed version of the edge that points in the
%      "direction" of the target.
%      Next, the removal of decreasingly sorted edges by their weights
%      produces a queing structure containing only edges that are "sufficiently"
%      horizontal.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

