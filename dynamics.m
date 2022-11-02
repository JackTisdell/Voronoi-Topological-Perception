%% Parameters
N = 400;
L = 1;
nu = 2.5;
tmax = 1000;

display=true;
fixframe = false;
frameinrad = 50;

fdim = 1;
M0 = L;
if fdim==2
    M0 = pi*L^2/2;
end

filename='data.mat';
erase='';

%% Initial conditions
ic_rad=.5*sqrt(N*pi/4/.91);
rng(2);
X = ic_rad*(2*rand(N,2) - 1);
rng(18);
ang = 2*pi*rand(N,1);
U = [cos(ang) sin(ang)];

% ang = 2*pi*rand(N,1);
% n1 = rand(N,2);
% n2 = pi*rand(N,2);
% X = ic_rad*( [cos(ang) sin(ang)] + .2*n1 );
% U = [-sin(ang) cos(ang)] + .25*n2;

% sep = 10;
% X = zeros(N,2);
% X(:,2) = 4*rand(N,1) - 2;
% X(1:N/2,1) = -sep;
% X(N/2+1:end,1) = sep;
% U = zeros(N,2);
% U(1:N/2,1) = 1;
% U(N/2+1:end,1) = -1;

% rad = 3;
% sep = 4;
% r = rad*rand(N,1);
% a = 2*pi*rand(N,1);
% X = sqrt(r).*[cos(a) sin(a)];
% X(1:N/2,1) = X(1:N/2,1) + sep;
% X(N/2+1:end,1) = X(N/2+1:end,1) - sep;
% U = [-sin(a) cos(a)];

% X = zeros(N,2);
% X(1:N/2,2) = 5;
% X(N/2+1:end,2) = -5;
% X(:,1) = 100*rand(N,1) - 25;
% U = zeros(N,2);
% U(1:N/2,1) = 1;
% U(N/2+1,end) = -1;

%% Initialize targets

% sR = [sep sep; -4 4]';
% sL = [-sep -sep; -4 4]';
% R = [.8 .8 .6 .5; .7 .9 .9 .6]';
% tar = Target({sL,sR},[1 2]);
ids = ones(N,1);
ids(1:N/2) = 2;
% tar = Target(1*[1 -1 -1 1; 1 1 -1 -1]');
% tar = Target([0 0]);
% tar = Target({[-w/2 0],[w/2 0]});
tar_rad = 4;
tar_ang = 2*pi/3;
tar = Target({tar_rad*[1 0],tar_rad*[cos(tar_ang) sin(tar_ang)],tar_rad*[cos(2*tar_ang) sin(2*tar_ang)]});
% tar = Target();

%% Preallocation of some intermediate variables
[r,a,h0,h,s,l,rho] = deal(zeros(N,1));
U1 = zeros(N,2);
Q = zeros(N,2);

%% Preallocation of time series variables
U_t = zeros(N,2,tmax);
DT_t = cell(tmax,1);


for t = 1:tmax
%% Compute step

DT = delaunayTriangulation(X);
DT_t{t} = DT;
U_t(:,:,t) = U;
[nbhd, nearest, d] = neighborhoods(DT);
fun = @(x) transition(x,'expReciprocal');
s = arrayfun(fun, d/L);

% repulsion
r = X - X(nearest,:);
rnorm = vecnorm(r,2,2);
r = s .* r./rnorm;  % rescale to length s

% alignment
a = alignTo(U,nbhd,'expReciprocal');

% homing
h0 = homeToTarget(tar,X);
h = (1-s) .* h0./vecnorm(h0,2,2);   % rescale to length 1-s
h(isnan(h)) = 0;                    % protect 0 entries

% direction
U1 = (r + h + nu*a)./(1 + nu);

% speed
[~,l] = voronoiProjectToBoundary(DT,U1);
M = l;
if fdim==2
    M = voronoiForwardArea(DT,U1);
end

%% plot
if display==true
    clipToVoronoi=false;
    
    fig = figure(1);
    scatter(X(:,1), X(:,2), '.k');
    hold on
    % triplot(DT,'b');
    % text(X(:,1),X(:,2),num2cell(1:N));
    
    % vor = voronoi(DT);
    % vor(2).Color = 'k';
    
    qs = 1;
    % quiver(X(:,1),X(:,2),qs*r(:,1),qs*r(:,2),'off', 'r');
    % quiver(X(:,1),X(:,2),qs*h(:,1),qs*h(:,2),'off', 'g');
    % quiver(X(:,1),X(:,2),qs*a(:,1),qs*a(:,2), 'off', 'b');
    
    % quiver(X(:,1),X(:,2),qs*U(:,1),qs*U(:,2),'off', 'b');
    quiver(X(:,1),X(:,2),qs*U1(:,1),qs*U1(:,2), 'b');
    
    % plot([X(:,1) Q(:,1)]', [X(:,2) Q(:,2)]', 'r:', 'LineWidth', 1);
    % s = min(l);
    % quiver(X(:,1), X(:,2), s*U(:,1), s*U(:,2), 'off', 'r', 'ShowArrowHead', 1, 'LineWidth', 1);
    % scatter(Q(:,1), Q(:,2), '*b');
    
    C = tar.Components;
    for m=1:size(C,1)
        if isa(C{m},'polyshape')
            plot(C{m},'EdgeColor','none');
        end
    end
    
    
    hold off
    
    clip = [min(X,[],1); max(X,[],1)];
    if clipToVoronoi==true
        V = voronoiDiagram(DT);
        clip = [min([V; X],[],1); max([V(2:end,:); X],[],1)];
    end
    clip_width = 1.2*max(clip(2,:)-clip(1,:));
    clip_center = mean(clip,1);
    clip_lim(1:4) = clip_center(ceil((1:4)/2)) - rescale(mod(1:4,2),-1,1)*.5*clip_width;
    
    axis equal
    
    
    frame = [-frameinrad frameinrad -frameinrad frameinrad];
    if ~fixframe
        com = mean(X);
        rmed = sqrt(median((X(:,1)-com(1)).^2+(X(:,2)-com(2)).^2));
        frame = [com(1)-3*rmed com(1)+3*rmed com(2)-3*rmed com(2)+3*rmed];
    end
    axis(frame);
    title(['$t=',num2str(t),'$'], 'Interpreter', 'Latex');
end

%% Update

U = tanh(M/M0).*U1;
X = X + U;

%% display progress

barwidth = 60;
frac = floor(t/tmax * barwidth);
fracpct = t/tmax*100;
sofar = repmat(sprintf('\x0023'),1,frac);
yet = repmat('-',1,barwidth-frac);
progbar = ['[' sofar yet ']' ' ' sprintf('%04.1f\t%s',fracpct,filename) newline];
fprintf([erase progbar]);
erase = repmat(sprintf('\b'),1,length(progbar));
end


save("data.mat",'DT_t','U_t');
