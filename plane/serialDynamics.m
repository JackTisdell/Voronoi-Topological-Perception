function serialDynamics(filename,N,L,nu,tmax,position_seed,angle_seed)

arguments
    filename (1,:) char
    N (1,1) double {mustBeInteger,mustBePositive}
    L (1,1) double {mustBeNonnegative}
    nu (1,1) double {mustBePositive}
    tmax (1,1) double {mustBeInteger,mustBePositive}
    position_seed (1,1) double {mustBeInteger} = 0
    angle_seed (1,1) double {mustBeInteger} = 0
end

%% Parameters
erase='';

fdim = 1;
M0 = L;
if fdim==2
    M0 = pi*L^2/2;
end

%% Initial conditions
ic_rad=.5*sqrt(N*pi/4/.91);
rng(position_seed);
X = ic_rad*(2*rand(N,2) - 1);
rng(angle_seed);
ang = 2*pi*rand(N,1);
U = [cos(ang) sin(ang)];

% ang = 2*pi*rand(N,1);
% n1 = rand(N,2);
% n2 = pi*rand(N,2);
% X = ic_rad*( [cos(ang) sin(ang)] + .2*n1 );
% U = [-sin(ang) cos(ang)] + .25*n2;

%% Initialize targets
tar = Target([0 0]);

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


save(filename,'N','L','nu','position_seed','angle_seed','DT_t','U_t');

end
