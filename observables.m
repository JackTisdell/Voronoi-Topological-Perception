numruns=9;
tmax = 500;

% nu = zeros(numruns,1);
% nustr = ["04","08","12","16"];
N = zeros(numruns,1);
Nstr = ["400","600","800"];
seeds = ["_s1","_s2","_s3"];
ids = Nstr'+("_04"+seeds);

pols = zeros(tmax,numruns);
angmoms = zeros(tmax,numruns);
absangmoms = zeros(tmax,numruns);
Ps = zeros(500,numruns);
As = zeros(500,numruns);
Ms = zeros(500,numruns);
ms = zeros(500,numruns);
rmed = zeros(500,numruns);
tri_root_area = zeros(500,numruns);

for i=1:numel(ids)
    filename = strcat('data',ids(i),'.mat');
    load(filename);

%     pol_t = zeros(tmax,1);
%     angmom_t = zeros(tmax,1);
%     absangmom_t = zeros(tmax,1);
%     A_t = zeros(tmax,1);
%     
    for t=1:tmax
        DT = DT_t{t};
        X = DT.Points;
        U = U_t(:,:,t);
        
%         pol_t(t) = polarization(U);
%         [angmom_t(t),absangmom_t(t)] = angularMomentum(X,U);
%         A_t(t) = inwardTotalArea(DT);
%         com = mean(X);
%         rmed(t,i) = median(vecnorm(X-com,2,2));
%         tri = pointLocation(DT,com);
%         verts = DT.ConnectivityList(tri,:);
%         coords = X(verts,:);
%         tri_root_area(t,i) = sqrt(poly_area(coords(:,1),coords(:,2)));
        Ps(t,i) = voronoiPressure(DT,1);
    end

    rmm = mean(rmed);
    ra = mean(tri_root_area);
% 
%     pols(:,i) = pol_t;
%     angmoms(:,i) = angmom_t;
%     absangmoms(:,i) = absangmom_t;
%     As(:,i) = A_t;



end

% figure(1)
% plot(Ms)
% legend(strcat("nu=",["4" "8" "12" "16"]));
% 
% figure(2)
% plot([4 8 12 16]', [mean(Ms)' mean(ms)'])



