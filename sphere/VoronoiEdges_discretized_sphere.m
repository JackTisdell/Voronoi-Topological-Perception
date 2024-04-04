function [Vx_discrete,Vy_discrete,Vz_discrete]=VoronoiEdges_discretized_sphere(V,C,resolution)
%%
% resolution=2*pi/180;
N=length(C);
Vx_discrete=cell(1,N); Vy_discrete=Vx_discrete; Vz_discrete=Vx_discrete;

for i=1:N
    nNeigh=length(C{i});
    VV=[V(C{i},1)'; V(C{i},2)'; V(C{i},3)']; %extracting vertice's coordinates of Voro cell i (3xnNeigh array)
    VV=[VV VV(:,1)]; %making the list of vertices cyclic (first and last vertices match each other)
    
    for j=1:nNeigh
        
        G=Arc(VV(:,j),VV(:,j+1),resolution);
        Vx_discrete{i}=[Vx_discrete{i}; G(1,1:end-1)'];
        Vy_discrete{i}=[Vy_discrete{i}; G(2,1:end-1)'];
        Vz_discrete{i}=[Vz_discrete{i}; G(3,1:end-1)'];
        
    end
    
end

end

function G = Arc(A,B,resolution)
% return an discretized arc between two points A and B

AxB = cross(A,B);
AdB = dot(A,B);
Ap = cross(AxB, A);
Ap = Ap/norm(Ap);
theta = atan2(sqrt(sum(AxB.^2,1)), AdB); % > 0
npnts = max(ceil(theta/resolution),2); % at least 2 points
theta = linspace(0, theta, npnts);
G = A*cos(theta) + Ap*sin(theta);
end % Arc