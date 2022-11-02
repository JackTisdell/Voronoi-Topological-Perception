function P = voronoiPressure(DT,p)

arguments
    DT delaunayTriangulation
    p (1,1) double = 1
end

P = 0;
X = DT.Points;
[V,C] = voronoiDiagram(DT);
N = size(X,1);

for i = 1:N
    cell = C{i};
    verts = V(cell,:);
    A = poly_area(verts(:,1),verts(:,2));
    if isnan(A)
        A = Inf;
    end
    P = P+1/A^p;
end

P = P^(1/p);