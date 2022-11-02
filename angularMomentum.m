function [L,Labs] = angularMomentum(X,U,c)
arguments
    X (:,2) double
    U (:,2) double {mustBeSameSize(U,X)}
    c (1,2) double {mustBeFinite} = mean(X)
end

Y = X-c;   % recenter X about c
YxU = Y(:,1).*U(:,2) - Y(:,2).*U(:,1);
normProd = sum(vecnorm(Y,2,2).*vecnorm(U,2,2));
L = abs(sum(YxU))/normProd;
Labs = sum(abs(YxU))/normProd;

end

function mustBeSameSize(a,b)
if ~isequal(size(a),size(b))
    eid = 'Size:sameSize';
    msg = 'Fist two inputs must have the same size.';
    throwAsCaller(MException(eid,msg));
end
end