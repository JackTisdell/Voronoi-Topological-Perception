function rmm = ringDists(X,c) 

arguments
    X (:,2) double {mustBeFinite}
    c (1,2) double {mustBeFinite} = mean(X)
end

rmm = median(vecnorm(X-c,2,2));

end
