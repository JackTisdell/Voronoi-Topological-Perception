function pol = polarization(U)
arguments
    U (:,2) double
end

pol = norm(sum(U),2)/sum(vecnorm(U,2,2));

end