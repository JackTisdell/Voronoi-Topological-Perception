function y = nearestOnSegment(S,x)

arguments
    S (2,2) double
    x (2,1) double
end

p = S(:,1);
q = S(:,2);
a = -atan2(p(2)-q(2),p(1)-q(1));

R = [ cos(a) -sin(a); sin(a) cos(a) ];

X = R*[p q x];

while true
    [~,j] = min(X(1,1:2));
    if X(1,3) <= X(1,j)
        y0 = X(:,j);
        break;
    end

    [~,j] = max(X(1,1:2));
    if X(1,3) >= X(1,j)
        y0 = X(:,j);
        break;
    end
    
    y0 = [ X(1,3); X(2,1) ];
    break;
end

y = R'*y0;

end


