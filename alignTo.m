function a = alignTo(U,nbhd,transition_func)

arguments
    U (:,2) double
    nbhd (:,1) cell {mustBeEqualSize(nbhd,U)}
    transition_func (1,:) char
end

N = size(U,1);
a = zeros(N,2);
g = @(x) transition(x,transition_func);
w = @(x) g(1/pi*acos(x));
for i=1:N
    ui = U(i,:);
    ui = ui/norm(ui,2);
    ui(isnan(ui)) = 0;
    J = nbhd{i};
    for j=1:size(J,2)
        uj = U(J(j),:);
        uj = uj/norm(uj,2);
        uj(isnan(uj)) = 0;
        s = dot(ui,uj);
        % due to numerical error, it can happen that |s|>1, if this is not
        % corrected, w(s) (which calls acos) will attempt to input an
        % imaginary value to g. The complex-valued acos function is
        % continuous at ±1 so this correction is well-behaved.
        if abs(s) > 1
            s = sign(s);
        end
        a(i,:) = a(i,:) + w(s)*uj;
    end
end
a = 1/6*a;

end

% validation
function mustBeEqualSize(a,b)
    if size(a,1) ~= size(b,1)
        eid = 'Size:equalSize';
        msg = 'second input must have same length as first';
        throwAsCaller(MException(eid,msg));
    end
end