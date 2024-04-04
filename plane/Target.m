%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TARGET
% Defines the target class for implementation of VTP
%
% tar = Target()
%           defines an empty Target, that is, tar.Components and tar.ids
%           are the empty cell array and empty array, respectively.
% tar = Target(P)
%           where P is a 2-column matrix defines a Target with a single
%           component whose veritices are the rows of P. P is taken to
%           represent a point, segment, or (closed) polygon, depending on
%           the number of rows.
% tar = Target(C) 
%           where C is a cell array of size k where each C{:} is a 2-column
%           matrix defines a target with k components specified by the C{:}
%           as above.
% tar = Target(_,ids)
%           where ids is a vector sets tar.Components exactly as above but
%           also sets tar.Ids to ids. the length of ids must match the
%           number of components. The purpose is to associate distinct
%           components into unions effectively creating subtargets which
%           may be sought by distinct subpopulations. (See the homeToTarget
%           method.)
% Notes:
%   - If P = C{i} in the above syntaxes, then the i-th cell of
%     tar.Components contains:
%       - P itself if P has only one or two rows;
%       - polyshape(P) if P has at least three rows.
%   - homeToTarget may give unexpected results if the target has nonconvex
%     components.
%   - homeToTarget may give unexpected results if the target has
%     overlapping components with the same id.

classdef Target
    
properties
    Components
    Ids
end

methods
    function tar = Target(C,ids)
        arguments
            C {mustBeA(C,["cell","double"])} = cell(0)
            ids double {mustHaveSameNum(ids,C)} = ones(numel(C),1)
        end

        if isempty(C)
            tar.Components = cell(0);
            tar.Ids = [];
            return;
        end
        if isa(C,'double')
            C = {C};
        end

        k = numel(C);
        K = cell(k,1);

        for m=1:k
            P = C{m};
            if size(P,1) < 3
                K{m} = P;
                continue; 
            end

            pgon = polyshape(P);
            conv = convhull(pgon);
            if size(pgon.Vertices,1) ~= size(conv.Vertices,1)
                warning('target:nonconvex',...
                    'Possible nonconvex target region (component %d). May have unexpected results.', m);
            end
            K{m} = pgon;
        end

        tar.Components = K;
        tar.Ids = ids;
    end
    
    function tar = setComponent(tar,k,P)
        % Sets the k-th component of tar to (the point/segment/polygon)
        % specified by P. Throws a warning if the new component has
        % differnt combinatorial properties than the original.
        % A possible use case is moving targets.
        arguments
            tar (1,1) Target
            k (1,1) double {mustBeInteger,mustBePositive}
            P (:,2) double
        end
        K = tar.Components;
        if k > size(K,1)
            error('Second input must not exceed number of components in %s.', inputname(1));
        end
        Q = K{k};
        bothPointsOrSegs = isa(Q,'double') && isa(P,'double') && size(Q,1)==size(P,1);
        bothNGons = isa(Q,'polyshape') && isa(P,'polyshape')...
            && ( size(Q.Vertices,1)==size(P.Vertices,1)...
            || Q.numRegions == P.numRegions || Q.numHoles == P.numHoles );
        if ~(bothPointsOrSegs || bothNGons)
            warning('Changed combinatorial properties of target (component %d)', k);
        end
        if size(P,1) >= 2
            P = polyshape(P);
        end
        K{k} = P;
        tar.Components = K;
    end

    function h = homeToTarget(tar, X, target_idx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % h = homeToTarget(tar, X)
        %       where X is a 2-column matrix returns a matrix h of the same
        %       size as X whose i-th row h(i,:), the _homing vector_ is the 
        %       vector of minimum length such that X(i,:) + h(i,:) is in 
        %       the target tar. If X(i,:) is interior to a polygonal
        %       component of the target, then h(i,:) is the zero vector,
        %       naturally.
        % h = homeToTarget(tar, X, target_idx)
        %       where X is a 2-column matrix and target_idx is a vector
        %       whose length is the number of rows of X is like the
        %       previous syntax except that the i-th homing vector h(i,:) 
        %       is now the minimum legth vector such that X(i,:) + h(i,:)
        %       is contained in the union of components of target tar whose
        %       id (as specified by tar.Ids) matches target_idx(i). If
        %       target_idx(i) has no match among tar.Ids, then h(i,:) is
        %       the zero vector; this subpopulation effectively has no
        %       target.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        arguments
            tar (1,1) Target
            X (:,2) double
            target_idx (:,1) double {mustBeSameSize(target_idx,X)} = ones(size(X,1),1);
        end

        N = size(X,1);
        h = zeros(N,2);

        K = tar.Components;
        k = size(K,1);
        if k == 0, return; end

        ids = tar.Ids;

        d = Inf(N,1);

        for l=1:k
            P = K{l};
            J = (target_idx==ids(l));
            h0 = zeros(N,2);
            if isa(P, 'double') == 1
                if size(P,1) == 1   % P is a point
                    h0(J,:) = P - X(J,:);
                else                % P is a segment
                    for i=1:N
                        if J(i)==false, continue; end
                        y = nearestOnSegment(P',X(i,:)');
                        h0(i,:) = y'-X(i,:);
                    end
                end
            else                    % P is a polygon
                v = P.Vertices;
                for i=1:N
                    if J(i)==false || isinterior(P,X(i,:))
                        continue;
                    end
                    r = Inf;
                    for m=1:size(v,1)
                        seg = [v(m,:)' v(1+mod(m-2,size(v,1)),:)'];
                        y = nearestOnSegment(seg, X(i,:)');
                        e = y'-X(i,:);
                        r0 = dot(e,e);
                        if r0 < r
                            r = r0;
                            h0(i,:) = e;
                        end
                    end
                end
            end
            d0 = dot(h0,h0,2);
            I = J & (d0 < d);
            d(I) = d0(I);
            h(I,:) = h0(I,:);
        end    
    end

end

end
    
% validation
function mustHaveSameNum(C,I)
    if numel(I) ~= numel(C)
        eid = 'Size:mustBeEqual';
        msg = 'Inputs must have the same number of elements.';
        throwAsCaller(MException(eid,msg));
    end
end 
function mustBeSameSize(a,b)
    if size(a,1) ~= size(b,1)
        eid = 'Size:equalSize';
        msg = 'Third input must have as many entries as the number of ros of the first input';
        throwAsCaller(MException(eid,msg));
    end
end