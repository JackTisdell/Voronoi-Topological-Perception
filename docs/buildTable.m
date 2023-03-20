function S = buildTable(datafiles,tmax)

arguments
    datafiles (1,:) string
    tmax (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end

erase = '';

numruns = size(datafiles,2);

data = cell(1,numruns);
for r = 1:numruns
    msg = [sprintf('reading file %1.0f of %1.0f',r,numruns) newline];
    fprintf([erase msg]);
    erase = repmat(sprintf('\b'),1,length(msg));

    data{r} = load(datafiles(r),'-mat');
end


if tmax==0
    datum = data{1};
    tmax = size(datum.DT_t,1);
end

% initialize angular momentum, pressure, energy, median radius, root-area
[L, P, rmed, rra] = deal(nan(tmax,numruns));

for r = 1:numruns
    run = data{r};
    DT_t = run.DT_t;
    U_t = run.U_t;
    tmax_r = size(DT_t,1);
    if tmax_r > tmax
        warning('%s includes %d iterations. Call writeTable(_,%d) to include all samples',...
            datafiles(r),tmax_r,tmax_r);
    end
    for t = 1:tmax
        if t > tmax_r
            continue;
        end
        DT = DT_t{t};
        U = U_t(:,:,t);
        L(t,r) = angularMomentum(DT.Points,U);
        P(t,r) = voronoiPressure(DT);
        [rmedt,rrat] = ringDists(DT);
        rmed(t,r) = rmedt;
        rra(t,r) = rrat;

        %% display progress

        barwidth = 60;
        frac = floor(t/tmax * barwidth);
        fracpct = t/tmax*100;
        sofar = repmat(sprintf('\x0023'),1,frac);
        yet = repmat('-',1,barwidth-frac);
        progbar = [sprintf('processing file %1.0f of %1.0f\t',r,numruns) '[' sofar yet ']' ' ' sprintf('%04.1f',fracpct) newline];
        fprintf([erase progbar]);
        erase = repmat(sprintf('\b'),1,length(progbar));

    end
end

S.data = data;
S.angmom = L;
S.pressure = P;
S.medianRadius = rmed;
S.rootArea = rra;

end





