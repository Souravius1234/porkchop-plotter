function G = porkchop_plotter(depName, arrName, t1_start, t1_end, t2_start, t2_end, varargin)
% PORKCHOP_PLOTTER  Porkchop with C3 colormap + TOF, vinf+ (dep), vinf- (arr) contours
% Call your kernels once before using this function.
%
% Example:
% G = porkchop_plotter('EARTH','MARS', ...
%       '2026-08-01','2027-04-01', '2026-10-01','2028-01-01', ...
%       'Ndep',141,'Narr',141,'LW','both','TOFminDays',50,'TOFmaxDays',500);
% 
% Author: Sourav Ghosh, Intelligent Space Systems Laboratory, The University of Tokyo.
% 

% ---------- Options ----------
p = inputParser; p.FunctionName = mfilename;
addParameter(p,'Ndep',121,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'Narr',121,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'LW','both',@(x)any(strcmpi(x,{'short','long','both'})));
addParameter(p,'Nrev',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Branch','l',@(x)any(strcmpi(x,{'l','r'})));
addParameter(p,'TOFminDays',10,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'TOFmaxDays',inf,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'SaveFig','',@(x)ischar(x)||isstring(x));
parse(p,varargin{:});
opt = p.Results;

% ---------- Time grids ----------
t1v = linspace_local(t1_start, t1_end, opt.Ndep);  % datetime row (UTC)
t2v = linspace_local(t2_start, t2_end, opt.Narr);
% SPICE likes: 'YYYY MON DD HH:MM:SS UTC' (uppercase)
t1str = cellstr(upper(string(t1v, 'yyyy MMM dd HH:mm:ss') + " UTC"));
t2str = cellstr(upper(string(t2v, 'yyyy MMM dd HH:mm:ss') + " UTC"));
t1et  = cellfun(@cspice_str2et, t1str);
t2et  = cellfun(@cspice_str2et, t2str);

% ---------- Constants from loaded kernels ----------
muSun = cspice_bodvrd('SUN','GM',1);    muSun = muSun(1);   % km^3/s^2
muDep = cspice_bodvrd(upper(depName),'GM',1); muDep = muDep(1);
muArr = cspice_bodvrd(upper(arrName),'GM',1); muArr = muArr(1);
radDep = mean(cspice_bodvrd(upper(depName),'RADII',3));     % km
radArr = mean(cspice_bodvrd(upper(arrName),'RADII',3));     % km
% (rLEO/rCap only needed if you later contour DV; not used for C3 plot)

% ---------- Planet heliocentric states ----------
r1s = zeros(3,opt.Ndep); v1p = zeros(3,opt.Ndep);
for i = 1:opt.Ndep, [r1s(:,i), v1p(:,i)] = pcp_get_rv(depName, t1et(i)); end
r2s = zeros(3,opt.Narr); v2p = zeros(3,opt.Narr);
for j = 1:opt.Narr, [r2s(:,j), v2p(:,j)] = pcp_get_rv(arrName, t2et(j)); end

% ---------- Allocate ----------
C3     = nan(opt.Ndep, opt.Narr);   % colormap
VinfD  = nan(opt.Ndep, opt.Narr);   % v∞+ at departure
VinfA  = nan(opt.Ndep, opt.Narr);   % v∞- at arrival
TOF_d  = nan(opt.Ndep, opt.Narr);
fail   = false(opt.Ndep, opt.Narr);

% ---------- Grid search ----------
for i = 1:opt.Ndep
    r1 = r1s(:,i); vp1 = v1p(:,i); et1 = t1et(i);
    for j = 1:opt.Narr
        et2 = t2et(j); dt = et2 - et1;
        if dt <= 0, fail(i,j) = true; continue; end
        TOFdays = dt/86400;
        if TOFdays < opt.TOFminDays || TOFdays > opt.TOFmaxDays
            fail(i,j) = true; continue;
        end
        r2 = r2s(:,j); vp2 = v2p(:,j);

        % Solve Lambert: try short/long as requested
        switch lower(opt.LW)
            case 'short', lwSet = 0;
            case 'long',  lwSet = 1;
            otherwise,    lwSet = [0 1];
        end

        best = []; bestC3 = inf;
        for lw = lwSet
            try
                [v1,v2,~,~,~,iter] = lambertI(r1, r2, dt, muSun, lw, 0, 'l');
                if isnan(v1(1)) || iter < 0, continue; end
                vinf_dep = v1 - vp1;  vinf_arr = v2 - vp2;
                C3cand   = dot(vinf_dep,vinf_dep);
                if C3cand < bestC3
                    bestC3 = C3cand;
                    best.vinfDep = norm(vinf_dep);
                    best.vinfArr = norm(vinf_arr);
                end
            catch
                % skip failures
            end
        end

        if isempty(best)
            fail(i,j) = true; continue;
        end

        C3(i,j)    = bestC3;
        VinfD(i,j) = best.vinfDep;
        VinfA(i,j) = best.vinfArr;
        TOF_d(i,j) = TOFdays;
    end
end

% ---------- Plot (C3 colormap) ----------
X = datenum(t1v); Y = datenum(t2v);
Z = C3; Z(fail) = NaN;

figure('Color','w'); hold on
[CX,CY] = meshgrid(X,Y);
contourf(CX,CY, Z.', 30, 'LineStyle','none');
cb = colorbar; colormap(parula);
ylabel(cb, 'C_3 (km^2/s^2)');

% TOF contours (black)
handles = []; labels = {};
if any(~isnan(TOF_d(:)))
    levTOF = nice_levels(TOF_d(~isnan(TOF_d)), 6);
    [Ctof,hTof] = contour(CX,CY, TOF_d.', levTOF, 'k', 'LineWidth', 0.9);
    clabel(Ctof,hTof,'Color','k','FontSize',8,'LabelSpacing',400);
    handles(end+1) = hTof; labels{end+1} = 'TOF (days)';
end

% v∞+ (departure) contours (red, solid)
if any(~isnan(VinfD(:)))
    levVD = nice_levels_v(VinfD(~isnan(VinfD)), 5);
    [Cvd,hVd] = contour(CX,CY, VinfD.', levVD, 'r-', 'LineWidth', 1.0);
    clabel(Cvd,hVd,'Color','r','FontSize',8,'LabelSpacing',400);
    handles(end+1) = hVd; labels{end+1} = 'v_\infty^+ (km/s)';
end

% v∞- (arrival) contours (green, dashed)
if any(~isnan(VinfA(:)))
    levVA = nice_levels_v(VinfA(~isnan(VinfA)), 5);
    [Cva,hVa] = contour(CX,CY, VinfA.', levVA, 'g--', 'LineWidth', 1.0);
    clabel(Cva,hVa,'Color','g','FontSize',8,'LabelSpacing',400);
    handles(end+1) = hVa; labels{end+1} = 'v_\infty^- (km/s)';
end

% Axes cosmetics
datetick('x','yyyy-mm','keeplimits'); datetick('y','yyyy-mm','keeplimits');
xlabel('Departure date (UTC)'); ylabel('Arrival date (UTC)');
title(sprintf('%s \\rightarrow %s porkchop (C_3 w/ TOF, v_\\infty^\\pm)', upper(depName), upper(arrName)));
grid on; box on; set(gca,'Layer','top');

% Global minimum marker (on C3)
[minVal,idx] = min(Z(:));
if ~isnan(minVal)
    [ii,jj] = ind2sub(size(Z),idx);
    plot(X(ii),Y(jj),'kp','MarkerSize',12,'MarkerFaceColor','y','LineWidth',1.2);
end

% Legend (only for what exists)
if ~isempty(handles)
    legend(handles, labels, 'Location','northeastoutside');
end

if ~isempty(opt.SaveFig), saveas(gcf, char(opt.SaveFig)); end

% ---------- Output ----------
G.dep_dates = t1v; G.arr_dates = t2v;
G.C3 = C3; G.vinf_dep = VinfD; G.vinf_arr = VinfA; G.TOF_days = TOF_d; G.fail = fail; G.options = opt;

end

% ===== helpers =====

function tvec = linspace_local(tstart, tend, N)
    if ~isa(tstart,'datetime'), tstart = datetime(tstart,'TimeZone','UTC'); else, tstart.TimeZone='UTC'; end
    if ~isa(tend,'datetime'),   tend   = datetime(tend,  'TimeZone','UTC'); else, tend.TimeZone='UTC'; end
    tvec = tstart + (0:N-1).*(tend - tstart)/(N-1);
end

function lev = nice_levels(vals, n)
    vals = vals(:); vals = vals(isfinite(vals));
    if isempty(vals), lev = []; return; end
    vmin = prctile(vals,5); vmax = prctile(vals,95);
    if ~isfinite(vmin) || ~isfinite(vmax) || vmin==vmax
        lev = unique(round(vals)); return;
    end
    lev = linspace(vmin, vmax, n);
    lev = round(lev/10)*10;      % days, nearest 10
    lev = unique(lev);
end

function lev = nice_levels_v(vals, n)
    vals = vals(:); vals = vals(isfinite(vals));
    if isempty(vals), lev = []; return; end
    vmin = prctile(vals,10); vmax = prctile(vals,90);
    if ~isfinite(vmin) || ~isfinite(vmax) || vmin==vmax
        lev = unique(round(vals,2)); return;
    end
    lev = linspace(vmin, vmax, n);
    lev = round(lev,2);          % km/s, ~0.01 precision (adjust if you like)
    lev = unique(lev);
end
