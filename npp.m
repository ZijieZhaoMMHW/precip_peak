function varargout = npp(precip, date_used, varargin)
%IDENTIFY_NPP_FROM_PRECIP Identify precipitation peak number from daily data.
%
%   OUT = IDENTIFY_NPP_FROM_PRECIP(PRECIP, DATE_USED) applies the NPP
%   detection algorithm to a gridded daily precipitation array PRECIP
%   with dimensions [nrow, ncol, time]. DATE_USED can be a datevec array,
%   datetime array, or datenum vector with one date per time step. Leap days
%   are removed before fitting.
%
%   [NUMP, ORDERM, RSM, RM31, RM61, RM91, HARMONICSS] =
%   IDENTIFY_NPP_FROM_PRECIP(...) returns the original output variables used
%   by the CPC subsection workflow.
%
%   Optional name-value arguments:
%       SmoothWindow          Gaussian smoothing window in days (default: 31)
%       Period                Harmonic period in days (default: 365)
%       PeakDropThreshold     Minimum peak-to-trough decrease (default: 0.2)
%       MergeDaysStrict       Always merge peaks closer than this (default: 30)
%       MergeDaysLoose        Merge close and similar peaks below this (default: 45)
%       MergePeakDiffThreshold Similarity threshold for close peaks (default: 0.1)
%       RainyPeriodThresholds Minimum rainy-period durations (default: [31 61 91])
%       MissingValueThreshold Absolute values above this are set to zero (default: 1e10)
%       MinNonzeroCount       Minimum nonzero daily samples per grid cell (default: 10)
%       Verbose               Print progress (default: false)
%


p = inputParser;
p.addParameter('SmoothWindow', 31, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Period', 365, @(x) isnumeric(x) && isscalar(x));
p.addParameter('PeakDropThreshold', 0.2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MergeDaysStrict', 30, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MergeDaysLoose', 45, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MergePeakDiffThreshold', 0.1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('RainyPeriodThresholds', [31 61 91], @(x) isnumeric(x) && isvector(x));
p.addParameter('MissingValueThreshold', 1e10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinNonzeroCount', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Verbose', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

if ndims(precip) ~= 3
    error('PRECIP must be a 3-D array with dimensions [nrow, ncol, time].');
end

if nargin < 2 || isempty(date_used)
    date_used = [];
else
    date_used = normalize_dates(date_used);
    if size(date_used, 1) ~= size(precip, 3)
        error('DATE_USED must have one row/date per time step in PRECIP.');
    end
    idx_rm = date_used(:, 2) == 2 & date_used(:, 3) == 29;
    precip(:, :, idx_rm) = [];
    date_used(idx_rm, :) = [];
end

precip(abs(precip) > opt.MissingValueThreshold) = 0;

nrow = size(precip, 1);
ncol = size(precip, 2);

orderm = NaN(nrow, ncol);
nump = NaN(nrow, ncol);
rsm = NaN(nrow, ncol);
harmonicss = NaN(nrow, ncol, opt.Period);

rainy_thresholds = opt.RainyPeriodThresholds(:)';
rainy_counts = NaN(nrow, ncol, numel(rainy_thresholds));

[x, y] = find(nanmean(precip, 3) ~= 0 & nansum(precip ~= 0, 3) >= opt.MinNonzeroCount);

if opt.Verbose
    fprintf('Processing %d valid grid cells...\n', length(x));
end

for i = 1:length(x)
    if opt.Verbose && (mod(i, 1000) == 0 || i == 1 || i == length(x))
        fprintf('Grid cell %d of %d\n', i, length(x));
    end

    ts_here = squeeze(smoothdata(squeeze(precip(x(i), y(i), :)), 'gaussian', opt.SmoothWindow));

    order = 1;
    X = design_harmonics_local(1:length(ts_here), opt.Period, order);
    [b, bint, ~, ~, stats] = regress(ts_here, X);
    n = length(ts_here);
    rss = nansum((X * b - ts_here).^2);
    aics = n * log(rss / n) + 2 * size(X, 2);
    rs = stats(1);
    orders = order;
    bint = bint(:, 1) .* bint(:, 2);

    while any(bint((end - 1):end) > 0) || order == 1
        order = order + 1;
        X = design_harmonics_local(1:length(ts_here), opt.Period, order);
        [b, bint, ~, ~, stats] = regress(ts_here, X);
        bint = bint(:, 1) .* bint(:, 2);
        n = length(ts_here);
        rss = nansum((X * b - ts_here).^2);
        aics = [aics; n * log(rss / n) + 2 * size(X, 2)];
        orders = [orders; order];
        rs = [rs; stats(1)];
    end

    if length(aics) > 1
        aics = aics(1:(end - 1));
        orders = orders(1:(end - 1));
        rs = rs(1:(end - 1));
    end

    [~, laic] = nanmin(aics);
    orders = orders(laic);
    rs = rs(laic);
    orderm(x(i), y(i)) = orders;
    rsm(x(i), y(i)) = rs;

    X = design_harmonics_local(1:length(ts_here), opt.Period, orders);
    coefs = X \ ts_here;
    harmonics = X * coefs(:);
    harmonics = harmonics(1:opt.Period);
    harmonicss(x(i), y(i), :) = harmonics;

    bw = bwconncomp(harmonics >= nanmean(harmonics));
    if ~isempty(bw.PixelIdxList) && bw.PixelIdxList{1}(1) == 1 && bw.PixelIdxList{end}(end) == opt.Period
        bw.PixelIdxList{1} = [bw.PixelIdxList{1}; bw.PixelIdxList{end}];
        bw.PixelIdxList{end} = [];
    end
    for rt = 1:numel(rainy_thresholds)
        rainy_counts(x(i), y(i), rt) = nansum(cellfun(@length, bw.PixelIdxList) >= rainy_thresholds(rt));
    end

    [pks, locs] = findpeaks(harmonics);
    idx_here = pks >= nanmean(harmonics);
    pks = pks(idx_here);
    locs = locs(idx_here);

    if harmonics(end) >= harmonics(end - 1) && harmonics(end) >= harmonics(1)
        pks = [pks; harmonics(end)];
        locs = [locs; opt.Period];
    elseif harmonics(1) >= harmonics(end) && harmonics(1) >= harmonics(2)
        pks = [harmonics(1); pks];
        locs = [1; locs];
    end
    idx_here = pks >= nanmean(harmonics);
    pks = pks(idx_here);
    locs = locs(idx_here);

    [~, locs_l] = findpeaks(-harmonics);
    if harmonics(end) <= harmonics(end - 1) && harmonics(end) <= harmonics(1)
        locs_l = [locs_l; opt.Period];
    elseif harmonics(1) <= harmonics(end) && harmonics(1) <= harmonics(2)
        locs_l = [1; locs_l];
    end

    troughs = NaN(length(locs), 1);
    for lc = 1:length(locs)
        dist_lh = locs_l - locs(lc);
        dist_lh_copy = dist_lh;
        dist_lh(abs(dist_lh) > 182 & dist_lh > 0) = dist_lh(abs(dist_lh) > 182 & dist_lh > 0) - opt.Period;
        dist_lh(abs(dist_lh) > 182 & dist_lh < 0) = opt.Period + dist_lh(abs(dist_lh) > 182 & dist_lh < 0);

        if length(dist_lh) > 1
            if any(dist_lh > 0) && any(dist_lh < 0)
                loc_b = locs_l(dist_lh == nanmax(dist_lh(dist_lh < 0)));
                loc_a = locs_l(dist_lh == nanmin(dist_lh(dist_lh > 0)));
            else
                [~, os] = sort(abs(dist_lh_copy));
                loc_b = locs_l(os(1));
                loc_a = locs_l(os(2));
            end

            troughb_n = (harmonics(locs(lc)) - harmonics(loc_b)) ./ (harmonics(locs(lc)));
            trougha_n = (harmonics(locs(lc)) - harmonics(loc_a)) ./ (harmonics(locs(lc)));
            troughs(lc) = nanmax([troughb_n trougha_n]);
        else
            troughs(lc) = (harmonics(locs(lc)) - harmonics(locs_l)) ./ (harmonics(locs(lc)));
        end
    end

    locs(troughs <= opt.PeakDropThreshold) = [];
    pks(troughs <= opt.PeakDropThreshold) = [];
    troughs(troughs <= opt.PeakDropThreshold) = [];

    idx_rm = [];
    for l = 1:length(pks)
        locs_here = (locs(l) - 15):(locs(l) + 15);
        locs_here = mod(locs_here - 1, opt.Period) + 1;

        if pks(l) < nanmax(harmonics(locs_here))
            idx_rm = [idx_rm; l];
        end
    end

    pks(idx_rm) = [];
    locs(idx_rm) = [];
    troughs(idx_rm) = [];

    pks = pks(:);
    locs = locs(:);
    troughs = troughs(:);

    if length(locs) > 1
        keep_merging = true;

        while keep_merging && length(locs) > 1
            [locs, r] = sort(locs);
            pks = pks(r);
            troughs = troughs(r);

            m = length(locs);
            dif_locs = diff([locs; locs(1) + opt.Period]);
            dif_pks = NaN(m, 1);

            for l = 1:m
                if l < m
                    seg = harmonics(locs(l):locs(l + 1));
                else
                    seg = [harmonics(locs(l):opt.Period); harmonics(1:locs(1))];
                end

                maxp = nanmax(seg);
                minp = nanmin(seg);
                dif_pks(l) = (maxp - minp) ./ maxp;
            end

            idx_pot = find((dif_locs < opt.MergeDaysLoose & dif_pks <= opt.MergePeakDiffThreshold) | ...
                (dif_locs < opt.MergeDaysStrict), 1, 'first');

            if isempty(idx_pot)
                keep_merging = false;
            else
                l1 = idx_pot;
                l2 = mod(idx_pot, m) + 1;

                pks_add = (pks(l1) + pks(l2)) / 2;
                troughs_add = (troughs(l1) + troughs(l2)) / 2;

                if l1 < m
                    locs_add = round((locs(l1) + locs(l2)) / 2);
                else
                    locs_add = mod(round((locs(l1) + locs(l2) + opt.Period) / 2) - 1, opt.Period) + 1;
                end

                rm_idx = sort([l1, l2], 'descend');
                pks(rm_idx) = [];
                locs(rm_idx) = [];
                troughs(rm_idx) = [];

                pks = [pks; pks_add];
                locs = [locs; locs_add];
                troughs = [troughs; troughs_add];
            end
        end
    end

    [locs, r] = sort(locs);
    pks = pks(r);
    troughs = troughs(r);

    nump(x(i), y(i)) = length(pks);
end
nump(nump==0)=1;

out = struct();
out.nump = nump;
out.orderm = orderm;
out.rsm = rsm;
out.harmonicss = harmonicss;
out.rainy_counts = rainy_counts;
out.rainy_thresholds = rainy_thresholds;
out.valid_rows = x;
out.valid_cols = y;
out.date_used = date_used;

rm31 = get_rainy_count(rainy_counts, rainy_thresholds, 31);
rm61 = get_rainy_count(rainy_counts, rainy_thresholds, 61);
rm91 = get_rainy_count(rainy_counts, rainy_thresholds, 91);
out.rm31 = rm31;
out.rm61 = rm61;
out.rm91 = rm91;

if nargout <= 1
    varargout = {out};
else
    varargout = {nump, orderm, rsm, rm31, rm61, rm91, harmonicss};
end

end

function X = design_harmonics_local(t, period, H)
t = double(t(:));
w = 2 * pi / period;
X = ones(length(t), 1 + 2 * H);

for k = 1:H
    X(:, 2 * k) = cos(k * w * t);
    X(:, 2 * k + 1) = sin(k * w * t);
end
end

function date_vec = normalize_dates(date_used)
if isdatetime(date_used)
    date_vec = datevec(date_used);
elseif isvector(date_used) && size(date_used, 2) ~= 6
    date_vec = datevec(date_used(:));
else
    date_vec = date_used;
end
end

function rm = get_rainy_count(rainy_counts, rainy_thresholds, threshold)
idx = find(rainy_thresholds == threshold, 1);
if isempty(idx)
    rm = NaN(size(rainy_counts, 1), size(rainy_counts, 2));
else
    rm = rainy_counts(:, :, idx);
end
end

