function parsdata = MERGE_PARSDATA_ROWS(parsdata, rows_to_merge)
% "MERGE_PARS_ROWS" merges selected rows of parsdata with the row after
%   ...them and remove duplicate aggregates. It is a tool to combine...
%   ...snapshots of interest from a 2nd-stage aggregation library
% ----------------------------------------------------------------------- %
%   
%   Inputs/Output:
%   - parsdata: struct array of aggregate info with fields dpp, da,...
%       ...n_hyb, pp, etc.
%   - rows_to_merge: vector of indices (in the CURRENT parsdata) such that
%       ...each i in rows_to_merge is merged with i+1.
%
%   Example: rows_to_merge = [1 3] merges row 1 with 2 and row 3 with 4.
% ----------------------------------------------------------------------- %

% work from the bottom up so index shifting doesn't break things
rows_to_merge = sort(rows_to_merge(:), 'descend');

for kk = 1:numel(rows_to_merge)
    i = rows_to_merge(kk);
    j = i + 1;

    if j > numel(parsdata)
        error('merge_pars_rows:badIndex', ...
              'Row %d has no row after it to merge with.', i);
    end

    % concatenate aggregates from rows i and j
    merged = parsdata(i); % start from row i

    % assuming these fields exist in parsdata; modify if needed
    merged.dpp   = [parsdata(i).dpp; parsdata(j).dpp];
    merged.sigmapp   = [parsdata(i).sigmapp; parsdata(j).sigmapp];
    merged.da    = [parsdata(i).da; parsdata(j).da];
    merged.dm   = [parsdata(i).dm; parsdata(j).dm];
    merged.dg   = [parsdata(i).dg; parsdata(j).dg];
    merged.n_hyb = [parsdata(i).n_hyb; parsdata(j).n_hyb];
    merged.pp    = [parsdata(i).pp; parsdata(j).pp];
    merged.npp    = [parsdata(i).npp; parsdata(j).npp];

    nagg = numel(merged.pp);

    % find duplicates based on primary-particle IDs
    ind_flt = false(nagg,1);

    for a = 1:nagg-1
        if ind_flt(a), continue; end
        pa = sort(unique(merged.pp{a}(:,1)));

        for b = a+1:nagg
            if ind_flt(b), continue; end
            pb = sort(unique(merged.pp{b}(:,1)));

            if isequal(pa, pb)
                ind_flt(b) = true; % mark as duplicate
            end
        end
    end

    keep = ~ind_flt;

    % apply filter to each per-aggregate field
    merged.dpp = merged.dpp(keep);
    merged.sigmapp = merged.sigmapp(keep);
    merged.da = merged.da(keep);
    merged.dm = merged.dm(keep);
    merged.dg = merged.dg(keep);
    merged.n_hyb = merged.n_hyb(keep);
    merged.pp = merged.pp(keep);
    merged.npp = merged.npp(keep);

    % put merged row back into parsdata
    parsdata(i) = merged; % replace row i with merged version
    parsdata(j) = []; % delete row j

end

end
