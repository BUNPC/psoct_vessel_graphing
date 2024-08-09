function [ad, cte, nc] = organize_metrics(metrics, subid, region, param)
%ORGANIZE_METRICS Separate metrics struct into local variables
%   The metrics struct has several substructures. The first layer is the
%   subject ID, the second layer is the region (tissue, wm, gm, etc.), and
%   the third layer is the actual metric. This script will retrieve the
%   desired metric across each subject and organize it into an array based
%   on the group (AD, CTE, NC).
% INPUTS:
%   metrics (struct): stores all metrics for all subject & regions
%   subid (cell array): list of subject IDs
%   region (string): tissue region to select from
%   param (string): parameter name to retrieve
%   ad_idcs (array): indices of the AD substructs
%   cte_idcs (array): indices of the CTE substructs 
%   nc_idcs (array): indices of the NC substructs 
% OUTPUTS:
%   ad (array): array of the metric from all AD subjects for this region
%   cte (array): array of the metric from all CTE subjects for this region
%   nc (array): array of the metric from all NC subjects for this region

%% Initialize indices
% Find indices of each group in the structure
f = fields(metrics);
ad_idcs = find(contains(f,'AD'));
cte_idcs = find(contains(f,'CTE'));
% Initialize arrays
ad = [];
cte = [];
nc = [];

%% Concatenate the AD, CTE, and NC

% The diameter is a column vector and must be transposed
if strcmp(param,'diameter') || strcmp(param,'tortuosity')
    for ii=1:length(subid)
        if ismember(ii, ad_idcs)
            ad = [ad; mean(metrics.(subid{ii}).(region).(param))];
        elseif ismember(ii, cte_idcs)
            cte = [cte; mean(metrics.(subid{ii}).(region).(param))];
        else
            nc = [nc; mean(metrics.(subid{ii}).(region).(param))];
        end
    end
else
    % Iterate through each subject
    for ii=1:length(subid)
        if ismember(ii, ad_idcs)
            ad = [ad; metrics.(subid{ii}).(region).(param)];
        elseif ismember(ii, cte_idcs)
            cte = [cte; metrics.(subid{ii}).(region).(param)];
        else
            nc = [nc; metrics.(subid{ii}).(region).(param)];
        end
    end
end

end