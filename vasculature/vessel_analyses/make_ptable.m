function [ptable] = make_ptable(varargin)
%MAKE_PTABLE organize p-values into table
%   Organize the pairwise hypothesis test p-values into a table. The first
%   input is the pairs of comparisons, e.g. 'AD vs. CTE'. The remaining
%   arguments are the metrics, e.g. 'LengthDensity', 'BranchDensity', etc.
%
%   INPUTS:
%       varargin.Pairs (cell array): names of pairwise comparisons
%   OUTPUTS:
%       ptable (table): p-values organized into a table structure

% Retrieve the names of each pairwise comparison
Pairs = varargin{1};

% Rerieve the metrics that occur for all comparisons
LengthDensity = varargin{2};
BranchDensity = varargin{3};
VolumeFraction = varargin{4};
TortOutliers = varargin{5};
Diameter = varargin{6};

% Inlcude tortuosity if present
if length(varargin) == 7
    Tortuosity = varargin{7};
    % Create p-value table with tortuosity
    ptable = table(Pairs, LengthDensity, BranchDensity, VolumeFraction,...
                    TortOutliers, Diameter,Tortuosity);
else
    % Create p-value table without tortuosity
    ptable = table(Pairs, LengthDensity, BranchDensity, VolumeFraction,...
                    TortOutliers, Diameter);
end

end