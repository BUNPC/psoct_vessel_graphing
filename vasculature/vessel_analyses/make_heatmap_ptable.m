function [ptable] = make_heatmap_ptable(varargin)
%MAKE_PTABLE organize p-values into table
%   Organize the pairwise hypothesis test p-values into a table. The first
%   input is the pairs of comparisons, e.g. 'AD vs. CTE'. The remaining
%   arguments are the metrics, e.g. 'LengthDensity', 'BranchDensity', etc.
%   
%   In addition, include the summary statistics (mean, median, variance)
%   for each metric.
%
%   INPUTS:
%       varargin.Pairs (cell array): names of pairwise comparisons
%   OUTPUTS:
%       ptable (table): p-values organized into a table structure

% Retrieve the names of each pairwise comparison
Pairs = varargin{1};

% Rerieve the metrics that occur for all comparisons
VolumeFraction = varargin{2};
LengthDensity = varargin{3};
BranchDensity = varargin{4};
Tortuosity = varargin{5};
Diameter = varargin{6};

% Create summary statistics table
ptable = table(Pairs, VolumeFraction, LengthDensity, BranchDensity,...
    Tortuosity, Diameter);

end