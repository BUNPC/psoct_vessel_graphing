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

% TODO:
%{
- determine correct way to have variable number of input arguments for the
metrics
- 
%}

% Retrieve the names of each pairwise comparison
pairs = varargin.Pairs;

% TODO Rerieve the metrics
ptable = table(pairs, varargin(2:end));

% TODO Generate table

end