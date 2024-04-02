%% Function to save during parallelization
function save_ref(fout, vol)
% Save the stacked ref matrix
% INPUTS:
%   fout (string): the filepath for the output file to save
%   vol (double matrix): the stack of images to save

% Save the output
save(fout,'vol','-v7.3')

end