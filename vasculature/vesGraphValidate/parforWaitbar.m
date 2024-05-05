function parforWaitbar(waitbarHandle,iterations)
% Display a wait bar to track progress of loop removal
    persistent count h N
    
    if nargin == 2
        % Initialize
        
        count = 0;
        h = waitbarHandle;
        N = iterations;
    else
        % Update the waitbar
        
        % Check whether the handle is a reference to a deleted object
        if isvalid(h)
            count = count + 1;
            waitbar(count / N,h);
            % Print the iteration to the console
            fprintf('\nFinished %i out of %i subgraphs\n',count, N)
        end
    end
end