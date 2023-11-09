function [triggerUpdateFcn] = waitbarParfor(totalLoops, varargin)
%waitbarParfor Waitbar implementation for parfor loops.
%   [triggerUpdateFcn] = waitbarParfor(totalLoops) creates a waitbar for use in parfor loops. totalLoops indicates the number of loops until
%       completion. Use the triggerUpdateFcn to increment the waitbar count. Make sure to call this function in the outer parfor loop.
%   waitbarParfor(___, message) displays the specified message on the waitbar.
%   waitbarParfor(___, Name, Value) accepts all Name-Value pairs for the waitbar function. Use the "Name" option to set the window title.
%
%     Example
%         nLoops = 100;
% 
%         updateWaitbar = waitbarParfor(nLoops, "Calculation in progress...");
%         parfor loopCnt = 1:nLoops
%             A = rand(5000);
%             updateWaitbar(); %#ok<PFBNS>
%         end
%
% Author: Girmi Schouten (girmi.schouten@uantwerpen.be), 2019. 
% Written in MATLAB 2019b, tested on Ubuntu 18.04 & Windows 10.
    %% Parse input arguments
    argParser = inputParser();
    argParser.KeepUnmatched = true;
    addRequired(argParser, "totalLoops", @isscalar);
    defaultWaitbarMessage = "Please wait ...";
    addOptional(argParser, "waitbarMessage", defaultWaitbarMessage, @(str) isstring(str) || ischar(str));
    addOptional( argParser, "kappa", 0, @(x) isscalar(x));
    parse(argParser, totalLoops, varargin{:});
    totalLoops = argParser.Results.totalLoops;
    waitbarMessage = argParser.Results.waitbarMessage;
    kappa = argParser.Results.kappa;

    %% Initialize waitbar requirements
    
    parellelDataQueue = parallel.pool.DataQueue;
    afterEach(parellelDataQueue, @updateWaitbar);
    
%     waitbarHandle = waitbar(0, waitbarMessage);
    % Pass unmatched parameters to the waitbar
%     if ~isempty(fieldnames(argParser.Unmatched))
%         waitbarOptions = namedargs2cell(argParser.Unmatched);
%         set(waitbarHandle, waitbarOptions{:});
%     end
    
    triggerUpdateFcn = @updateProxy;
    loopCnt = 1;
    ratio = 1 / totalLoops;
    limit = ratio;
    lastsize = fprintf('Starting Psi MCMC updates...\n');
    
    %% Helper functions
    
    function updateWaitbar(~)
        
        progress = loopCnt / totalLoops;

        if( limit >= 0.01 )
            
            fprintf(repmat('\b', 1, lastsize));
            lastsize = fprintf( 'Psi MCMC updates %s%% done. %s\n', num2str( 100 * progress), kappa);
            limit = 0;
        end

        if( progress == 1 )
            fprintf(repmat('\b', 1, lastsize));
        end

        loopCnt = loopCnt + 1;
        limit = limit + ratio;
    end
    function updateProxy()
        send(parellelDataQueue, []);
    end
end

