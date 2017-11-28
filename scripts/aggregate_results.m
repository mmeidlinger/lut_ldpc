function [results, result_names ] = aggregate_results(results_dir, results_prefix,  varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Initialize output
results = struct([]);
result_names = cell.empty;

% If nothing is parsed compress everything in results folder
if isempty(varargin)
    D = dir( [results_dir '/' results_prefix '*']);
    noSim = length(D);
    for dd = 1:noSim
        result_names{dd} = D(dd).name;
    end
else
    result_names = varargin{1};
    % Find how many simulations we have
    noSim = length(result_names);
end



% Analyze each simulation separately
for ii = 1:noSim
    
	% Find how many files there are in the folder
    D = dir(sprintf('%s/%s/%s*.it', results_dir, result_names{ii},results_prefix));
    noFiles(ii) = length(D(not([D.isdir])));
    
    % Load each file and import results
    for jj = 1:noFiles(ii)
        
        itload(sprintf('%s/%s/%s', results_dir, result_names{ii}, D(jj).name));
        
        if( jj == 1 )
            % Initialize FER, BER and SNRdB
            results(ii).sim_SNRdB = sim_SNRdB;
            results(ii).sim_Nframes = sim_Nframes;
            results(ii).sim_Ndatabits = sim_Ndatabits;
            results(ii).sim_frame_errors = sim_frame_errors;
            results(ii).sim_data_bit_errors = sim_data_bit_errors;
            results(ii).sim_uncoded_bit_errors = sim_uncoded_bit_errors;
            
            results(ii).ldpc_nvar = ldpc_nvar;
            results(ii).ldpc_nchk = ldpc_nchk;
            results(ii).ldpc_code_rate = ldpc_code_rate;
            
            if( exist('runtime', 'var') )
                results(ii).runtime = runtime;
            end
            
            if( exist('gitversion', 'var') )
                results(ii).gitversion = gitversion;
            end
            
        else
            % Check consistency of SNR range and skip if inconsistent
            if( length(results(ii).sim_SNRdB) ~= length(sim_SNRdB) || sum(results(ii).sim_SNRdB == sim_SNRdB) ~= length(results(ii).sim_SNRdB) )
                warning('Inconsistent SNR ranges, skipping!');
                continue;
            end
            if(results(ii).ldpc_nchk ~= ldpc_nchk || ...
               results(ii).ldpc_nvar ~= ldpc_nvar || ...
               results(ii).ldpc_code_rate ~= ldpc_code_rate)
                    warning('Inconsistent LDPC code parameters!')
            end
            if( isfield('gitversion', results(ii)) && exist('gitversion', 'var') )
                if( strcmp(results(ii).gitversion, gitversion) ~= 0)
                    warning('Result files produced with different gitversions!');
                end
            end
            % Update FER and BER
            results(ii).sim_Nframes = results(ii).sim_Nframes + sim_Nframes;
            results(ii).sim_Ndatabits = results(ii).sim_Ndatabits + sim_Ndatabits;
            results(ii).sim_frame_errors = results(ii).sim_frame_errors + sim_frame_errors;
            results(ii).sim_data_bit_errors = results(ii).sim_data_bit_errors + sim_data_bit_errors;
            results(ii).sim_uncoded_bit_errors = results(ii).sim_uncoded_bit_errors + sim_uncoded_bit_errors;
            
            if( exist('runtime', 'var') )
                results(ii).runtime = results(ii).runtime + runtime;
            end
        end
        
    end        
end

ii=1;
while(ii<=noSim)
    if(noFiles(ii)==0)
        noFiles(ii)=[];
        results(ii)=[];
        result_names(ii)=[];
        noSim = noSim-1;
    else
        ii=ii+1;
    end
    
    
end

end

