function compress_results(results_dir, varargin)

RESULTS_PREFIX = 'RES';
EPS = 1e-6;
[results, result_names ] = aggregate_results(results_dir, RESULTS_PREFIX,  varargin{:});


for ii=1:length(results)
    % Get  result files
    D = dir(sprintf('%s/%s/%s*_rseed*.it',results_dir,result_names{ii},RESULTS_PREFIX));
    numFiles = length(D(not([D.isdir])));
    
    if(numFiles == 0)
        fprintf('%s/%s no results files, skipping...\n',results_dir,result_names{ii})
    elseif(numFiles == 1)
        fprintf('%s/%s seems to be compressed already, skipping...\n',results_dir,result_names{ii})
        continue;
    end
    
    outfilename = sprintf('%s/%s/%s_%dJobs_rseed0000.it',results_dir,result_names{ii},result_names{ii}, numFiles);
    
    % Because MATLAB is stupid, result field member need to be hardcoded
    % the following way. We add EPS so that integers are saved as
    % doubles, which eliminates overflows for integers larger than 2^31
    sim_SNRdB              = results(ii).sim_SNRdB;
    sim_Nframes            = results(ii).sim_Nframes;
    sim_Ndatabits          = results(ii).sim_Ndatabits;
    sim_frame_errors       = results(ii).sim_frame_errors;
    sim_data_bit_errors    = results(ii).sim_data_bit_errors;
    sim_uncoded_bit_errors = results(ii).sim_uncoded_bit_errors;
    
    ldpc_nchk           = results(ii).ldpc_nchk;
    ldpc_nvar           = results(ii).ldpc_nvar;
    ldpc_code_rate      = results(ii).ldpc_code_rate;
    
    sim_SNRdB(sim_SNRdB > 2^31) = sim_SNRdB(sim_SNRdB > 2^31) + EPS;
    sim_Nframes(sim_Nframes > 2^31) = sim_Nframes(sim_Nframes > 2^31) + EPS;
    sim_Ndatabits(sim_Ndatabits > 2^31) = sim_Ndatabits(sim_Ndatabits > 2^31) + EPS;
    sim_frame_errors(sim_frame_errors > 2^31) = sim_frame_errors(sim_frame_errors > 2^31) + EPS;
    sim_data_bit_errors(sim_data_bit_errors > 2^31) = sim_data_bit_errors(sim_data_bit_errors > 2^31) + EPS;
    sim_uncoded_bit_errors(sim_uncoded_bit_errors > 2^31) = sim_uncoded_bit_errors(sim_uncoded_bit_errors > 2^31) + EPS;
    ldpc_nchk(ldpc_nchk > 2^31) = ldpc_nchk(ldpc_nchk > 2^31) + EPS;
    ldpc_nvar(ldpc_nvar > 2^31) = ldpc_nvar(ldpc_nvar > 2^31) + EPS;
    ldpc_code_rate(ldpc_code_rate > 2^31) = ldpc_code_rate(ldpc_code_rate > 2^31) + EPS;
 
    
    itsave(outfilename, sim_SNRdB, ...
                        sim_Nframes, ...
                        sim_Ndatabits, ...
                        sim_frame_errors, ...
                        sim_data_bit_errors, ...                 
                        sim_uncoded_bit_errors, ...
                        ldpc_nchk, ...
                        ldpc_nvar, ...
                        ldpc_code_rate );
                    
    if(isfield('gitversion', results(ii)))
        gitversion = results(ii).gitversion;
        itsave(outfilename, gitversion);
    end
    
    % Delete old results
    for jj=1:numFiles
        delete(sprintf('%s/%s/%s', results_dir, result_names{ii}, D(jj).name));
    end

end

