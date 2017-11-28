function [  ] = generate_csv( source_foldername, outname, varargin)
%this function extracts the BER FER and generates CSV files in the form
% X; Y; Y_confidence_min, Y_confidence_max

% Assumptions: 
% - source_foldername is the name to the folder (relative to the current working directory
%   and without a trailing slash!). In the folder  there should be a single, compressed .mat file
% - The file names are relative to the current working directory.
% - Optionally, the condidence interval can be configured with the third argument. default is 0.95


%% check if there is one .mat file in the folder
D=dir([source_foldername '/RES_*.it']);

if length(D)>1
    error('There is more than one mat file in the source folder. Did you compress your results?')
end

%% Load the file
filename = [ source_foldername, '/' , D(1).name];
itload(filename);

%% Compute confidence intervalls:
% set confidence
alpha = 0.95;
if length(varargin)==1
    alpha = varargin{1};
elseif length(varargin)>1
   error('At most one optional argument allowed!')
end



nonzero_bitrate_idx = (sim_Ndatabits>0) & (sim_data_bit_errors>0);
nonzero_framerate_idx = (sim_Nframes>0) & (sim_frame_errors>0);

frame_errors = round(sim_frame_errors(nonzero_framerate_idx));
frames_num = round(sim_Nframes(nonzero_framerate_idx));
frames_snr = sim_SNRdB(nonzero_framerate_idx);

bit_errors = round(sim_data_bit_errors(nonzero_bitrate_idx));
bit_num = round(sim_Ndatabits(nonzero_bitrate_idx));
bit_snr = sim_SNRdB(nonzero_bitrate_idx);

[BER,BER_ci] = berconfint(bit_errors, bit_num, alpha);
[FER,FER_ci] = berconfint(frame_errors, frames_num, alpha);

% pgfplots likes to have the confidence intervals relative to the datapoints
BER_ci = [BER-BER_ci(:,1), BER_ci(:,2)-BER];
FER_ci = [FER-FER_ci(:,1), FER_ci(:,2)-FER];

%% write files
outnameBER = [outname '_BER.csv'];
outnameFER = [outname '_FER.csv'];



BER_data = [bit_snr,  BER, BER_ci];
FER_data = [frames_snr,  FER, FER_ci];

csvwrite(outnameBER, BER_data);
csvwrite(outnameFER, FER_data);



end

