
function [SNRdB, FER, BER] = analyze_results(results_dir, varargin)

RESULTS_PREFIX = 'RES_';

% Close previous plots
close all

% Change default axes fonts.
fontSize = 12;
set(0,'DefaultAxesFontName', 'Times')
set(0,'DefaultAxesFontSize', fontSize)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times')
set(0,'DefaultTextFontSize', fontSize)

% Define plot colors
c = {'b', 'r', 'm', 'c', 'g', 'k', [0.9294 0.6941 0.1255], [0 0.4980 0], [0.2 0.6 1], [0.4941    0.1843    0.5569] };
% Define plot markers
m = {'s', 'o', 'v', '^',  '<', '>', 'h', 'p', '*', 'd'};

% Define linewidth
lineWidth = 2;

[results, result_names ] = aggregate_results(results_dir, RESULTS_PREFIX,  varargin{:});


for ii=1:length(results)
    FER{ii} = results(ii).sim_frame_errors./results(ii).sim_Nframes;
    BER{ii} = results(ii).sim_data_bit_errors./results(ii).sim_Ndatabits;
    SNRdB{ii} = results(ii).sim_SNRdB;
           
    % Plot results  
    figure(1)
    semilogy(SNRdB{ii}, FER{ii}, 'Color', c{mod(ii,length(c))+1}, 'Marker', m{mod(ii,length(m))+1}, 'LineWidth', lineWidth)
    legendStr{ii} = result_names{ii};
    hold on

    figure(2)
    semilogy(SNRdB{ii}, BER{ii}, 'Color', c{mod(ii,length(c))+1}, 'Marker', m{mod(ii,length(m))+1}, 'LineWidth', lineWidth)
    hold on     
end

% Beautify plots
figure(1)
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('Eb/N0 (dB)', 'FontSize', fontSize)
ylabel('FER', 'FontSize', fontSize)
grid on


figure(1);
legend(legendStr, 'Location', 'Best', 'Interpreter', 'none');
figure(2);
legendStr = {legendStr{1:end}, 'BIAWGN limit'};
legend(legendStr, 'Location', 'Best', 'Interpreter', 'none');

% Print Runtime results
for ii=1:length(results)
    if( isfield(results(ii), 'runtime') )
        fprintf(' Average runtime for simulation "%s" = %g s / frame\n', result_names{ii}, results(ii).runtime / sum(results(ii).sim_Nframes) )
    end
end


%% Error rate limit curve
% for n -> infty, the bit error probability is boubded as 
%(cf. Fig. 1.19 in Information theory, inference, and learning algorithms (Version 7.2 Cambridge University Press 2003), David J. Mackay,  )
% Pb >  H_2_inv(1 - C/R)
% For a Code of fixed rate R, we can thus compute the capacity C for
% various channel parameters such that R > C and then plot a BER limit
% curve.

Rate = results(1).ldpc_code_rate;

% Set capacity as function of channel parameter
c_func = @(sig) c_biawgn(sig);
%c_func = @(sig) c_awgn(sig);

H_2     = @(p) -p.*log2(p) - (1-p).*log2(1-p);
H_2_inv = @(p) fzero( @(x)(H_2(x)-p), [1e-16, .5] );


snr_min = -.01;

sig_max = fzero( @(x)(c_func(x)-Rate), 10^(-snr_min/20)/sqrt(2*Rate) );

snr_max = -20*log10(sig_max* sqrt(2*Rate));


snr_range = linspace(snr_min, snr_max, 1e2);

Pb_bound = zeros(1,length(snr_range));
for(ss=1:length(snr_range)-1)
    sig = 10^(-snr_range(ss)/20)/sqrt(2*Rate);
    Pb_bound(ss) = H_2_inv(1- c_func(sig)/Rate);
end
Pb_bound(end) = 1e-7;



figure(2);
semilogy(snr_range, Pb_bound);
% legend(legendStr, 'Location', 'Best', 'Interpreter', 'none');

end



function C = c_biawgn(sig)
    phi_sig_x = @(x,sig) 1/sqrt(8*pi*sig^2) *( exp(-(x+1).^2/(2*sig^2)) + exp(-(x-1).^2/(2*sig^2)) );
    x_range = linspace(-20*sig,20*sig,1e5);
    delta = x_range(2)-x_range(1);
    C = -delta* trapz( phi_sig_x(x_range, sig).*log2(phi_sig_x(x_range, sig)) ) -.5*log2(2*pi*exp(1)*sig^2);
end

function C = c_awgn(sig)
    C = .5*log2(1+1/sig^2);
end

