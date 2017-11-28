function plot_de_sims(plotdir, varargin)

RESULTS_SUFFIX = 'minLUT0_iters_*bits.txt';

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
c = {'bs-', 'ro-', 'mv-', 'c^-', 'g<-', 'k>-'};

% Define linewidth
lineWidth = 2;


if isempty(varargin)
    D = dir( [plotdir '/' 'rate' '*' RESULTS_SUFFIX]);
    noSim = length(D);
    for dd = 1:noSim
        result_names{dd} = D(dd).name;
    end
else
    result_names = varargin{1};
    % Find how many simulations we have
    noSim = length(result_names);
end

ii_res = 1;
ii_iters = 1;
legendStr_res={};
legendStr_iters={};
for ii=1:noSim
    % Read file:
    fid = fopen([plotdir '/' result_names{ii}],'rt');
    C = textscan(fid,'%s','Delimiter','','endofline','');
    C = C{1}{1};
    fclose(fid);

    % Get thresholds
    thresholds = regexp(C,'Threshold\(s\) found = \[(.*)\]', 'tokens', 'dotexceptnewline');
    % Prior to commit fda5137341699ae0784ce89fa28a25203b162d8a, result
    % files were inconsistent.
    if(size(thresholds)==0)
        thresholds = regexp(C,'Thresholds = \[(.*)\]', 'tokens', 'dotexceptnewline');
    end
    if(size(thresholds)==0)
        thresholds = regexp(C,'Threshold found =\s([-+]?[0-9]*\.?[0-9]*)\s', 'tokens', 'dotexceptnewline');
    end
%     result_names{ii}
    thresholds = str2num(thresholds{1}{1});

    % Get iterations
    iters = regexp(C,'Maximum Number of message passing iterations = \[(.*)\]', 'tokens', 'dotexceptnewline');
    
    if(size(iters)==0)
        iters = regexp(C,'Maximum Number of message passing iterations =\s(\d+)\s', 'tokens', 'dotexceptnewline');
    end
    iters = str2num(iters{1}{1});
    
    % Get resolutions
    res = regexp(C,'Resolutions \[channel bits, message bits; \.\.\.\] = \[\[(.*)\]\]', 'tokens');
    if(size(res)~=0)
        res = str2num( [ '[' res{1}{1} ']' ]).';
    else
        res = regexp(C,'Resolution of discrete pmfs =\s(\d+)\sbit', 'tokens', 'dotexceptnewline');
        res = str2num(res{1}{1});
    end
    
    % Plot results 
    if(length(res) == length(thresholds) && length(thresholds)>1)
        figure(1)
        plot(res(2,:), thresholds, c{mod(ii,length(c))+1}, 'LineWidth', lineWidth)
        legendStr_res{ii_res} = result_names{ii};
        ii_res = ii_res +1;
        hold on
    elseif(length(iters) == length(thresholds) && length(thresholds)>1)
        figure(2)
        plot(iters, thresholds, c{mod(ii,length(c))+1}, 'LineWidth', lineWidth)
        legendStr_iters{ii_iters} = result_names{ii};
        ii_iters = ii_iters + 1;
        hold on
    else
        warning([plotdir '/' result_names{ii} ' seems to neither contain a resolution nor an iteration sweep. Skipping...'])
        continue;
    end
     
    
end

figure(1)
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('Resolution', 'FontSize', fontSize)
ylabel('BIAWGN threshold', 'FontSize', fontSize)
grid on
legend(legendStr_res, 'Location', 'Best', 'Interpreter', 'none');

figure(2)
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('# Iterations', 'FontSize', fontSize)
ylabel('BIAWGN threshold', 'FontSize', fontSize)
grid on
legend(legendStr_iters, 'Location', 'Best', 'Interpreter', 'none');
end