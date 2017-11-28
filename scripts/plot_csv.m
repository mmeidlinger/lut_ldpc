function  plot_csv(csv_dir, suffix, varargin)


% If nothing is parsed plot everything in csv_dir folder
if isempty(varargin)
    D = dir( [csv_dir '/' '*' suffix '.csv']);
    num_files = length(D);
    for dd = 1:num_files
        csv_names{dd} = D(dd).name;
    end
else
    csv_names = varargin{1};
    % Find how many simulations we have
    num_files = length(csv_names);
end

% Close previous plots
close all

% Change default axes fonts.
fontSize = 12;
set(0,'DefaultAxesFontName', 'Times')
set(0,'DefaultAxesFontSize', fontSize)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times')
set(0,'DefaultTextFontSize', fontSize)

% Define linewidth
lineWidth = 2;


% Define plot colors
c = {'bs-', 'ro-', 'mv-', 'c^-', 'g<-', 'k>-'};

for ii = 1:num_files
    BER_data = csvread(csv_names{ii});
    figure(1)
    semilogy(BER_data(:,1), BER_data(:,2), c{mod(ii,length(c))+1}, 'LineWidth', lineWidth)
    legendStr{ii} = csv_names{ii};
    hold on

end

% Beautify plots
figure(1)
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('Eb/N0 (dB)', 'FontSize', fontSize)
ylabel(suffix , 'FontSize', fontSize);
grid on


figure(1);
legend(legendStr, 'Location', 'Best', 'Interpreter', 'none');

end

