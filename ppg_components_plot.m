function ppg_components_plot
% ppg_components_plot creates a plot of the PPG signal illustrating the
% origins of light attenuation (arterial blood, venous blood and other
% tissues).
%
%               ppg_components_plot
%
%	This file creates an image adapted from:
%           Peter H Charlton. (2018, August). Capitalising on Smart
%           Wearables to Improve Health Monitoring. Zenodo.
%           DOI: http://doi.org/10.5281/zenodo.1406011  
%   Please cite this publication when using this image.
%   
%   Output:
%       image files in the same folder as this script
%     
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Licence:
%       please see the accompanying file named "LICENSE"
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function specifies that data will be downloaded to, and the image 
% files created in, the folder containing this script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
up = setup_up;

create_folders(up);

download_data(up);

extract_data(up);

plot_data(up);

end

function up = setup_up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The path of the folder in which to store the data and images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[filepath,~,~] = fileparts(mfilename('fullpath'));
up.paths.folders.root = [filepath, filesep];

% remaining paths
up.paths.folders.data = [up.paths.folders.root, 'raw_data', filesep];
up.paths.folders.plot = [up.paths.folders.root, 'plots', filesep];
up.paths.converted_data = [up.paths.folders.data, 'converted_data'];

% download paths
up.files.web = {'https://www.physionet.org/physiobank/database/bidmc/bidmc05'};
up.files.times.deb = [3];
up.files.times.fin = [8];
up.files.sigs.names = {{'ppg'}};
up.files.sigs.nos = {[2]};

close all

% Check that the WFDB Matlab Toolbox has been installed
if ~exist('getWfdbClass', 'file')
    error('Couldn''t find the WFDB Matlab Toolbox. Please install as described at the top of the file.')
end

% Check that the folder in which to store downloaded data exists
if ~exist(up.paths.folders.root, 'dir')
    error('Couldn''t find the folder in which to store the data.')
end

% The script will need to download data in order to run
do_download = 1;
if ~do_download
    error('The script will need to download data from PhysioNet in order to run')
end

end

function create_folders(up)

folders = fieldnames(up.paths.folders);
for folder_no = 1 : length(folders)
    eval(['curr_folder = up.paths.folders.' folders{folder_no} ';'])
    
    if ~exist(curr_folder)
        mkdir(curr_folder)
    end
   
end

end

function download_data(up)

exts = {'.hea', '.dat'};

for file_no = 1 : length(up.files.web)
    
    url = up.files.web{file_no};
    
    temp = strfind(url, '/');
    filename = url(temp(end)+1:end);
    filepath =  [up.paths.folders.data, filename];
    
    for ext_no = 1 : length(exts)
        curr_ext = exts{ext_no};
        expected_outfilename = [filepath,curr_ext];
        if ~exist(expected_outfilename, 'file')
            outfilename = websave([filepath,curr_ext],[url, curr_ext]);
        end
    end
    
end

end

function extract_data(up)

if exist([up.paths.converted_data, '.mat'])
    return
end

counter_no = 0;
for file_no = 1 : length(up.files.web)
    
    % Extract data from downloaded file
    url = up.files.web{file_no};
    temp = strfind(url, '/');
    filename = url(temp(end)+1:end);
    filepath =  [up.paths.folders.data, filename];
    cd(up.paths.folders.data)
    [signal,Fs,tm]=rdsamp(filename, [],[]);
    
    % Extract required data
    no_sigs = length(up.files.sigs.names{file_no});
    for sig_no = 1 : no_sigs
        curr_sig_no = up.files.sigs.nos{file_no}(sig_no);
        curr_times = [up.files.times.deb(file_no), up.files.times.fin(file_no)];
        rel_els = find(tm>= curr_times(1) & tm <= curr_times(2));
        curr_sig.v = signal(rel_els,curr_sig_no);
        curr_sig.t = [0:length(curr_sig.v)-1]./Fs;
        
        % store data
        counter_no = counter_no+1;
        data(counter_no).sig = curr_sig.v;
        data(counter_no).t = curr_sig.t;
        data(counter_no).db = strrep(url(temp(end-1)+1:temp(end)-1), 'db', '');
        data(counter_no).name = up.files.sigs.names{file_no}{sig_no};
    end
    
end

save(up.paths.converted_data, 'data')

end

function plot_data(up)


%% Make short PPG Plot

% setup
load(up.paths.converted_data)
ftsize = 24;
pos_long = [20,20,500,250];
pos_ppg = [20,20,600,600];
pos_short = [20,20,1000,470];

% Extract data
if exist('extractfield')
    temp = extractfield(data, 'name');
    sig_no = find(strcmp(temp, 'ppg')); clear temp
else
    sig_no = 1;
end
curr = data(sig_no);

curr.fix = linspace(0.9,1,length(curr.t));
curr.ven = curr.fix + 0.7 + 0.1*sin(curr.t*2*pi/4);
curr.sig = (curr.sig-min(curr.sig))/(max(curr.sig)-min(curr.sig));
curr.ppg = curr.ven + 0.5 + curr.sig(:)';

% detect peaks
S.v = curr.ppg;   % extract PPG data
S.fs = 125;
beat_detector = 'MSPTD';     % Select Incremental-Merge Segmentation beat detector
[peaks, onsets, mid_amps] = detect_ppg_beats(S, beat_detector);     % detect beats in PPG

%plot(curr.t, curr.fix), hold on, plot(curr.t, curr.ven), plot(curr.t, curr.ppg)

% - Make figure (PPG and points)
figure('Position', pos_short)

% truncate signal to only portion with detected beats
tol = 0.3; % in secs, time either side of labelled portion of signal to include in plot
rel_els = find(curr.t>=curr.t(onsets(1))-tol & curr.t<=curr.t(peaks(end))+tol);
curr.t = curr.t(rel_els);
curr.ppg = curr.ppg(rel_els);
curr.t = curr.t - curr.t(1);
peaks = peaks-rel_els(1)+1;
onsets = onsets-rel_els(1)+1;
mid_amps = mid_amps-rel_els(1)+1;

% plot
plot(curr.t, curr.ppg, 'k', 'LineWidth', 2)
hold on
% plot peaks, onsets, and mid-pts
rel_pts = peaks; do_line = false; curr_color = 'r';
plot_pts(curr, rel_pts, do_line, curr_color);
rel_pts = onsets; do_line = false; curr_color = 'k';
plot_pts(curr, rel_pts, do_line, curr_color);
rel_pts = mid_amps; do_line = true; curr_color = 'b';
plot_pts(curr, rel_pts, do_line, curr_color);

ylim([0 3.5])
xlims = xlim;
set(gca, 'FontSize', ftsize, 'YTick', [], 'XTick', xlims(1):xlims(2))
xlabel('Time (s)', 'FontSize', ftsize)
ylabel('PPG (au)', 'FontSize', ftsize)
box off

% annotate dc offset
rel_pt = mid_amps(end);
p1 = [curr.t(rel_pt) 0];
p2 = [curr.t(rel_pt) curr.ppg(rel_pt)];
labeltxt = 'DC amplitude'; linecolor = 'b'; before_log = false;
annotate_point(rel_pt, curr, ftsize, labeltxt, p1, p2, linecolor, before_log);

% annotate ac component
rel_pt = mid_amps(1);
p1 = [curr.t(onsets(1)) curr.ppg(onsets(1))];
p2 = [curr.t(onsets(1)) curr.ppg(peaks(1))];
labeltxt = 'AC amplitude'; linecolor = 'r'; before_log = true;
annotate_point(rel_pt, curr, ftsize, labeltxt, p1, p2, linecolor, before_log);

% legend
legend({'', 'peaks', 'onsets', 'mid-points', ''}, 'Location', 'southwest', 'FontSize', ftsize-2)

% save
plot_handle = gcf;
filename = 'short_raw_ppg_init_pts';
save_plot(plot_handle, filename, up)

close all

% - Make figure (just PPG)
figure('Position', pos_short)
plot(curr.t, curr.ppg, 'k', 'LineWidth', 2)
ylim([0 3.31])
set(gca, 'FontSize', ftsize, 'YTick', [])
xlabel('Time (s)', 'FontSize', ftsize)
title('Raw Photoplethysmogram', 'FontSize', ftsize)
box off
plot_handle = gcf;
filename = 'short_raw_ppg_init';
save_plot(plot_handle, filename, up)

close all

% - Make figure (components)

% fill in
figure('Position', pos_short)
h = area(curr.t, curr.ppg, 'LineStyle', 'none'); hold on
h.FaceColor = [1,0,0];
plot(curr.t, curr.ppg, 'k', 'LineWidth', 2)
%plot(curr.t, curr.ven, '--k', 'LineWidth', 2)
h = area(curr.t, curr.ven, 'LineStyle', 'none'); hold on
h.FaceColor = [0,0,1];
plot(curr.t, curr.ven, '--k', 'LineWidth', 2)
h = area(curr.t, curr.fix, 'LineStyle', 'none'); hold on
h.FaceColor = [0.8,0.8,0.8];
plot(curr.t, curr.fix, '--k', 'LineWidth', 2)
ylim([0 3.31])
set(gca, 'FontSize', ftsize, 'YTick', [])
xlabel('Time (s)', 'FontSize', ftsize)
ylab = ylabel('PPG', 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'position', get(ylab,'position')-[0.2,0.1,0]);
%title('Raw Photoplethysmogram', 'FontSize', ftsize)
box off

% annotate
dim = [.2 .19 .1 .1];
str = 'Other Tissues';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'BackgroundColor', [1,1,1], 'FontSize', ftsize, 'HorizontalAlignment', 'center');
dim = [.2 .40 .1 .1];
str = 'Venous Blood';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'BackgroundColor', [1,1,1], 'FontSize', ftsize, 'HorizontalAlignment', 'center');
dim = [.2 .58 .1 .1];
str = 'Arterial Blood';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'BackgroundColor', [1,1,1], 'FontSize', ftsize, 'HorizontalAlignment', 'center');

plot_handle = gcf;
filename = 'short_raw_ppg';
save_plot(plot_handle, filename, up)
close all

end

function annotate_point(rel_pt, curr, ftsize, labeltxt,p1,p2, linecolor, before_log)
dp = p2-p1;
quiver(p1(1),p1(2),dp(1),dp(2),0,'Color', linecolor, 'LineWidth', 2)
if before_log
    x_offset = -0.1;
    y_val = p1(2)-0.1*p1(2);
else
    x_offset = 0.1;
    y_val = mean([p1(2), p2(2)]);
end
text(p1(1)+x_offset,y_val, labeltxt, 'Color', linecolor, 'FontSize', ftsize)
end

function plot_pts(curr, rel_pts, do_line, curr_color)
if nargin<3
    do_line = false;
end
if do_line
    markertype = 'o-';
else
    markertype = 'o';
end
plot(curr.t(rel_pts), curr.ppg(rel_pts), markertype, 'Color', curr_color, 'MarkerSize', 6, 'LineWidth', 2)
end

function save_plot(plot_handle, filename, up)

save_path = [up.paths.folders.plot, filename];
print(plot_handle, save_path, '-depsc')
print(plot_handle, save_path, '-dpng')
print(plot_handle, '-dsvg', save_path)
fid = fopen([save_path, '.txt'], 'w');
fprintf(fid, ['Created using ' mfilename, ', ', date]);
fclose(fid);

end