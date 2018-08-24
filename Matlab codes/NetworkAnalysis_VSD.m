%%
% % %  Set up paths and empty matrices to read in and process data % % % 

% enter the folder where the laser stim data is stored
the_folder='C:\Users\PsychiatryLab\Desktop\Diana\Neurophotonics\For website\Data\';

% list the tif files associated with the first rep of photostimulation
the_stim_list={'M2L1.tif', ... %Secondary Motor Cortex - Left Hem.
    'BCR1.tif', ... %Primary Barrel Cortex - Right Hem.
    'RSL1.tif', ... %Retrosplenial Cortex - Left Hem.
    'MFR1.tif', ... %Forelimb Associated Motor Cortex - Right Hem.
    'BCL1.tif', ... %Primary Barrel Cortex - Left Hem.
    'HLR1.tif', ... %Primary Hind limb Cortex - Right Hem.
    'MBL1.tif', ... %Whisker assoicated motor Cortex - Left Hem.
    'V2R1.tif', ... %Secondary visual cortex - Right Hem.
    'FLL1.tif', ... %Forelimb cortex - Left Hem.
    'PTR1.tif', ... %Parietal association cortex - Right Hem.
    'HLL1.tif', ... %Hindlimb cortex - Left Hem.
    'M2R1.tif', ... %Secondary Motor Cortex - Right Hem.
    'PTL1.tif', ... %Parietal association cortex - Left Hem.
    'FLR1.tif', ... %Forelimb cortex - Right Hem.
    'V2L1.tif', ... %Secondary visual cortex - Left Hem.
    'MBR1.tif', ... %Whisker assoicated motor Cortex - Right Hem.
    'V1L1.tif', ... %Primary Visual cortex - Left Hem.
    'V1R1.tif', ... %Primary Visual cortex - Right Hem.
    'MFL1.tif', ... %Forelimb Associated Motor Cortex - Left Hem.
    'RSR1.tif'}; %Retrosplenial Cortex - Right Hem. 

% every two stims has a no stim file for normalization
the_no_list={'NOA1.tif','NOB1.tif','NOC1.tif','NOD1.tif'...
    'NOE1.tif','NOF1.tif','NOG1.tif','NOH1.tif','NOI1.tif','NOJ1.tif'};

% How many trials per site?
nreps=2;

% which time point is the stim?
tstim=31;

% Determine the number of stims we have
nstim=numel(the_stim_list);

% How many frames? 
framenum = 108; 

% read in the first file to get the dimensions
the_file=[the_folder the_stim_list{1}];
tmp=imreadalltiff(the_file,framenum);
w=size(tmp,1); %width
l=size(tmp,2); %length
% framenum = size(tmp,3); 

s = [w,l,framenum]; %3D size
ss = [w,l,framenum,nreps]; %4D size
clear tmp

%% SECTION A
% % %  Create averaged VSD responses (dF/F0) % % % 

% make a matrix to store all the data
img_01 = zeros(ss); img_no = zeros(ss);
img_01_av = zeros(s); img_no_av = zeros(s);
img_norm = zeros(ss);
df_stim = zeros(w,l,framenum,nreps,nstim); % store df per stim site
baseline_from = 3; baseline_to = 28; % choose frames to take baseline from

% Read in files
for i = 1:nstim
    for j = 1:nreps
        % read in the stim file
        if j == 1
            the_file=[the_folder the_stim_list{i}];
        else
            prefix=the_stim_list{i};
            prefix=prefix(1:end-5);
            suffix=the_stim_list{i};
            suffix=suffix(end-3:end);
            the_file=[the_folder prefix num2str(j) suffix];
        end
        img_01(:,:,:,j)=imreadalltiff(the_file,108);
        
        % read in the associated NO trial.
        the_no_file=[the_folder the_no_list{floor((j-1)/2)+1}];
        img_no(:,:,:,j)=imreadalltiff(the_no_file,108);
        
        % divide stim trial by associated no trial to reduce bleaching
        img_norm(:,:,:,j) = img_01(:,:,:,j)./img_no(:,:,:,j);
       
        % isolate baseline and calculate df/F0 
        baseline_01=((mean(img_norm(:,:,baseline_from:baseline_to,j),3)));
        baseline_01= repmat(baseline_01,1,framenum); %divide by same baseline for each frame
        baseline_01 = reshape(baseline_01, [128 128 108]); 
        img_norm(:,:,:,j)=((img_norm(:,:,:,j)-baseline_01)./baseline_01)*100; % calculate df/F0 as percent
        
        % optional filtering step 
        sigma = 2.5; filtsize = 5; % set filter type and size
        H = fspecial('gaussian',filtsize,sigma);
        img_gauss = zeros(ss);
        for jj = 1:framenum
            img_gauss(:,:,jj,j) = imfilter(squeeze(img_norm(:,:,jj,j)),H);
        end
        img_norm = img_gauss;
        
    end
    
    df_stim(:,:,:,:,i) = img_norm; % df_stim[w.l,framenum,nreps,nstim]

end

% Average VSD response over number of reps    
img_av = squeeze(mean(df_stim,4)); 
    
%% SECTION B                
% % % Calculate the mean VSD response from each region of interest (ROI) at one time point % % % 

% read in the ROI positions (pixels)
load('pos.mat')
npos=numel(pos)/2;
roisize=5;

% make a matrix to store all the data
dat_roi = zeros(npos,framenum); 
resp = zeros(npos,1);
resp_mat = zeros(npos,nstim); 

% Show ROI positions on masked image (for visualization purposes only)
m = imreadalltiff('Mask.tif'); m = double(m); 
figure;imshow((img_av(:,:,1,1).*m),[]);axis off; colormap('gray');hold on;
title('ROI placement'); 
for i = 1:npos
    rectangle('Position',[pos(i,1),pos(i,2),roisize,roisize],'LineWidth',1,'EdgeColor','r');
end


% calculate average VSD response at each ROI position
for ii = 1:nstim %number of unique dF/F0 matrices 
    img = img_av(:,:,:,ii); 
    
    % delete frame which has the stimulus artifact
    img(:,:,32) = []; 
    img(:,:,108) = img(:,:,107);
    
    for i = 1:framenum % number of frames 
        for j = 1:npos % number of VSD response sites (ROIs)
            dat_roi(j,i) = MeanFromROI(img(:,:,i),pos(j,1),pos(j,2),roisize);
        end
    end
 
    for j = 1:npos
        base_av = mean(dat_roi(j,20:28)); %calculate baseline on these frames only
        the_sd= std(dat_roi(j,20:28));
        the_max=max(dat_roi(j,(tstim+1):(tstim+11))); %look for max response within these frames 
        
        % Threshold to avoid capturing noise as a response
        if the_max > (the_sd * 2.5)
            resp(j,1) = (sum(dat_roi(j,33:35))) - (base_av*length(33:35)); % Take response from sum of first three frames, approx 20ms
        else
            resp(j,1) = 0;
        end
        
        % Replace any negative values with zeros
        if resp(j,1) < 0
            resp(j,1) = 0;
        end
                
    end
    
    
    % store responses in matrix
    resp_mat(:,ii) = resp(:); 
    
end

% remove diagonal (do not take responses from the stimulation site)
for k = 1:20
    resp_mat(k,k) = 0;
end




%% SECTION C 
% % % Create the image of the connectivity matrix % % %

% change the order of ROIs (for visualization purposes only)
the_order=[1, 7, 19, 9, 11, 5, 13, 3, 15, 17, 12, 16, 4, 14, 6, 2, 10, 20, 8, 18]; % NB: this order is roughly anterior to posterior, left hemisphere then right hemisphere
pos_order = pos(the_order,:); 
resp_mat = resp_mat(the_order,the_order);
figure;imagesc(resp_mat)
% set labels for figure
ax = {'M2L', 'MBL','MFL','FLL','HLL','BCL','PTL','RSL','V2L','V1L','M2R','MBR','MFR','FLR','HLR', 'BCR','PTR','RSR','V2R','V1R'}; % x is the point of stimulation
set(gca,'xtick',1:nstim, 'xticklabel', ax); set(gca,'ytick',1:npos,'yticklabel', ax);
xlabel ('Stimulation site'); ylabel('Response site');
title('Connectivity Matrix'); 


%% SECTION D
% % % % % Determine threshold for network diagram % % %
% % % NB: This sections requires functions from the brain connectivity
% % % toolbox (Rubinov and Sporns, 2010)

% Determine multiple threshold levels to test
max_resp = max(max(resp_mat)); 
thr_max=(max_resp.*0.8); %Set boundaries for thresholds to try
n_thr= 30; %number of thresholds to test
thr_mat=linspace(0,thr_max,n_thr); %Create 30 thresholds to try 

% Display number of connections as a function of threshold level
num_conn = zeros(1,20); 
for i = 1:n_thr
    thr_abs = threshold_absolute(resp_mat,thr_mat(i)); % NB: threshold_absolute from brain connectivity toolbox
    num_conn(i) = length(nonzeros(thr_abs));
end

% Calculate efficiency and characteristic path length (L) as a function of threshold
% NB: Distance_wei(for a weighted matrix)- function from brain connectivity toolbox
% NB: Charpath function from brain connectivity toolbox. 
%   Outputs: lambda (characteristic path length); efficiency (global); ecc
%   (eccentricity for each vertix); radius (of graph); diameter (of graph)

lambda_mat = zeros(1,20);
eff_mat = zeros(1,20);
for i=1:length(thr_mat)
    thr_abs = threshold_absolute(resp_mat, thr_mat(i));
    D=distance_wei(1./thr_abs); 
    [lambda,eff,ecc,radius,diameter] = charpath(D);
%     lambda_mat(i) = lambda; % characteristic path length 
    eff_mat(i) = eff; % inverse shortest path
end
% %  Uncomment to display plot of characteristic path length
% figure; scatter(thr_mat, lambda_mat); xlabel('threshold'); ylabel('char path length'); 


% Display plots for number of connections and global efficiency.
% The vertical line indicates a threshold level of 0.2577, which resulted
% in a 5% decrease in global efficiency and removed 31.6% of the
% connections
figure; ax1 = subplot(2,1,1); 
scatter(ax1,thr_mat,num_conn); ylabel('Number of Connections'); xlabel('threshold');
axis([0 max(thr_mat) 0 max(num_conn*1.1)]); % Set axes
line([0.2577 0.2577],[0 (max(num_conn)*1.1)],'Marker','.','LineStyle','-')% Show threshold level at 0.2577
ax2 = subplot(2,1,2);
scatter(ax2,thr_mat,eff_mat); ylabel('Efficiency'); xlabel('threshold');
axis([0 max(thr_mat) 0 max(eff_mat*1.1)]); % Set axes
line([0.2577 0.2577],[0 (max(eff_mat)*1.1)],'Marker','.','LineStyle','-')

% % Optional: Calculate clustering coefficient (C) as a function of threshold
% % NB: clustering_coef_wd (for a weighted matrix)- function from brain connectivity toolbox
% C_mat = zeros(20,20);
% for i=1:length(thr_mat)
%     thr_abs = threshold_absolute(resp_mat, thr_mat(i));
%     C=clustering_coef_wd(thr_abs); 
%     C_mat(:,i) = C;  
% end
% figure; scatter(thr_mat, mean(C_mat,1)); xlabel('threshold'); ylabel('clustering coefficient');

clear('thr_abs'); 



%% SECTION E
% % %  Create the image of the network diagram using bioconnectivity toolbox % % %

thresh = 0.2577; % Set threshold for network diagram based on efficiency in previous steps
thr_abs = threshold_absolute(resp_mat, thresh); % apply threshold to reponse matrix
bg=biograph(thr_abs);
bg.NodeAutoSize='off';
dolayout(bg)

% % % % % Set custom node properties % % % % % 
% % Node color
num_nodes = nstim; 
node_color = zeros(3,num_nodes);
% Left hemisphere nodes
node_color(:,1) = [255 0 255]; %M2 - pink
node_color(:,2) = [255 204 51]; %MB - light orange
node_color(:,3) = [0 136 255]; %MF - light blue
node_color(:,4) = [46 49 146]; %FLS1 - dark blue
node_color(:,5) = [209 211 212]; %HLS1 - grey
node_color(:,6) = [245 130 32]; %BCS1 - orange
node_color(:,7) = [236 28 36]; %PTA - red
node_color(:,8) = [202 108 56]; %RS - brown
node_color(:,9) = [0 227 0]; %V2 - light green
node_color(:,10) = [0 153 0]; %V1 -  green
% Right hemisphere nodes
node_color(:,11) = [255 0 255]; %M2
node_color(:,12) = [255 204 51]; %MB
node_color(:,13) = [0 136 255]; %MF
node_color(:,14) = [46 49 146]; %FLS1
node_color(:,15) = [209 211 212]; %HLS1
node_color(:,16) = [245 130 32]; %BCS1
node_color(:,17) = [236 28 36]; %PTA
node_color(:,18) = [202 108 56]; %RS
node_color(:,19) = [0 227 0]; %V2M
node_color(:,20) = [0 153 0]; %V1
node_color(:,:) = node_color./255;

% Node names
bg.nodes(1).ID='M2L';
bg.nodes(2).ID='MBL';
bg.nodes(3).ID='MFL';
bg.nodes(4).ID='FLL';
bg.nodes(5).ID='HLL';
bg.nodes(6).ID='BCL';
bg.nodes(7).ID='PTL';
bg.nodes(8).ID='RSL';
bg.nodes(9).ID='V2L';
bg.nodes(10).ID='V1L';
bg.nodes(11).ID='M2R';
bg.nodes(12).ID='MBR';
bg.nodes(13).ID='MFR';
bg.nodes(14).ID='FLR';
bg.nodes(15).ID='HLR';
bg.nodes(16).ID='BCR';
bg.nodes(17).ID='PTR';
bg.nodes(18).ID='RSR';
bg.nodes(19).ID='V2R';
bg.nodes(20).ID='V1R';

% Node shape
for i=1:num_nodes
    bg.nodes(i).Shape='circle';
    bg.nodes(i).color= node_color(:,i);
end

% Node position - based on pos.mat (functional and stereotaxic coordinates of nodes)
% This positions the nodes in approximate anatomical position but has no
% effect on network properties 
coor_tmp((1:10),:) = pos_order(1:10,:); 
coor_tmp((11:20),:) = pos_order(11:20,:);
x=coor_tmp(:,1); y=coor_tmp(:,2)*-1; %invert Y axis
x=(x-min(x))/(max(x)-min(x))*350; y=(y-min(y))/(max(y)-min(y))*350; % rescale for biograph window

for i=1:num_nodes
    tmp=[x(i) y(i)];
    bg.nodes(i).Position = tmp;
end
dolayout(bg, 'Pathsonly', true)
% view(bg)

% Node size
node_weight=sum(thr_abs,2); % calculate node out strength
max_node_weight=max(node_weight);
node_weight(~node_weight) = inf; % replace any zeros w/ inf to calculate nonzero min
min_node_weight=min(node_weight);
node_weight(isinf(node_weight)) = min_node_weight; % replace infs with minimum node weight

% Calculate slope and intercept for linear function
min_size=15.; max_size=40.;
m=(max_size-min_size)/(max_node_weight-min_node_weight);
b= max_size-m*max_node_weight;

% Use linear function to set node size proporational to node weight
for i=1:num_nodes
      tmp = [b+m*node_weight(i), b+m*node_weight(i)];  % linear function F(x)= mx+b
      bg.nodes(i).Size = tmp; 
end

% Font size (also using linear function)
min_fsize=4; max_fsize=12.;
m=(max_fsize-min_fsize)/(max_node_weight-min_node_weight);
b= max_fsize-m*max_node_weight;

%  Set font size proportional to the node weight.
for i=1:num_nodes
    tmp = round(b+m*node_weight(i)); %Needs full integers
    bg.nodes(i).FontSize=tmp;
end


% % Line/edge width
min_lwidth=1.; max_lwidth=8.;
max_w=max(max(resp_mat));
tmp = resp_mat; tmp(~tmp) = inf; % replace any zeros w/ inf to calculate nonzero min
min_w=min(min(tmp));

m=(max_lwidth-min_lwidth)/(max_w-min_w);
b=max_lwidth-m*max_w;

% Set line widths proportional to edge weights - use power to increase the
% range of widths
for i=1:length(bg.edges)
    bg.edges(i).LineWidth=min_w+m*(((bg.edges(i).Weight)/b)^2);
end


% Set edges to match node color
for i=1:num_nodes
    for j=1:num_nodes
        path = [i j];
        edges = getedgesbynodeid(bg,get(bg.Nodes(path),'ID'));
        set(edges,'LineColor',node_color(:,i)) 
    end
end

dolayout(bg, 'Pathsonly', true)
view(bg)


