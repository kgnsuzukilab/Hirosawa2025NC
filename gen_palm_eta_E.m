function gen_palm_eta
% Matlab eta version of PALM video image generation program gen-palm
% TRead the CSV file output by ThunderSTORM and reconstruct the temporal variation of the structure.
% In this case, we assume a Gaussian distribution on the spatial axis and a uniform distribution
% on the temporal axis.
% This program uses the common calculation engine gen_palm_zeta.m, as DoCAnalysis.
% Parallel computing : use. 
% GPU : no use.
% 2021.9.11 Coded by Y.Yokota 

% Parameters and Condition Setting
PALMf = 5; % PALM factor
D = 6; % Input number of frames to shift PALM image
W = 1002; % Input number of frames for each PALM image
R = 10; % nm / pixel for PALM image
FR = 30; % Frame rate
A = 8; % Parameters to determine kernel width
flag1 = true; % If true, save the generated PALM image as a mat file
flag3 = true; % If true, save the generated PALM image as an avi movie file
flag4 = true;  % If true, binarize the generated PALM image and save it as an avi movie file.

% Specify and read CSV files output by ThunderSTORM
[csvfilename,csvpath] = uigetfile('*.csv',' Specify csv file output by ThunderSTORM.'); % GUI use
% csvpath = 'E:\PALM\project5\210712 Test movie\CountSpotInStructure\sample1\'; % direct specification
% csvfilename = 'hc6crop-188.csv';  % direct specification
Tbl = readtable(fullfile(csvpath,csvfilename),'PreserveVariableNames',true);

% Specify the video file to be analyzed by ThunderSTORM to determine the size of the video image
[originalavifilename,originalavipath] = uigetfile(fullfile(csvpath,'*.avi'),' Specify the video file to be analyzed by ThunderSTORM.'); % GUI use
% originalavipath = 'E:\PALM\project5\210712 Test movie\CountSpotInStructure\sample1\'; % direct specification
% originalavifilename = 'hc6crop-188.avi'; % direct specification
v = VideoReader(fullfile(originalavipath,originalavifilename));
Nx = PALMf * v.Width; % Horizontal size of PALM video image
Ny = PALMf * v.Height; % Vertical size of PALM video image
Nframe = v.NumFrames; % Frame size of PALM video image

% Specify the folder where the PALM images are to be saved
resultpathname = csvpath; % Save in the same folder as the loading folder
% resultpathname = uigetdir('',' Specify the folder where the PALM images are to be saved');

% Read the csv file output by ThunderSTORM and extract the positional information of the molecules
t1 = Tbl.('frame'); % Time (frame)
x1 = Tbl.('x [nm]')/R+1; % x coordinate (Consider upper left as coordinate (0,0))
y1 = Tbl.('y [nm]')/R+1; % y coordinate (Consider upper left as coordinate (0,0))
s1 = Tbl.('uncertainty [nm]')/R;
spot = single([x1 y1 t1 s1]);
clear Tbl

% The function gen_palm_zeta is called to create a palm movie
[img2,spot,tse] = gen_palm_zeta(spot,Nx,Ny,Nframe,A,D,W);
img1 = img2(:,:,1);
Th = 1;
Nt = size(img2,3); % The number of frames in PALM video image

% Save result as a mat file
if flag1
    save([resultpathname '\PALM.mat'],'img1','img2','Th','spot','FR','PALMf','tse','R','A','W','D','Nx','Ny','Nframe','Nt','-v7.3')
end

% Save grayscale PALM video image as avi file
if flag3
    % Integer to the range [0, 255]
    ind = round(linspace(1,Nt,round(Nt/50)));
    img3 = img2(:,:,ind);
    [h,v] = histcounts(img3(:),200);
    r = v(find(cumsum(h)/sum(h)<0.9999,1,'last'));
    r = 255/r;
    img3 = uint8(img2*r);
    v = VideoWriter([resultpathname '\palm.avi'],'Grayscale AVI');
    v.FrameRate = FR;
    open(v);
    writeVideo(v,img3);
    close(v);
end

% Binarize PALM video using the binary method and save the result as an avi file
if flag4
    img3 = uint8(img2>1)*255;
    vb = VideoWriter([resultpathname '\palmbw.avi'],'Grayscale AVI');
    vb.FrameRate = FR;
    open(vb);
    writeVideo(vb,img3);
    close(vb);
    disp('Check the binarized movie palmbw.avi and if the binarization threshold is inappropriate,')
    disp('run setthresholdmovie.mlapp to manually set the threshold.')
end
