function wintld_gamma
% Matlab gamma version of WinTLD
% The PALM movie (mat file) output by gen_palm_delta.m is read to track the motion of the structure.
% In addition, read the trajectory data of each particle obtained by WinATR.
% The distance between the trajectory of each particle and the nearest structure is calculated.
% Parallel computing : use. 
% GPU : no use.
% 2019.7.2 Coded by Y.Yokota 
% 2020.1.23 Added a program (flag3) to output a table of co-localization distances by Isogai
% 2020.1.29 Modified wintld.subfuc2 by Yokota
% 2020.2.1 Added co-localization analysis of random spt and PALM images.
%          Changed to measure the distance between all structures and spt by Isogai

% Parameters and Condition Setting
FR = 30; % Frame rate of videos
BinSize = 5; % Histogram bin width of the distances (pixel) 
BinStart = -50; % The minimal value of the histogram (pixel) 
flag1 = false; % If true, save structure tracking results as text file
flag2 = false; % If true, graph the distance to each spt and proximity structure
flag3 = true; % If true, make a frequency distribution table of the distance between each SPT
              % and the neighboring structure and the probability of existence of each spt.
flag4 = true; % If true, make a pseudo color movie of all structures and all SPT particles
              % and save it to an avi movie file
flag5 = false; % If true, make a pseudo-color movie of the structure and only the SPT particles
               % in proximity to it and save it in an avi file.
               % To do this, the proximity must be defined and the entire analysis must be performed.
flag6 = false; % If true, in addition to option flag5, write the ID of the spt particles and structures
               % into the video.

RndN = 100; % The number of random spt generated per frame
flagR3 = true; % If true, make a frequency distribution table of the distance
               % between the random spt and the nearby structure
               % and the probability of the existence of the random spt
flagR4 = false; % If true, all structures and all random spt are made into a pseudo-color video
                % and saved in an avi file

% Loads images of numbers to write IDs in the video image.
S = load('suuji','suuji');
suuji = S.suuji;
clear S

% Specify the PALM movie (mat file) output by gen_palm_gamma.m
[PALMfilename,PALMpath] = uigetfile('*.mat','Specify the PALM movie (mat file) output by gen_palm_gamma.m'); % GUI use
% PALMpath = '..\project2\sample data\PALM images\'; % direct specification
% PALMfilename = 'PALM.mat';

% Specify the folder where the SPT trajectory saved by WinATR
sptpathname = uigetdir('','Specify the folder where the SPT trajectory saved by WinATR'); % GUI use
% sptpathname = '..\project2\sample data\spt'; % direct specification
files = dir([sptpathname '\*.dat']);

% Specify the foler where results are saved
resultspathname = uigetdir('','Specify the folder where results are saved'); % GUI use
% resultspathname = '..\project2\sample data\results of wintld'; % direct specification

% load PALM movie
load(fullfile(PALMpath,PALMfilename),'img2');

% Structure information acquisition and spt co-localization analysis
Nx = size(img2,2);
Ny = size(img2,1);
T = size(img2,3);
q = cell(length(files),1);
M = length(files); % The number of SPT particles
for n=1:M
    str = fileread(fullfile(files(n).folder,files(n).name));
    str = splitlines(str);
    q{n} = zeros(length(str),3);
    for k=1:length(str)
        temp = sscanf(str{k},'%f');
        if isempty(temp)
            q{n}(k:end,:) = [];
            break
        end
        q{n}(k,:) = temp(1:3)';
    end
    q{n}(:,1) = q{n}(:,1)+1; % Because WinATR output starts with the first frame at 0
    q{n}(:,2:3) = q{n}(:,2:3) * 5+1; % To match the resolution of PALM, increase it by a factor of 5
end

% Binarization of PALM video
load(fullfile(PALMpath,PALMfilename),'Th');
bw1 = img2>Th;
clear img2

% Identification of each structure in a PALM movie
C = bwconncomp(bw1);
N = C.NumObjects;
clear bw1

% To improve memory efficiency during parallel processing
PixelIdxList = cell(N,1);
for n=1:N
    PixelIdxList{n} = C.PixelIdxList{n};
end

% Examine the distance between each particle and each structure that is closer than 50 pixels
clear r
parfor k=1:N % for each structure
    for n=1:M % for each SPT particle
        r{n,k} = wintld_beta_subfunc2(q{n},PixelIdxList{k},Nx,Ny,T);
    end
end

% Examine the structure in close proximity to each SPT.
s = false(size(r));
M = size(r,1);
N = size(r,2);
parfor n=1:M % for each SPT particle
    for k=1:N % for each structure
        if isempty(r{n,k})
            s(n,k) = false;
        else
            s(n,k) = true;
        end
    end
end

save([resultspathname '\wintld.mat'],'q','C','r','s','FR')

% The information of each structure (frame number, x,y coordinates of the center of gravity,
% area, major axis length, minor axis length) is saved in a text file.When an region is divided
% into two or more areas, a line corresponding to the number of regions is written 
% for each frame number. The unit is pixel. The decimal point is shown to two decimal places, 
% but the significant digits are above the decimal point.
if flag1
    parfor n=1:N
        filename1 = [resultspathname '\PALM_' dec2base(n,10,4) '.dat'];
        fid = fopen(filename1,'w');
        [i1,i2,i3] = ind2sub([Ny Nx T],PixelIdxList{n});
        for k = i3(1):i3(end)
            ind = i3==k;
            bw2 = full(sparse(i1(ind),i2(ind),true,Ny,Nx));
            stats = regionprops(bw2,'Area','Centroid','MajorAxisLength','MinorAxisLength');
            for j=1:length(stats)
                fprintf(fid,'%5d %7.2f %7.2f %7.2f %7.2f %7.2f\n',k-1,stats(j).Centroid(1)-1,...
                    stats(j).Centroid(2)-1,stats(j).Area,stats(j).MajorAxisLength,stats(j).MinorAxisLength);
            end
        end
        fclose(fid);
    end
end

% Draw and store the change in distance between each SPT and the neighboring structure
if flag2
    for n=1:M % for each spt particle
        if ~any(s(n,:))
            continue
        end
        for k=1:N % for each structure
            if ~s(n,k)
                continue
            end
            ind = ~isnan(r{n,k}(:,2));
            plot(r{n,k}(ind,1),r{n,k}(ind,2),'k.-')
            text(min(r{n,k}(ind,1))+5,max(r{n,k}(ind,2))+2,num2str(k))
            xlim([1 T])
            ylim([-50 50])
            hold on
        end
        hold off
        title(['Distance to the approaching structure for SPT No.' num2str(n)])
        ylabel('Distance (pixel)')
        set(gcf,'Position',[10 50 600 500])
        drawnow
        [~,name1,~] = fileparts(files(n).name);
        outputfilename = fullfile(resultspathname,[name1 '.jpg']);
        print('-djpeg','-r300',outputfilename);
    end
    close all
end

% The distance between each SPT and the neighboring structure and the probability of existence
% of each spt are stored in a frequency distribution table.
if flag3
    tspt = [[],[]];
    for n=1:M % for each spt particle
        if ~any(s(n,:))
            continue
        end
        tsptk = [[],[]];
        for k=1:N % for each structure
            if ~s(n,k)
                continue
            end
            ind = ~isnan(r{n,k}(:,2));
            tsptn = [r{n,k}(ind,1),r{n,k}(ind,2)];
            tsptk = [tsptk; tsptn];
        end
        tspt = [tspt; tsptk];
    end
    outtspt = Hist_table(tspt(:,2), BinSize, BinStart);
    outputouttsptname = fullfile(resultspathname, 'Raw_ColoDistance_Hist.xls');
    writecell(outtspt,outputouttsptname)
    close all
end

% All structures and all SPT particles in pseudo color video and saved in avi file
if flag4
    L = labelmatrix(C);
    cmap = [zeros(1,3,'uint8'); uint8(rand(N,3)*255.5); repmat(uint8(255),1,3)];
    parfor k=1:size(L,3)
        for n=1:M % for each SPT particle
            idx = find(q{n}(:,1)==k,1,'first'); % SPT particles present at time k
            if isempty(idx)
                continue
            end
            L(:,:,k) = addspt(L(:,:,k),q{n}(idx,2:3),10,N+1);
        end
    end
    v = VideoWriter([resultspathname '\wintld1.avi']);
    v.Quality = 100;
    v.FrameRate = FR;
    open(v);
    for k=1:ceil(T/100)
        img = ind2rgb_y(L(:,:,(k-1)*100+1:min(k*100,T)),cmap);
        writeVideo(v,img);
    end
    close(v)
end

% Make a pseudo-color movie of the structure and only the SPT particles
% in proximity to it and save it in an avi file.
% To do this, the proximity must be defined and the entire analysis must be performed.
if flag5
    ind_spt = find(any(s,2)'); % SPT particle ID of interest
    Nspt = length(ind_spt);
    ind_str = find(any(s,1)); % structure ID of interest
    Nstr = length(ind_str);
    CC = C;
    CC.PixelIdxList = C.PixelIdxList(ind_str);
    CC.NumObjects = Nstr;
    L = labelmatrix(CC);
    cmap = [zeros(1,3,'uint8'); uint8(rand(Nstr,3)*255.5); repmat(uint8(255),1,3)];
    parfor k=1:size(L,3)
        for n=ind_spt % for each spt particle
            idx = find(q{n}(:,1)==k,1,'first'); % SPT particles present at time k
            if isempty(idx)
                continue
            end
            L(:,:,k) = addspt(L(:,:,k),q{n}(idx,2:3),10,Nstr+1);
        end
    end
    v = VideoWriter([resultspathname '\wintld2.avi']);
    v.Quality = 100;
    v.FrameRate = FR;
    open(v);
    for k=1:ceil(T/100)
        img = ind2rgb_y(L(:,:,(k-1)*100+1:min(k*100,T)),cmap);
        writeVideo(v,img);
    end
    close(v)
end

% If true, in addition to option flag5, write the ID of the spt particles and structures
% into the video.
if flag6
    ind_spt = find(any(s,2)'); % SPT particle ID of interest
    Nspt = length(ind_spt);
    ind_str = find(any(s,1)); % structure ID of interest
    Nstr = length(ind_str);
    CC = C;
    CC.PixelIdxList = C.PixelIdxList(ind_str);
    CC.NumObjects = Nstr;
    L = labelmatrix(CC);
    cmap = [zeros(1,3,'uint8'); uint8(rand(Nstr,3)*255.5); repmat(uint8(255),1,3)];
    parfor k=1:size(L,3)
        for n=ind_spt % for each spt particle
            idx = find(q{n}(:,1)==k,1,'first'); % SPT particles present at time k
            if isempty(idx)
                continue
            end
            L(:,:,k) = addspt(L(:,:,k),q{n}(idx,2:3),10,Nstr+1);
            ix = q{n}(idx,2);
            iy = q{n}(idx,3);
            if iy>30
                iy = iy-30;
            else
                iy = iy+30;
            end
            L(:,:,k) = addnum(L(:,:,k),[ix iy],n,20,Nstr+1,suuji);
        end
        for n=1:Nstr % for each structure
            temp = L(:,:,k)==n;
            if ~any(temp(:)) % If there is no structure at that time, go to next structure
                continue
            end
            % write structure ID
            ix = find(any(temp,1),1,'first');
            iy = find(any(temp,2),1,'first');
            if iy>20
                iy = iy-20;
            else
                iy = find(any(temp,2),1,'last');
            end
            L(:,:,k) = addnum(L(:,:,k),[ix iy],ind_str(n),20,Nstr+1,suuji);
        end
    end
    v = VideoWriter([resultspathname '\wintld3.avi']);
    v.Quality = 100;
    v.FrameRate = FR;
    open(v);
    for k=1:ceil(T/100)
        img = ind2rgb_y(L(:,:,(k-1)*100+1:min(k*100,T)),cmap);
        writeVideo(v,img);
    end
    close(v)
end

% Co-localization analysis of random spt and PALM images generate random spt
% coordinates for 1-668 frame
qR = cell(RndN,1);
parfor n=1:RndN
    qR{n} = zeros(T,3);
    Rndt = 1:T;
    qR{n} = Rndt.';
    Rndx = 0+(Nx-0).*rand(T,1);
    qR{n} = cat(2,qR{n},Rndx);
    Rndy = 0+(Ny-0).*rand(T,1);
    qR{n} = cat(2,qR{n},Rndy);
end
   
% Examine the distance of each structure approaching a random spt.
clear rR
parfor k=1:N % for each structure
    for n=1:RndN % for each random spt
        rR{n,k} = wintld_beta_subfuncR(qR{n},PixelIdxList{k},Nx,Ny,T);
    end
end

% Make a frequency distribution table of the distance between the random spt and
% the nearby structure and the probability of the existence of the random spt
if flagR3
    trnd = [[],[]];
    for n=1:RndN % for each random spt(1-668 frame)
        trndk = [[],[]];
        for k=1:N % for each structure
            ind = ~isnan(rR{n,k}(:,2));
            trndn = [rR{n,k}(ind,1),rR{n,k}(ind,2)];
            trndk = [trndk; trndn];
        end
        trnd = [trnd; trndk];
    end
    outtrnd = Hist_table(trnd(:,2), BinSize, BinStart);
    outputtrndname = fullfile(resultspathname,'Rnd_ColoDistance_Hist.xls');
    writecell(outtrnd,outputtrndname)
    close all
end

% All structures and all random spt are made into a pseudo-color video and saved in an avi file
if flagR4
    L = labelmatrix(C);
    cmap = [zeros(1,3,'uint8'); uint8(rand(N,3)*255.5); repmat(uint8(255),1,3)];
    parfor k=1:size(L,3)
        for n=1:RndN % for each spt particle
            idx = find(qR{n}(:,1)==k,1,'first'); % SPT particles present at time k
            if isempty(idx)
                continue
            end
            L(:,:,k) = addspt(L(:,:,k),qR{n}(idx,2:3),10,N+1);
        end
    end
    v = VideoWriter([resultspathname '\wintld1_random.avi']);
    v.Quality = 100;
    v.FrameRate = FR;
    open(v);
    for k=1:ceil(T/100)
        img = ind2rgb_y(L(:,:,(k-1)*100+1:min(k*100,T)),cmap);
        writeVideo(v,img);
    end
    close(v)
end

%% subfunctions

%% Calculate the distance between the structure and the SPT particle
function r = wintld_beta_subfunc2(q,PixelIdxList,Nx,Ny,T)

conn = 8; % Difinition of connectivity
          % 8:8 neighborhoods including diagonally, 
          % 4:4 neighborhoods not including oblique direction
% Coordinate of SPT particle
spt_y = q(:,3); % y coordinate
spt_x = q(:,2); % x coordinate
spt_t = q(:,1); % t coordinate

% Structure
[str_y,str_x,str_t] = ind2sub([Ny Nx T],PixelIdxList);

% Structures and spt particles that do not intersect in time are not analyzed.
k3 = intersect(spt_t,str_t);
if isempty(k3)
    r = [];
    return
end

r = zeros(length(k3),4);
r(:,1) = k3(:);

for k=1:length(k3)
    ind_spt = find(spt_t==k3(k));
    ind_str = find(str_t==k3(k));
    d = (str_y(ind_str)-round(spt_y(ind_spt))).^2+(str_x(ind_str)-round(spt_x(ind_spt))).^2;
    [r(k,2),idx] = min(d);
    r(k,2) = sqrt(r(k,2));
    r(k,3:4) = [str_x(ind_str(idx)) str_y(ind_str(idx))];
    % If the SPT particle enters the interior of the structure,
    % the distance to the structure is -1 times the distance to the nearest boundary
    if r(k,2)==0
        bw2 = full(sparse(str_y(ind_str),str_x(ind_str),true,Ny,Nx));
        B = bwboundaries(bw2,conn);
        d = (B{1}(:,2)-spt_x(ind_spt)).^2+(B{1}(:,1)-spt_y(ind_spt)).^2;
        [r(k,2),idx] = min(d);
        r(k,2) = -sqrt(r(k,2));
        r(k,3:4) = [B{1}(idx,1) B{1}(idx,2)];
    end
end

%% Calculate the distance between the structure and the random SPT
function rR = wintld_beta_subfuncR(qR,PixelIdxList,Nx,Ny,T)

conn = 8; % Difinition of connectivity
          % 8:8 neighborhoods including diagonally, 
          % 4:4 neighborhoods not including oblique direction
% Coordinate of SPT particle
spt_y = qR(:,3); % y coordinate
spt_x = qR(:,2); % x coordinate
spt_t = qR(:,1); % t coordinate

% structure
[str_y,str_x,str_t] = ind2sub([Ny Nx T],PixelIdxList);

% Structures and spt particles that do not intersect in time are not analyzed.
k3 = intersect(spt_t,str_t);
rR = zeros(length(k3),2);
rR(:,1) = k3(:);

for k=1:length(k3)
    ind_spt = find(spt_t==k3(k));
    ind_str = find(str_t==k3(k));
    d = (str_y(ind_str)-round(spt_y(ind_spt))).^2+(str_x(ind_str)-round(spt_x(ind_spt))).^2;
    [rR(k,2),idx] = min(d);
    rR(k,2) = sqrt(rR(k,2));
    rR(k,3:4) = [str_x(ind_str(idx)) str_y(ind_str(idx))];  
    % If the SPT particle enters the interior of the structure,
    % the distance to the structure is -1 times the distance to the nearest boundary
    if rR(k,2)==0
        bw2 = full(sparse(str_y(ind_str),str_x(ind_str),true,Ny,Nx));
        B = bwboundaries(bw2,conn);
        d = (B{1}(:,2)-spt_x(ind_spt)).^2+(B{1}(:,1)-spt_y(ind_spt)).^2;
        [rR(k,2),idx] = min(d);
        rR(k,2) = -sqrt(rR(k,2));
        rR(k,3:4) = [B{1}(idx,1) B{1}(idx,2)];
     end
    
end

%% Create a frequency distribution table
function Out_table = Hist_table(Data, BinSize, BinStart)

% Calculate Bin
BinEnd = max(Data);
Bin_N = ceil((BinEnd - BinStart) / BinSize) +1;
Bin_array(1) = BinStart;

% Create bin array
for k=2:Bin_N
    Bin_array(k) = Bin_array(k-1) + BinSize;
end

% Calculate the median for each bin
Bin_center = Bin_array(1:end-1) + BinSize/2;

% draw histogram
figure(1);
Hist = histogram(Data, Bin_array);

% create table
Out(:,1) = transpose(Bin_center);
Out(:,2) = transpose(Hist.BinCounts);
N = sum(Out(:,2));
Out(:,3) = Out(:,2) / N * 100;
Out = num2cell(Out);
Header = {'Bin center', 'Raw count', '% Count'};
Out_table = vertcat(Header, Out);

% delete histogram
F = gcf;
% Close current figure
close(figure(F.Number))

%% Write a number v with size r at coordinate p of grayscale image img
function img = addnum(img,p,v,ny,kido,suuji)
% img : Image to be written
% p : coordinate (x,y)
% v : the number to be written
% ny : the height of the number
% kido : the intensity of the number in the image
% suuji : the image data of number 0-9

if nargin<6
    load('suuji','suuji')
end
if nargin<5
    kido = 255;
end

cv = num2str(v);
for k=1:length(cv)
    nv(k) = str2num(cv(k));
    if nv(k)==0
        nv(k) = 10;
    end
end
suujiimg = suuji{nv(1)};
for k=2:length(nv)
    suujiimg = cat(2,suujiimg,suuji{nv(k)});
end
suujiimg = imresize(double(suujiimg),ny/size(suujiimg,1))>0.5;
if isa(img,'uint8')
    suujiimg = uint8(suujiimg)*uint8(kido);
elseif isa(img,'uint16')
    suujiimg = uint16(suujiimg)*uint16(kido);
end
nx = size(suujiimg,2);

Nx = size(img,2);
Ny = size(img,1);
p = round(p);

x1 = p(1);
x2 = p(1)+nx-1;
if x2>Nx
    x1 = x1-(x2-Nx);
    x2 = Nx;
end
y1 = p(2);
y2 = p(2)+ny-1;
if y2>Ny
    y1 = y1-(y2-Ny);
    y2 = Ny;
end

if ismatrix(img)
    img(y1:y2,x1:x2) = max(img(y1:y2,x1:x2),suujiimg);
elseif ndims(img)==3 && size(img,3)==3
    img(y1:y2,x1:x2,1) = max(img(y1:y2,x1:x2,1),suujiimg);
    img(y1:y2,x1:x2,2) = max(img(y1:y2,x1:x2,2),suujiimg);
    img(y1:y2,x1:x2,3) = max(img(y1:y2,x1:x2,3),suujiimg);
end

%% Write '+' of size n in the coordinates p of the grayscale image img with the intensity kido.
function img = addspt(img,p,n,kido)

Nx = size(img,2);
Ny = size(img,1);
p = round(p);

img(max(1,p(2)-1):min(Ny,p(2)+1),max(1,p(1)-n):min(Nx,p(1)+n)) = kido;
img(max(1,p(2)-n):min(Ny,p(2)+n),max(1,p(1)-1):min(Nx,p(1)+2)) = kido;


%% Matlab's original ind2rgb function can be applied to videos at high speed
function img2 = ind2rgb_y(img1,cmap)
Ny = size(img1,1);
Nx = size(img1,2);
Nt = size(img1,3);
img1 = img1(:)+1;
img_r = reshape(cmap(img1,1),[Ny Nx 1 Nt]);
img_g = reshape(cmap(img1,2),[Ny Nx 1 Nt]);
img_b = reshape(cmap(img1,3),[Ny Nx 1 Nt]);
img2 = cat(3,img_r,img_g,img_b);
