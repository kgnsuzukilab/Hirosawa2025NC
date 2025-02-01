function [img2,spot,tse] = gen_palm_zeta(spot,Nx,Ny,Nframe,A,D,W)
% This function is the calculation engine for PALM video image generation
% used in gen_palm_eta and DoCanalysis
% gen_palm_zeta reconstructs the temporal changes in the structure from the spot data
% output by ThunderSTORM.
% spot : Matrix consisting of [x y t s]
% Nx, Ny : Horizontal and vertical size of the original video and PALM video
% Nframe : Frame size of the original video
% A : Parameter to determine kernel width, default 8
% D : Shift width for making PALM videos
% W : Window width for making PALM videos
% img2 : Generated PALM video image
% spot : Removed spots that were out of range of Nx, Ny, Nframe
% tse : The range of times that make up each frame of a PALM video
% est_th1.m or est_th2.m is necessary to determine the threshold for binarization
% 2021.9,11 Coded by Y. Yokota

x1 = round(spot(:,1));
y1 = round(spot(:,2));
t1 = spot(:,3);
s = spot(:,4);

if Nframe<max(t1)
    disp('The number of frames in the original video is smaller than the maximum number of frames in spot.')
    disp('Resets the number of frames in the original video to the maximum number of frames in spot.')
    Nframe = max(t1);
end

% uncertainty
u = cumsum(histcounts(s,0:0.01:10))/length(s);
u = find(u>0.99,1,'first')*0.01;
u = linspace(0,u,21);
u(1) = [];
s = uint8(max(1,min(20,round(s/u(1)))));

% Pre-create a spatial kernel for each value of uncertainty
for k=1:length(u)
    b(k).size = ceil(3*u(k)*A);
    b(k).kernel = single(fspecial('gaussian',b(k).size*2+1,double(u(k)*A)));
end

if isinf(W) || isinf(D)
    M = 1;
    Nt = 1;
    T = Nframe;
    W = T;
    D = T;
    tse = [1 T];
else
    M = fix(Nframe/D); % The number of frames for intermediate PALM video
    Nt = M-W/D+1; % The number of frames of PALM video.
                  % W/D represents the number of intermediate frames par a PALM video
    tse = (0:Nt-1)'*D+[1 W]; % Upper and lower limits of time that make up each image in the PALM video
    T = tse(end,2); % Maximum time value of the target spot
end

% Remove spots whose coordinates are beyond the range of Nx, Ny, and T,
% and spots whose uncertainty is too small.
ind = x1>Nx | y1>Ny | t1>T | s<2;
spot(ind,:) = [];
x1(ind) = [];
y1(ind) = [];
t1(ind) = [];
s(ind) = [];

% Create a PALM video from spot
img2 = zeros(Ny,Nx,Nt,'single');
img1 = zeros(Ny,Nx,W/D,'single');
for n=1:M
    ind = find(t1'>(n-1)*D & t1'<=n*D);
    if isempty(ind)
        continue
    end
    img0 = zeros(Ny,Nx,'single');
    for k = ind
        j1 = y1(k)-b(s(k)).size;
        j3 = max(1,2-j1);
        j1 = max(1,j1);
        j2 = y1(k)+b(s(k)).size;
        j4 = max(0,j2-Ny);
        j2 = min(Ny,j2);
        i1 = x1(k)-b(s(k)).size;
        i3 = max(1,2-i1);
        i1 = max(1,i1);
        i2 = x1(k)+b(s(k)).size;
        i4 = max(0,i2-Nx);
        i2 = min(Nx,i2);
        img0(j1:j2,i1:i2) = img0(j1:j2,i1:i2) + b(s(k)).kernel(j3:end-j4,i3:end-i4);
    end
    img0 = img0/D;
    if n<=W/D
        img1(:,:,n) = img0;
    end
    if n==W/D
        img2(:,:,1) = sum(img1,3);
    end
    if n>W/D
        img2(:,:,n-W/D+1) = img2(:,:,n-W/D)+img0-img1(:,:,mod(n-1,W/D)+1);
        img1(:,:,mod(n-1,W/D)+1) = img0;
    end
end
clear img1 img0
img2 = img2/(W/D); % averaging


ind = round(linspace(1,Nt,min(Nt,10)));
img3 = img2(:,:,ind);
Th = est_th2(img3); % Automatic determination of threshold value for binarization 
img2 = img2/Th; % Standardize PALM video with determined thresholds.


%% Automatic determination of threshold value for binarizing image
function Th = est_th2(x,K)
% x : Data, image, or video image
% K : Resolucion. The threshold value is determined with a resolution of K divisions between the largest and smallest values of x
% 2021.8.19 Coded by Y.Yokota

% default setting
if nargin<2
    K = 1000;
end

x = double(x(:));

Fun = @(Th)fun(Th,x);
options = optimset('Display','off','TolX',(max(x)-min(x))/K);

Th = fminbnd(Fun,min(x),max(x),options);


%% evaluation function for est_th2
function fval = fun(Th,x)
ind = x>Th;
fval = sqrt((var(x(ind))*sum(ind)+var(x(~ind))*sum(~ind))/length(x));
% To use Otsu's binarization method, comment out the following
% fval = fval/sqrt(((mean(x(ind)-mean(x)).^2*sum(ind)+(mean(x(~ind))-mean(x)).^2*sum(~ind))/(length(x))));
