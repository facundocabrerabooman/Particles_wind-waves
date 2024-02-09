close all
clear
clc

% input
fpath0 = '/Volumes/landau2/Particle_wind-waves/base_case/';
cd(fpath0)
%Istart = 2; % better to ignore the 1st image
 
%%
fpath = [fpath0 filesep 'raw'];
% get the folder contents
addpath(genpath('/Users/fcb/Documents/GitHub/Particles_wind-waves'));
%dfolders = FunSubfolder(fpath);
dfolders = [];

%% preproc & split


    image_list = dir([fpath filesep 'Camera1' filesep '*.tiff']);
    img_num = size(image_list,1);
    % 
    % Istart = 2;
    % Iend = img_num-100;
    Istart = 2250;
    Iend = 2500;

    preproc_dirp = [fpath0 filesep 'preproc' filesep 'preproc_particle'];
    mkdir(preproc_dirp)
    
    %%% get Backgrounds
    %cam 1
    bkg1_originalSize = getBkg(fpath,'Camera1',Istart,Iend,5,[]);
    %cam2
    bkg2_originalSize = getBkg(fpath,'Camera2',Istart,Iend,5,[]);

    counter = 0;
   
    for k=Istart:Iend
        counter = counter+1;
        k/Iend

        %%% cam1
        fname=[fpath filesep 'Camera1' filesep 'frame_' num2str(k,'%06d') '.tiff'];
        Im1_originalSize = imread(fname);

        intensity_thr = 5e2;
        [Im01,Im1t,Im1p]=gray_preprocessing(Im1_originalSize,intensity_thr,7,2,bkg1_originalSize);

        fnameo = ['cam1_frame_preproc_' num2str(counter,'%06d') '.tiff'];
        imwrite(uint16(Im1p),[preproc_dirp filesep fnameo]);


 
        %%% cam2
        fname=[fpath filesep 'Camera2' filesep 'frame_'  num2str(k,'%06d') '.tiff'];
        Im2_originalSize = imread(fname);

        intensity_thr = 1e3;
        [Im02,~,Im2p]=gray_preprocessing(Im2_originalSize,intensity_thr,6,1,bkg2_originalSize);

        fnameo = ['cam2_frame_preproc_' num2str(counter,'%06d') '.tiff'];
        imwrite(uint16(Im2p),[preproc_dirp filesep fnameo]);

    end

stop
%% Test camera 1

if 1==pi
    figure;
    subplot(2,1,1)
    imshow(Im1_originalSize);axis equal;
    [C, ~] = imfindcircles(imbinarize(Im1p),[1 15]);
    viscircles(C,ones(size(C,1),1).*12,'Color','r')


    subplot(2,1,2)
    imagesc(Im1p);axis equal; hold on
    viscircles(C,ones(size(C,1),1).*12,'Color','r')

    figure;imagesc(Im1t);axis equal;hold on
    %viscircles(C,ones(size(C,1),1).*5,'Color','r')

end
%% Test camera 2
if 1==pi
    figure;
    subplot(2,1,1)
    imshow(Im2_originalSize);axis equal;
    [C, ~] = imfindcircles(imbinarize(Im2p),[1 15]);
    viscircles(C,ones(size(C,1),1).*12,'Color','r');


    subplot(2,1,2)
    imagesc(Im2p);axis equal; hold on
    viscircles(C,ones(size(C,1),1).*12,'Color','r');

    figure;imagesc(Im2t);axis equal;hold on

end
