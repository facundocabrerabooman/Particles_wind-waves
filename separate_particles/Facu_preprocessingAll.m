close all
clear
clc

% input
fpath0 = '/Users/fcb/AuxFiles/Particles_wind&waves/dec14/';
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
    
    Istart=3;
    Iend = 7000;

    subfolders = fpath;
    preproc_dirp = [fpath0 filesep 'preproc' filesep 'preproc_particle'];
    mkdir(preproc_dirp)
    
    %%% get Backgrounds
    %cam 1
    bkg1_originalSize = getBkg(subfolders,'Camera1',Istart,Iend,100,[]);
    %cam2
    bkg2_originalSize = getBkg(subfolders,'Camera2',Istart,Iend,100,[]);

    counter = 0;
    for k=5000:7000%Istart:Iend
    %for k=Istart:Iend
        counter = counter+1;
        k/Iend

        %%% cam1
        fname=[subfolders filesep 'Camera1' filesep 'frame_' num2str(k,'%06d') '.tiff'];
        Im1_originalSize = imread(fname);
        intensity_thr = 1e4;
        [Im01,Im1t,~]=Facu_preprocessing(Im1_originalSize,intensity_thr,7,1,bkg1_originalSize);
        fnameo = ['cam1_frame_preproc_' num2str(counter,'%06d') '.tiff'];
        %imwrite(uint16(Im1p),[preproc_dirt filesep fnameo];
        imwrite(uint16(Im1t),[preproc_dirp filesep fnameo]);



        %%% cam2
        fname=[subfolders filesep 'Camera2' filesep 'frame_'  num2str(k,'%06d') '.tiff'];
        Im2_originalSize = imread(fname);
        %         Im2 = cast(zeros(1080,1920),class(Im2_originalSize));
        %         Im2(285:796,321:1600) = Im2_originalSize;
        intensity_thr = 1e4;
        [Im02,~,Im2p]=Facu_preprocessing(Im2_originalSize,intensity_thr,10,1,bkg2_originalSize);
        fnameo = ['cam2_frame_preproc_' num2str(counter,'%06d') '.tiff'];
        imwrite(uint16(Im2p),[preproc_dirp filesep fnameo]);
%                figure(10);
%                  subplot(2,1,1);imagesc(Im2_originalSize);axis equal
%                  subplot(2,1,2);imagesc(Im2p);axis equal
%                  pause(0.05)


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
