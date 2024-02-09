%% Blackout image

cd('/Volumes/landau2/Particle_wind-waves/calibration_2_1_24/')

im1 = imread('/Volumes/landau2/Particle_wind-waves/calibration_2_1_24/cam1.tiff');
im1_black = uint16(zeros(1024,1280));
im1_black(1:700,540:end) = im1(1:700,540:end);

imwrite(uint16(im1_black),['Camera1_' num2str(1,'%05d') '.tif'])
clf;imshow(im1_black)


im2 = imread('/Volumes/landau2/Particle_wind-waves/calibration_2_1_24/cam2.tiff');
im2_black = uint16(zeros(1024,1280));
im2_black(300:750,600:1020) = im2(300:750,600:1020);

imwrite(uint16(im2_black),['Camera2_' num2str(1,'%05d') '.tif'])
clf;imshow(im2_black)



%% Average and filter calibration images
clfc

cam = 3;
for frame=1:7

%%%%%%%%%
pathcalib = ['/Users/fcb/AuxFiles/calib_12-13-23/Camera' num2str(cam) filesep num2str(frame)];

pathout = '/Users/fcb/AuxFiles/calib_12-13-23/';

fname = ['cam' num2str(cam) '_frame_preproc_00000' num2str(frame)];

mkdir(pathout)
cd(pathout)
image_list = dir([pathcalib filesep '*.tiff']);


imsum = zeros(512,1280);

 im = imread([pathcalib filesep image_list(5).name]);
 %imshow(im)
% h = imfreehand;

for i=5:numel(image_list)
    
    im = imread([pathcalib filesep image_list(i).name]);
    
    im=imadjust(im);

    imsum = imadd(double(im),imsum);
    
end

imaver = uint16(imsum./numel(image_list));
%imaver(imaver>2e4)=0;

st=strel('disk',1);
imaver=imopen(imaver,st);
imaver=imadjust(imaver);

imaver = imaver.*(5.3e4/65500);

imshow(imaver)

%cam1 = 65526
%cam2 = 53000
%cam3 = 53000
imwrite(uint16(imaver),[pathout filesep fname '.tiff'])
end
stop
%% Filter averaged frames

pathcalib = ['/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop ' ...
    'Tower Multiphase Flow Project/Data/bigbiginertial/calibs/calib_27nov_b/Camera1/'];

pathout = ['/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop ' ...
    'Tower Multiphase Flow Project/Data/bigbiginertial/calibs/calib_27nov_b/mean_filtered/'];
mkdir(pathout)

fname = 'cam1_frame_preproc_000001';


im = imread([pathcalib filesep fname '.tiff']);
im(im<1500)=0;

st1 = strel('disk', 10);
st2 = strel('disk', 5);

imo=imopen(im,st2);
imo=imclose(imo,st1);



%% Cut calibration frames

pathcalib = '/Users/FC/Aux_files/Calib_26oct23/';
pathout = '/Users/FC/Aux_files/Calib_26oct23/cut/';
mkdir(pathout)
image_list = dir([pathcalib filesep '*.tiff']);

for i=1:numel(image_list)
    
    im = imread([pathcalib filesep image_list(i).name]);
    im = im(285:796,321:1600); 
    imwrite(uint16(im),[pathout filesep image_list(i).name])


end

