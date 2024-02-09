function [Imt,Imp] = remove_part(Im,PartRadius)

st = strel('disk',4);
stdil = strel('disk',round(PartRadius/3));

Imo = imopen(Im,st);
%Imd = imdilate(imbinarize(imopen(Im,st)),stdil);
Imd = imdilate(Imo,stdil);

%% refine particle image
Rp=sqrt(numel(find(Imd~=0)))/pi;
%Rp=PartRadius;

Rmin = max(1*Rp,4);
Rmax = max(1.6*Rp,15);
[C, R] = imfindcircles(Imd, cast([Rmin Rmax],class(Im)));

if size(C,1) == 1
    roi=images.roi.Circle('Center',C,'Radius',R);
    mask = createMask(roi,size(Imd,1),size(Imd,2));
else
    mask = Imd;
end
Imp = immultiply(Im,cast(mask,class(Im)));
Imt = immultiply(Im,cast(abs(1-mask),class(Im)));

%Imp = imgaussfilt(Imp,1);
st = strel('disk',4);
Imt = imopen(Imt,st);
Imt = imgaussfilt(Imt,1);

st = strel('disk',8);
Imt = imopen(Imt,st);
Imp(Imp<1e4)=0;
Imp = imgaussfilt(Imp,1);














