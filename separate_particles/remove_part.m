function [Imt,Im] = remove_part(Im,PartRadius)

% st = strel('disk',4);
% stdil = strel('disk',round(PartRadius/3));
% 
% Imo = imopen(Im,st);
% %Imd = imdilate(imbinarize(imopen(Im,st)),stdil);
% Imd = imdilate(Imo,stdil);
% 
% %% refine particle image
% Rp=sqrt(numel(find(Imd~=0)))/pi;
% %Rp=PartRadius;
% 
% Rmin = max(1*Rp,4);
% Rmax = max(1.6*Rp,15);
% [C, R] = imfindcircles(Imd, cast([Rmin Rmax],class(Im)));
% 
% if size(C,1) == 1
%     roi=images.roi.Circle('Center',C,'Radius',R);
%     mask = createMask(roi,size(Imd,1),size(Imd,2));
% else
%     mask = Imd;
% end
% Imp = immultiply(Im,cast(mask,class(Im)));
% Imt = immultiply(Im,cast(abs(1-mask),class(Im)));

st = strel('disk',3);
Im = imopen(Im,st);
Im(Im<1e4)=0;
Im = imgaussfilt(Im,1);

Imt=Im;














