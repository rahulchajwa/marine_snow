%for file format see https://doi.org/10.5061/dryad.v15dv4253
%this code takes as an input the output of the PIV analysis in the repository on PIV-gravity-machine-on-oceans
currentDir = pwd;
addpath('/home/rahul/Documents/Maine analysis/analysis codes/Colormaps/Colormaps (5)/Colormaps');
addpath('/home/rahul/Documents/Maine analysis/analysis codes/smoothn.m');
Dd = dir; %list of subfolders
l =1;
for i=3:length(Dd)
%for ii=6:length(D)
cd(Dd(i).name);
%cd('DF1');
cd('registered_output')
a = dir('*.tif');
nn = numel(a);
lis = 1:1:nn;
Flist = dir('*.mat'); %list only data files
n = length(Flist);
PIV_flname = sprintf(Flist(1).name);
load(PIV_flname)

%frame_n1 = sscanf(PIV_flname(13),'%d');
%Image_flname = sprintf(a(frame_n1).name);
%frame_n2 = sscanf(Image_flname,'%d .tif ');

%frame_n1 = sscanf(PIV_flname(13),'%d');
Image_flname = sprintf(a(1).name);
frame_n2 = sscanf(Image_flname,'%d .tif ');

I = imread(Image_flname); %reading the same tif image as
clearvars Flist
cd .. 
Flist = dir('track*.csv'); %list only data files
n=length(Flist);
X1=readtable(Flist(n).name);
cd('registered_output')
tt = X1.DF1(:,1);
index = find(~cellfun(@isempty,tt(:,1)'));
triggerTime = X1.Time(index,1);
delT = triggerTime(frame_n2 +1+1) - triggerTime(frame_n2+1); %adding +1 for the frames starts with zero
cd ..
%% colobar of vertical velocity
[r, c] = size(v_filt{1,1});
x = 0:1:r-1;
y = 0:1:c-1;
[X,Y] = meshgrid(y,x);
imshow(I)
clearvars U V
U = u_filt{1,1} *0.828*0.1*0.864/delT;   %in meters/day
V = v_filt{1,1} *0.828*0.1*0.864/delT;   %in meters/day
%U(U>=0) = 0;
V(V>=0) = 0;
imagesc(y,x,V); hold on
li = streamslice(X,Y,u_filt{1,1}*0.828*0.1*0.864/delT,v_filt{1,1}*0.828*0.1*0.864/delT); hold on
colormap(viridis)
%clim([0 60*0.828*0.1*0.864/delT]);
set(li,'LineWidth',3)
set(li,'Color','w');
h = colorbar;
set(h, 'ylim', [-50 5])
%set(h, 'ylim', [-50 0])
axis off
colorbar off;
daspect([1 1 1])
%set(gca, 'units', 'normalized'); %Just making sure it's normalized
%Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                 %[Left Bottom Right Top] spacing
%NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
%NewPos = [0    0  1   1];
%set(gca, 'Position', NewPos);

set(gca, 'color', 'none');
exportgraphics(gca,'flow_field.png','BackgroundColor','none');

%% masking the particle on it
BW = imbinarize(I);
%imshow(BW)
BW2 = imfill(BW,'holes');
%imshow(BW2)
I1 = medfilt2(BW2);
I5 = bwareaopen(I1, 500);
sz = sqrt((bwarea(I5)*0.828*0.828*10^-8)/pi); %equivalent circle in cm^2
%F = getframe(gcf);
%[A, Map] = frame2im(F);
A = imread('flow_field.png');
close
%imshow(A)
I2 = imcomplement(I1);
%imshow(I2)
%imshow(imerode(I2))
%imshow(I2)
%imshow(Z)
Z= immultiply(I1,I);
[m,n] = size(BW2);
A2 = imresize (A,[m,n]);
Z1(:,:,1)= immultiply(I2,A2(:,:,1));
Z1(:,:,2)= immultiply(I2,A2(:,:,2));
Z1(:,:,3)= immultiply(I2,A2(:,:,3));
%imshow(Z1)
%imshow(Z)
Z2(:,:,1) = imadd(Z1(:,:,1),Z);
Z2(:,:,2) = imadd(Z1(:,:,2),Z);
Z2(:,:,3) = imadd(Z1(:,:,3),Z);
imshow(Z2)
%M = v_component;
M = v_filt{1,1};
M(M<-5) = 0;
M(~M == 0) = 1;
Mbw = imresize (M,[m,n]);
%I = imread('0002510.tif');
%BW = imbinarize(I);
%imshow(BW)
%BW2 = imfill(BW,'holes');
s = regionprops(Mbw,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
if(isempty(s)==0)
% Calculate the ellipse line
theta = linspace(0,2*pi);
col = (s.MajorAxisLength/2)*cos(theta);
row = (s.MinorAxisLength/2)*sin(theta);
M = makehgtform('translate',[s.Centroid, 0],'zrotate',deg2rad(-1*s.Orientation));
Di = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];
% Visualize the result
%figure
%imshow(Mbw)
hold on
plot(Di(1,:),Di(2,:),'--','Color','r','LineWidth',3)
msz = (s.MajorAxisLength/2)*0.828*10^-4;

set(gca, 'color', 'none');
exportgraphics(gca,'mucus_halo.png','BackgroundColor','none')
%myfilename = sprintf('mucus_halo.png');
%print(gcf,myfilename,'-dpng','-r200');

image = imread('mucus_halo.png');
%figure
clearvars n lis r c
if (frame_n2 == 0)
    lwindow = 0;
else
lwindow = 20;
end
rwindow = 20; %chose the width to be 40
%%NOTE -- change the index of the index in the fitting function
%plot(X1.Time(:,1)-X1.Time(1,1),(X1.Z_objStage(:,1)-X1.Z_objStage(1,1))*0.1,'o'); hold on
f=fit(X1.Time(:,1)-X1.Time(1,1),(X1.Z_objStage(:,1)-X1.Z_objStage(1,1))*0.1,'poly1');
%f=fit(X1.Time(index(frame_n2+1)-lwindow:index(frame_n2+1)+rwindow,1)-X1.Time(index(frame_n2+1)-lwindow,1),(X1.Z_objStage(index(frame_n2+1)-lwindow:index(frame_n2+1)+rwindow,1)-X1.Z_objStage(index(frame_n2+1) - lwindow,1))*0.1,'poly1');
%v_temp = f.p1 *864; %fitted velocity in m/day
v_temp = f.p1;
XX1(1,1) = v_temp; XX1(1,2) = sz; XX1(1,3) = (s.MajorAxisLength/2)*0.828*10^-4;
XX1(1,4) = (s.MinorAxisLength/2)*0.828*10^-4; XX1(1,5) = s.Orientation ; %array for concatenate

prompt = "Is this data clean? yes = 1, No = 0 ";
Cldata = input(prompt);
close
%cd ..
%cd ..
%Flist = dir('track*.csv'); %list only data files
%n=length(Flist);
%m=0;
%for jj=1:n
%m=m+1;
%X=readtable(Flist(n).name);
%plot(X.Time(:,1)-X.Time(1,1),X.Z_objStage(:,1)-X.Z_objStage(1,1),'o'); hold on
%f=fit(X.Time(:,1)-X.Time(1,1),X.Z_objStage(:,1)-X.Z_objStage(1,1),'poly1');
%v_temp(m) = f.p1;
%end
%v(l) = mean(v_temp);
%clearvars X
%clearvars n lis
Stru.stokes{l} = XX1;
Stru.flow{l} = image;
Stru.clean{l} = Cldata;
Stru.imraw{l} = I;
clearvars -except cm_magma cm_inferno cm_plasma cm_viridis Stru l Dd i
cd ..
l = l+1;
else
    clearvars -except cm_magma cm_inferno cm_plasma cm_viridis Stru l Dd i
    close
    cd ..
end
end
