%% Main test code for the paper 'Semi-calibrated Near Field Photometric Stereo'. from CVPR 2017
% Author Fotios Logothetis fl302@cam.ac.uk
clear all 
close all

results_dir='./results/';
data_dir='./data/';

name='scene'; 
% name='hand'; 

data_f=[data_dir,name,'.mat'];  

%% Load data
%Data MUST contain: 
%--I: nrows x ncols x nimages . image raw data. values assumed to be
% NORMALISED in [0,1]
%--mask: nrows x ncols. forground mask
%--S: 3 x nimages vector containing light source position in pixels. It is
%VERY IMPORTANT to get this right otherwise the reconstruction will fail.
%run the test_axis script to see a demonstration of the convention assumed
%--f: camera focal length in pixels
%--cc: [x0,y0]. camera principal point
%--mean_distance: mean distance between camera plane and object (in
%pixels). Just a rough estimate (e.g. measure with a ruler) should be ok
%--mu: nimages x 1 vector containing light source radial attenuation
%params. If not sure just set the zeros(nimages,1). It is not very important
%--Sd: 3 x nimages vector containing light source maximum illumination
%directions. If mu=0 it does not matter (as long as it non-zero). if not
%sure use [0;0;1] for each source
%--mm_to_px:absolute world scale (optional)

load(data_f) 
%% SOME TESTS
[nrows,ncols,nimages] = size(I);
assert(size(S,1)==3);
assert(size(S,2)==nimages);
assert(size(Sd,1)==3);
assert(size(Sd,2)==nimages);
assert(size(mu,1)==nimages);

if ~exist('mm_to_px','var')
    mm_to_px=1;
end
%% RESIZE
%This helps reducing computation time and RAM
%The L1 optimiser will need 10+GB of RAM if running at full resolution
resize_ratio =4;
%% Get chanel colors for rho/proper resizing. This assumes BGGR BAYER filter
Ir=I(2:2:end,2:2:end,:);
Ig1=I(1:2:end,2:2:end,:);
Ig2=I(2:2:end,1:2:end,:);
Ib=I(1:2:end,1:2:end,:);
%do chanels seperatelly and re-merge
Ir=imresize(Ir,[nrows,ncols]/(2*resize_ratio));
Ig1=imresize(Ig1,[nrows,ncols]/(2*resize_ratio));
Ig2=imresize(Ig2,[nrows,ncols]/(2*resize_ratio));
Ib=imresize(Ib,[nrows,ncols]/(2*resize_ratio));
%re-merge
I=zeros(nrows/resize_ratio, ncols/resize_ratio, nimages);

I(2:2:end,2:2:end,:)=Ir;
I(1:2:end,2:2:end,:)=Ig1;
I(2:2:end,1:2:end,:)=Ig2;
I(1:2:end,1:2:end,:)=Ib;

I = max(0.01,I);
   
mask=imresize(mask,[nrows,ncols]/resize_ratio);

S=S/resize_ratio;

f = f/resize_ratio;
cc = cc/resize_ratio;
mean_distance =mean_distance/resize_ratio;
mm_to_px =mm_to_px/resize_ratio;

[nrows,ncols,nb_images] = size(I);
%% group vars
cam.f=f;
cam.cc=cc;

S_struct.S=S;
S_struct.Sd=Sd;
% S_struct.Phi=Phi;
S_struct.mu=mu;
%% Misk opts
C =1*ones(nrows,ncols);  %initialise C as being Lambertian.
% 
opts.refine_C=1;
opts.shadow_threshold = 0.03; 
opts.saturation_thress=0.99;
%this is the semi-calibrated part. if 0, the code performs as standard near
%field PS (Mecca et al., SIAM 2016). if 1, light source brightness are
%calculated (Phi). If 2, light attenuation is also calculated
opts.calibrate=2; 
opts.estimate_shadows=0; %ray-trace shadow maps. expensive and not really essential
opts.parallelise=1; %use parfor in some calculations (e.g. patch attenuations) that can be parallelised

[X,Y,Z,N, rho,C_refined,phi_new] = semi_calibrated_ps(I, mask, mean_distance, cam, C,S_struct,opts);
mask_out=mask;
mask_out(isnan(Z))=0;
%% rho colors
%make almost black pixels black and dont white balance them
imean=mean(I,3);
rho(imean<opts.shadow_threshold)=0;
rho=0.5*rho;%rho has mean value 1, just scale to have a descent color range
%white balancing params for our camera
rw=1.02;
gw=1.1;
bw=0.88;
%again assume BGGR BAYER filter
rho_r=rho(2:2:end,2:2:end)*rw;
rho_g1=rho(1:2:end,2:2:end)*gw;
rho_g2=rho(2:2:end,1:2:end)*gw;
rho_b=rho(1:2:end,1:2:end)*bw;

rho_rgb=cat(3, rho_r, 0.5*(rho_g1+rho_g2),rho_b);

imwrite(rho_rgb, [results_dir,'rho_',name,'.jpg']);
imwrite((N+1)/2, [results_dir,'norm_',name,'.jpg']);
%%
XYZ = cat(3,X,Y,Z)/mm_to_px;
export_ply(XYZ,mask_out,[results_dir,name,'.ply'],rho_rgb); 