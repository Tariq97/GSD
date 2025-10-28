clear
close all
clc


L = 20*1000;   % 2D field geometry in m
W = 10*1000; 

samp = [200;200];   % in m

lx = [0:samp(2):L];
lz = [0:samp(2):W];

% correlation length parameters for Von Karman/Gaussian
ax = 4; az = 4; H = 0.77; 
corr1 = [az ax H]; lmin =  10*25; lmax = L; 

realz = 458; % seed value

acf = 'ak';  % 'gs' for Gaussian and 'ak' for Von-Karman

field0 = SpecSyn3([W L],samp,[az*1000 ax*1000 H],acf,realz,lmin,lmax);

sigma = std(field0(:));
field0_old = field0;

figure;
imagesc(lx./1000,lz./1000,field0_old);
% colormap(slipcolor); colorbar; 
axis tight equal;
set(gcf,'position',[200,200,500,500])

% Taking the first line of this field (field0) as the reference trace for generating different
% fields....!!

reference_trace = field0(1,:);

%% 

num_realz = 10;  % Number of fields to generate

A_new = [];


A_new = zeros(size(field0));
A_new(1,:) = reference_trace;


% GSD parameters
ssd.neighbor_radius = 2; % in km; 
ssd.corr = corr1;   % correlation lengths of new fields.!!


% Takes around ~5min for Von karman example.
[A_new2,~,~] = GSD(A_new,ssd,[W,L],num_realz,samp./1000,acf); 
   
%% Plot fields


for i=1:num_realz
    figure(2)
    subplot(5,2,i)
    imagesc(lx./1000,lz./1000,A_new2{i,1});
    axis tight equal;
    set(gcf,'position',[200,200,500,900])




end