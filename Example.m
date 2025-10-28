 clear
close all
clc


L_min = 20*1000;   % Reference 2D field 
W_min = 14*1000; 

L = [ 25 30 35 40 ].*1000;   % Increased 2D fields for which VK fields are generated with reference field taken from above.
W = [ 14 14 14 14 ].*1000;
 

% -1 standard deviation.!!

samp = [200;200];

lx = [0:samp(2):L_min];
lz = [0:samp(2):W_min];


ax = 4; az = 4; H = 0.77; corr1 = [az ax H]; lmin =  10*25; lmax = L_min; 
realz = 458
acf = 'ak';  % 'gs' for Gaussian and 'ak' for Von-Karman


field0 = SpecSyn3([W_min L_min],samp,[az*1000 ax*1000 H],acf,realz,lmin,lmax);

sigma = std(field0(:));
field0_old = field0;

figure;
imagesc(lx./1000,lz./1000,field0_old);
% colormap(slipcolor); colorbar; 
axis tight equal;
set(gcf,'position',[200,200,500,500])



%% Generate pseudo-fields around the main field for different combinations of (L,W)
close all;

% first, generate entirely pseudo fields for increased L,W
[nz,nx] = size(field0);

for i=1:length(L)
    
   L1{i,1} = [0:samp(2):(L(i))]'./1000; 
   W1{i,1} = [0:samp(1):(W_min)]'./1000; 

   
   M = [ length(L1{i,1}), length(L1{i,1}) ];  % For now making both same
  
   
   field0 = field0_old;

   % Adding columns of zeros left/right
   if rem(M(2)-nx,2)==1
       field0 =  [ zeros([size(field0,1),(M(2)-nx-1)/2]) field0 zeros([size(field0,1),(M(2)-nx+1)/2]) ];
   else
       field0 =  [ zeros([size(field0,1),(M(2)-nx)/2]) field0 zeros([size(field0,1),(M(2)-nx)/2]) ];
   end
       
   
    A_new{i,1} = field0; % all same fields; new fields will be created.!!
      
    
end


ssd.neighbor_radius = 2; % in km; 
ssd.corr = corr1; 

n_stoch=5;  % number of realizations per each L-W combination...

% Takes around ~3min for Von karman,.
parfor i=1:length(L)
    
    [A_new2{i,1},A_stoch{i,1},realz_stoch{i,1}] = GSD(A_new{i,1},ssd,[W(i),L(i)],n_stoch,samp./1000,acf); 
   
end


for i=1:length(L) 
    for j=1:n_stoch
         slip_lw{i,j} = A_new2{i,1}{j,1}; 
    end
end




%% Plots new fields


close all
load('slipcolor.mat')

xx = 0.1; 
zz = 0.1*(W_min/L_min);

for pl=1 %1:10
    figure(pl)
    [ha, pos] = tight_subplot(2,3,[0 0],[0.1 0],[0, 0]);
    
    for i=1:length(ha) 
       ha(i).Position(3) = 0; ha(i).Position(4) = 0;
    end

% first row of plots.!!

axes(ha(1))
imagesc([0:0.2:L_min./1000],[0:0.2:W_min./1000],field0_old); colormap(slipcolor); colorbar;
title(['L_{min} = ' num2str(round(L_min/1000,1)) 'km ,W_{min} = ' num2str(round(W_min/1000,1)) 'km'])
ha(1).Position(1) = 0.05; ha(1).Position(2) = 0.8;
ha(1).Position(3) = xx;
ha(1).Position(4) = zz;


f = find(W./1000<11);
ff1 = [f L(f) W(f)];

k = 1;
axes(ha(2))

imagesc([L1{k,1}],[W1{k,1}],slip_lw{k,pl});  colormap(slipcolor);   colorbar;
ha(2).Position(1) = ha(1).Position(1) + ha(1).Position(3) + 0.05;
ha(2).Position(2) = ha(1).Position(2) ; %+ ha(1).Position(3);
ha(2).Position(3) = (xx)*(L(k)/(L_min)); %+ ha(1).Position(3);
ha(2).Position(4) = zz*(W(k)/(W_min));
title(['L = ' num2str(round(L(k)/1000,1)) 'km ,W = ' num2str(round(W(k)/1000,1)) 'km' ])
[row, col] = find(A_new{k,1} ~= 0);
hold on;
plot( [ L1{k,1}(min(col)) L1{k,1}(max(col)) L1{k,1}(max(col)) L1{k,1}(min(col)) L1{k,1}(min(col)) ],...
[ W1{k,1}(min(row)) W1{k,1}(min(row)) W1{k,1}(max(row)) W1{k,1}(max(row)) W1{k,1}(min(row)) ],'k','LineWidth',1);


k = 2;
axes(ha(3))
imagesc([L1{k,1}],[W1{k,1}],slip_lw{k,pl});  colormap(slipcolor);   colorbar;
ha(3).Position(1) = ha(2).Position(1) + ha(2).Position(3) + 0.05;
ha(3).Position(2) = ha(1).Position(2) ; %+ ha(1).Position(3);
ha(3).Position(3) = (xx)*(L(k)/(L_min)); %+ ha(1).Position(3);
ha(3).Position(4) = zz*(W(k)/(W_min));
title(['L = ' num2str(round(L(k)/1000,1)) 'km ,W = ' num2str(round(W(k)/1000,1)) 'km' ])
[row, col] = find(A_new{k,1} ~= 0);
hold on;
plot( [ L1{k,1}(min(col)) L1{k,1}(max(col)) L1{k,1}(max(col)) L1{k,1}(min(col)) L1{k,1}(min(col)) ],...
[ W1{k,1}(min(row)) W1{k,1}(min(row)) W1{k,1}(max(row)) W1{k,1}(max(row)) W1{k,1}(min(row)) ],'k','LineWidth',1);





%%%% second row of plots

k= 3;
axes(ha(4))
imagesc([L1{k,1}],[W1{k,1}],slip_lw{k,pl});  colormap(slipcolor);   colorbar;
ha(4).Position(1) = ha(1).Position(1);
ha(4).Position(2) = 0.65; %+ ha(1).Position(3);
ha(4).Position(3) = (xx)*(L(k)/(L_min)); %+ ha(1).Position(3);
ha(4).Position(4) = zz*(W(k)/(W_min));
title(['L = ' num2str(round(L(k)/1000,1)) 'km ,W = ' num2str(round(W(k)/1000,1)) 'km' ])
[row, col] = find(A_new{k,1} ~= 0);
hold on;
plot( [ L1{k,1}(min(col)) L1{k,1}(max(col)) L1{k,1}(max(col)) L1{k,1}(min(col)) L1{k,1}(min(col)) ],...
[ W1{k,1}(min(row)) W1{k,1}(min(row)) W1{k,1}(max(row)) W1{k,1}(max(row)) W1{k,1}(min(row)) ],'k','LineWidth',1);



k= 4;
axes(ha(5))
imagesc([L1{k,1}],[W1{k,1}],slip_lw{k,pl});  colormap(slipcolor);   colorbar;
ha(5).Position(1) = ha(4).Position(1) + ha(4).Position(3) + 0.1;
ha(5).Position(2) = ha(4).Position(2) ;%+ ha(1).Position(3);
ha(5).Position(3) = (xx)*(L(k)/(L_min)); %+ ha(1).Position(3);
ha(5).Position(4) = zz*(W(k)/(W_min));
title(['L = ' num2str(round(L(k)/1000,1)) 'km ,W = ' num2str(round(W(k)/1000,1)) 'km' ])
[row, col] = find(A_new{k,1} ~= 0);
hold on;
plot( [ L1{k,1}(min(col)) L1{k,1}(max(col)) L1{k,1}(max(col)) L1{k,1}(min(col)) L1{k,1}(min(col)) ],...
[ W1{k,1}(min(row)) W1{k,1}(min(row)) W1{k,1}(max(row)) W1{k,1}(max(row)) W1{k,1}(min(row)) ],'k','LineWidth',1);



end

set(gcf,'position',[0,0,1500,1500])

% saveas(gcf,['2D_GSD_gaussian.png'])


%%
% figure(2)
% [ha, pos] = tight_subplot(3,3,[0.05 0.05],[0.02 0.02],[0.03, 0.03]);
% set(gcf,'position',[50,50,1800,1300])
% 
% 
% for i=1:9
%       axes(ha(i))
%       imagesc([L1{i,1}],[W1{i,1}],slip_lw{i});  colormap(slipcolor); axis tight equal   
% end
