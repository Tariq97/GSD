function [A_new,A_stoch,realz_stoch] = GSD(A_new1,ssd,stoch,n_stoch,samp,acf)

%%% For a given 2D field with zeros around, the zeros are filled follwoing
%%% Generalized stochastic subdivision (Lewis 1985). Currently supports
%%% Gaussian and Von Karman ACFs

% Uses SpecSyn3, modified version of SpecSyn2 of Martin Mai.

% inputs
% 1. A_new -  slip field with zeros around it. Has to be a cell vector with
% each cell a matrix
% 2. ssd - structure with details of VK function.
% 3. A_stoch - pseudo fields.! size same as of A_new.
% 4. samp - sampling size of [dz,dx]


rng('default')
rng('shuffle')

ngh = ssd.neighbor_radius; corr = ssd.corr; %sigma = ssd.sigma; % incase random noise is added.
stoch = round(stoch);
for i=1:n_stoch
%     Rseed(1) = sum(100*clock);
%     rand('seed',Rseed(1));
    
   A_new{i,1} = A_new1;
   realz_stoch(i,1) = rand(1)*1000;
   
   % THe following part is added as SpecSyn2 has limitations due to number
   % of points.
   
if mod(stoch(2),samp(2)*1000) ~= 0 && mod(stoch(1),samp(1)*1000) == 0
    
    A_stoch{i,1} = SpecSyn3([stoch(1) stoch(2)+100],samp*1000,[corr(1)*1000 corr(2)*1000 corr(3)],acf,...
       realz_stoch(i),250,stoch(2)); 
 
elseif mod(stoch(1),samp(1)*1000) ~= 0 && mod(stoch(2),samp(2)*1000) == 0

     A_stoch{i,1} = SpecSyn3([stoch(1)+100 stoch(2)],samp*1000,[corr(1)*1000 corr(2)*1000 corr(3)],acf,...
       realz_stoch(i),250,stoch(2)); 
   
   
elseif mod(stoch(2),samp(2)*1000) ~= 0 && mod(stoch(1),samp(1)*1000) ~= 0
   
    A_stoch{i,1} = SpecSyn3([stoch(1)+100 stoch(2)+100],samp*1000,[corr(1)*1000 corr(2)*1000 corr(3)],acf,...
       realz_stoch(i),250,stoch(2)); 
else
    A_stoch{i,1} = SpecSyn3([stoch(1) stoch(2)],samp*1000,[corr(1)*1000 corr(2)*1000 corr(3)],acf,...
       realz_stoch(i),250,stoch(2)); 
   
end 
 
  % Balance the size of A_stoch to size of A_new
 if size(A_stoch{i,1},1) ~= size(A_new1,1)
    if size(A_stoch{i,1},1) > size(A_new1,1)        
        A_stoch{i,1}( size(A_stoch{i,1},1) -size(A_new1,1) ,:) = [];       
    else
         A_stoch{i,1} = [ A_stoch{i,1}; zeros(size(A_new1,1) - size(A_stoch{i,1},1),size(A_stoch{i,1},2))];       
    end   
 end

 if size(A_stoch{i,1},2) ~= size(A_new1,2)
   if size(A_stoch{i,1},2) > size(A_new1,2)  
       A_stoch{i,1}(:, size(A_stoch{i,1},2) -size(A_new1,2) ,:) = [];   
   else
        A_stoch{i,1} = [ A_stoch{i,1} (zeros(size(A_new1,2) - size(A_stoch{i,1},2),size(A_stoch{i,1},1)))'];  
   end
     
 end
   % Taper this A_stoch to avoid vertical artifacts in outputs.!!
   % But need A_stoch the size to that of A_new1
   
   A_stoch{i,1} = TaperSlip(A_stoch{i,1},[5,5,5],'hn');
   
   
end




% Step 1: Find all points with zero values (unknown points)
[unknown_i, unknown_j] = find(A_new{1,1} == 0);

% Find known points (non-zero points); (Make them just fault boundaries)
[known_i, known_j] = find(A_new{1,1} ~= 0);


min_i  = min(known_i); max_i = max(known_i);
min_j  = min(known_j); max_j = max(known_j);

known_i = [ [min_i:1:max_i]' ; ones((max_j-min_j),1).*max_i; [max_i-1:-1:min_i]'; ones((max_j-min_j),1).*min_i; ];
known_j = [ ones((max_i-min_i),1).*min_j; [min_j:1:max_j]'; ones((max_i-min_i),1).*max_j; [max_j-1:-1:min_j]' ]; 


lz = [0:samp(1):(size(A_new{1,1},1)-1)*(samp(1))]; 
lx = [0:samp(2):(size(A_new{1,1},2)-1)*(samp(2))];

[X,Z] = meshgrid(lx,lz);


%%%%%%%%%%%%%%%%%% Generalized stochastic subdivision loop below %%%%%%%%%%%%%%%%%%


az = corr(1); ax = corr(2); H = corr(3);
tic
k = 1;
while ~isempty(unknown_i)
    Rseed(1) = sum(100*clock);
    rand('seed',Rseed(1));
    %disp(length(unknown_i))

    coord_unknown = [ X(sub2ind(size(X), unknown_i, unknown_j))  Z(sub2ind(size(Z), unknown_i, unknown_j)) ];
    coord_known = [ X(sub2ind(size(X), known_i, known_j))  Z(sub2ind(size(Z), known_i, known_j)) ];
 
    
    % Create a KDTree for known points
    
    Mdl = KDTreeSearcher([coord_known(:,1), coord_known(:,2)]);
    
    % Instead of pdist2, use knnsearch to find the farthest points
    [idx, dists] = knnsearch(Mdl, coord_unknown);
    
    % Find the farthest index based on the distance
    [~, farthest_idx] = max(dists);



    p = [unknown_i(farthest_idx), unknown_j(farthest_idx) ]; % p is the farthest point

    %%%%%%%%%%%%%%%% Interpolation code below %%%%%%%%%

    % Step #1 Compute the distance matrix from point p    
    dist = sqrt( (X-X(p(1),p(2))).^2 + (Z-Z(p(1),p(2))).^2 );

    
    % Step #2 Find neighborhood point within threshold    
    [i,j] = find(dist<ngh);


    % Step 3: Select points with non-zero data
    validInd = A_new{1,1}(sub2ind(size(A_new{1,1}), i, j)) ~= 0;

    % Get the valid indices
    ind_i = i(validInd);
    ind_j = j(validInd);


    if isempty(ind_i)==1
       % disp('No points found in the neighborhood - Adding a random value')

        for ns=1:length(A_stoch)

           A_new{ns,1}(p(1),p(2)) = A_stoch{ns,1}(p(1),p(2)); %normrnd(0,sigma); %A_stoch{ns,1}(p(1),p(2));

        end

    else

       % disp([ num2str(length(ind_i)) ' points found in the neighborhood'])


    % Combine the indices into a single matrix
    ind = [ind_i, ind_j];

    % Step 4: Compute distances and b vector
    X_diff = (X(sub2ind(size(X), ind(:,1), ind(:,2))) - X(p(1), p(2))) / ax; 
    Z_diff = (Z(sub2ind(size(Z), ind(:,1), ind(:,2))) - Z(p(1), p(2))) / az;

    r = sqrt(X_diff.^2 + Z_diff.^2);

    if strcmp(acf,'gs')==1
        b = exp(-r.^2);
    elseif strcmp(acf,'ak')==1
     b = (2^(1 - H) / gamma(H)) .* (r .^ H) .* besselk(H, r);
     
    end
    b(r == 0) = 1;
    % Step 5: Compute pairwise differences in X and Z for all points (fully vectorized)
    X1a = (X(sub2ind(size(X), ind(:,1), ind(:,2))) - X(sub2ind(size(X), ind(:,1)', ind(:,2)'))) / ax; 
    Z1a = (Z(sub2ind(size(Z), ind(:,1), ind(:,2))) - Z(sub2ind(size(Z), ind(:,1)', ind(:,2)'))) / az; 


    % Step 6: Compute distance matrix and R (vectorized)
    
    r_matrix = sqrt(X1a.^2 + Z1a.^2);

    if strcmp(acf,'gs')==1
        R = exp(-r_matrix.^2);    % gaussian
    elseif strcmp(acf,'ak')==1
        R = (2^(1 - H) / gamma(H)) .* (r_matrix .^ H) .* besselk(H, r_matrix);
    end
    R(r_matrix == 0) = 1;

    a = pinv(R)*b;






    % 1) Symmetrize (kills roundoff-induced asymmetry)
    R = 0.5 * (R + R.');

    % 2) Add a small nugget for numerical stability (scale-aware)
    eps_reg = 1e-10 * max(1, max(abs(diag(R))));
    R = R + eps_reg * eye(size(R));

    % 3) Solve Ra = b stably (prefer Cholesky if possible)
    [LL,pchol] = chol(R,'lower');
    if pchol == 0
        % Cholesky succeeded: use it
        a = LL' \ (LL \ b);
        quad = b' * a;    % same as b' * (R \ b)
    else
        % Fall back to backslash if R not PD enough even after nugget
        a = R \ b;
        quad = b' * a;
    end
    
    % 4) Kriging variance; clamp tiny negatives to zero
    var_krig = max(0, 1 - quad);
    
    % 5) Standard deviation for normrnd
    sigma = sqrt(var_krig);
    
    
     % Add random noise to this.!!!
        % noise_variance = sqrt(1 - sum((a.*b)));   

        for ns=1:length(A_stoch)
            val1 = sum(a .* A_new{ns,1}(sub2ind(size(A_new{ns,1}), ind(:,1), ind(:,2)))); 

            A_new{ns,1}(p(1),p(2)) = val1 + normrnd(0,sigma);
            
        end

    end 

    %%%%%% Now add the calculated point to set of known points and remove from
    %%%%%% unknown points
    known_i = [known_i; p(1)];
    known_j = [known_j; p(2)];

    unknown_i(farthest_idx) = [];
    unknown_j(farthest_idx) = [];
    


end

toc

end