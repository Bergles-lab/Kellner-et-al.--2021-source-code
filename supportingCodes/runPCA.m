%Kellner et al., 2021
function mve = runPCA(X, saveDirName)

X = imresize(X, 0.5);

[dFoF,Fo, X] = normalizeImg(X,10,1);

[m,n,T] = size(X);
t = (1:T)/10; % frame rate of 10 fps
X = double(X);
%% set-up parameters. 'defaults' in parentheses
param = [];
param.numA = 25*4; % number of components (50)
param.numDb = [1 1]; % number of components in sinusoid (7) - first number is for sinusoids, second number (optional) is for exponentials
param.num_iter = 50; % number of iterations to run factorization for (50)
param.display = true; % display outputs as code is running (true)
param.tau = [.55 .12]; % decay constant in seconds. 2nd number (optional) is the rise time ([.55 .12] for gcamp6s)
param.kZ = [80 15 1]*.1; % parameters for spatial fit (in general, bigger number = smoother spatial components)
param.kA = [2 0 2]/20*.25; % parameters for temporal fit (in general, bigger number = less accurate fits to small transients)
param.lam = .01; % global scaling parameter that is multiplied times kZ and kA
param.dt = mean(diff(t)); % assumes frame time is constant throughout trace
param.A_init = zeros(T,param.numA);
for k = 1:param.numA, param.A_init(k:param.numA:T,k) = 1; end
%% run Ben Haeffele's structured factorization algorithm
tic; [A,Z,B,Da,Db,P] = SFC_setup(X,param); toc
%% at this point, you can save data (Da is usually large so we don't save it since it's easy to remake)
%%% put in desired save_name first %%%
dims = [m n T];
%save([dname save_name '.mat'],'t','A','Z','B','Db','P','dims','crop_inds','Xnormalizer','fn')
%save([processed_dir save_name '.mat'],'t','A','Z','B','Db','P','dims','fn','dname','processed_dir','fn')

%% create variables to be used for analysis
imgs = reshape(X,m*n,T)';
F_rec = Da*A*Z'; %F_rec = Da*A*Z'; %Estimated transients
B_rec = Db*B'; %B_rec = Db*B'; %Estimated background
%F_rec2 = A*Z'; %F_rec2 = A*Z'; %Estimated spike times
%% plot first 25 spatial components
%for j = 1:25, subplot(5,5,j); imagesc(reshape(Z(:,j),m,n)); title(num2str(j)); axis off; axis image; end

%%
mve = reshape(permute(F_rec,[2 1]),m,n,T);
writeTif(single(mve),[saveDirName 'PCA.tif'],32);

end