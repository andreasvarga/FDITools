% MATLAB commands to solve the MDP in Example 6.1
% Uses the Control System Toolbox and DSTOOLS V0.5 (and later)

clear variables
% Lateral aircraft model without faults
A = [-.4492 0.046 .0053 -.9926;
     0 0 1 0.0067;
     -50.8436 0 -5.2184 .722;
     16.4148 0 .0026 -.6627];
Bu = [0.0004 0.0011; 0 0; -1.4161 .2621; -0.0633 -0.1205]; 
C = eye(4); p = size(C,1); mu = size(Bu,2); 
% define the LOE faults Gamma_i
Gamma = 1 - [ 0  0 0 .5 .5 .5 1  1 1;
            0 .5 1  0 .5  1 0 .5 1 ]';
N = size(Gamma,1);
% define multiple physical fault model Gui = Gu*Gamma_i
sysu = ss(zeros(p,mu,N,1));
for i=1:N
    sysu(:,:,i,1) = ss(A,Bu*diag(Gamma(i,:)),C,0);
end

% setup initial full order model detector Q1i = [I -Gui]
Q1 = [eye(p) -sysu];

% solve a minimum dynamic cover problem
% select row 1 of h*Q1i to combine it with rows 1:4
% of Q1i; the result is a least order Qi = Q2i*Q1i
h = [ 0.7645   0.8848   0.5778   0.9026];
tol = 1.e-7;                     % set tolerance
Q = ss(zeros(1,p+mu,N,1));
for i = 1:N-1
  Q(:,:,i,1) = glmcover1([h;eye(p)]*Q1(:,:,i,1),1,tol);  
  %Q(:,:,i,1) = glcf(temp,struct('sdeg',-1)); 
end
Q(1,1:p,N,1) = Q(1,1:p,1,1);

% compute Rij and their norms
R = ss(zeros(1,mu,N,N));
for i = 1:N
  for j = 1:N
    temp = Q(:,:,i,1)*[sysu(:,:,j,1); eye(mu)];
    R(:,:,i,j) = gir(temp,tol);
  end
end

% scale Qi and Rij
distinf = norm(R,inf);
for i=1:N
  scale = 1/min(distinf(i,[1:i-1 i+1:N]));
  Q(:,:,i,1) = scale*Q(:,:,i,1);
  for j = 1:N
     R(:,:,i,j) = scale*R(:,:,i,j);
  end
end

%% MATLAB commands to generate the figures for Example 6.1

R1 = ss(zeros(p,mu,N,N));
for i = 1:N
    for j = 1:N
        R1(:,:,i,j) = Q1(:,:,i,1)*[sysu(:,:,j,1); eye(mu)];
    end
end

distinf1=norm(R1,inf);
figure,mesh(distinf1)
colormap hsv
% colormap winter
% colormap parula
title('Norms of initial residual models')
ylabel('Residual numbers')
xlabel('Model numbers')

%%
distinf=norm(R,inf);
figure,mesh(distinf)
colormap hsv
% colormap winter
% colormap parula
title('Norms of final residual models')
ylabel('Residual numbers')
xlabel('Model numbers')


%%
figure
k1 = 0;
for j = 1:N
  k1 = k1+1; 
  k = k1;
  for i=1:N 
    subplot(N,N,k) 
    [r,t]=step(R(:,:,j,i),4); 
    plot(t,r(:,:,1),t,r(:,:,2)), 
    if i == 1, title(['Model ',num2str(j)]), end
    if i == j, ylim([-1 1]), end 
    if j == 1, ylabel(['r^(^', num2str(i),'^)'],'FontWeight','bold'), end
    if i == N && j == 5, xlabel('Time (seconds)','FontWeight','bold'), end
    k = k+N;
  end
end



