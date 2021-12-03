% MATLAB commands to solve the AMDP in Example 6.2 using Kalman filters

% Uses only the Control System Toolbox 

% Lateral aircraft model without faults
A = [-.4492 0.046 .0053 -.9926;
     0 0 1 0.0067;
     -50.8436 0 -5.2184 .722;
     16.4148 0 .0026 -.6627];
n = size(A,1); p = 2; mw = n+p;
Bu = [0.0004 0.0011; 0 0; -1.4161 .2621; -0.0633 -0.1205]; Bw = eye(n,n+p); 
C = 180/pi*eye(p,n); mu = size(Bu,2); Du = zeros(p,mu); Dw = [zeros(p,n) eye(p)];
% define the LOE faults Gamma_i
Gamma = 1 - [ 0  0 0 .5 .5 .5 1  1 1;
            0 .5 1  0 .5  1 0 .5 1 ]';
N = size(Gamma,1);
% define multiple physical fault model Gui = Gu*Gamma_i
sysuw = ss(zeros(p,mu+n+p,N,1));
for i=1:N
    sysuw(:,:,i,1) = ss(A,[Bu*diag(Gamma(i,:)) Bw],C,[Du Dw]);
end

% Kalman estimator
Qw = 0.0001*eye(4); Rv = 0.04*eye(p); 
for i = 1:N
  [Kest(:,:,i,1),L,P ] = kalman(sysuw(:,1:mu+n,i,1),Qw,Rv);
  Q_KF(:,:,i,1) = Kest(1:p,[mu+1:mu+p 1:mu],i,1);
  Q_KF(:,:,i,1).d = -eye(p,p+mu); 
end

% compute Rij and their norms
R_KF = ss(zeros(p,mu+mw,N,N)); 
tol=1.e-7;
for i = 1:N
  for j = 1:N
    temp = Q_KF(:,:,i,1)*[sysuw(:,:,j,1); eye(mu,mu+mw)];
    R_KF(:,:,i,j) = minreal(temp,tol);
  end
end
distinf_KF = norm(R_KF(:,1:mu),inf); 

% scale Qi and Rij with beta_i
gamma_KF = diag(norm(R_KF(:,mu+1:mu+mw),inf));
beta_KF = zeros(N,1); gap_KF = zeros(N,1);
for i=1:N
  scale = min(distinf_KF(i,[1:i-1 i+1:N]));
  distinf_KF(i,:) = distinf_KF(i,:)/scale;
  Q_KF(:,:,i,1) = Q_KF(:,:,i,1)/scale;
  for j = 1:N
     R_KF(:,:,i,j) = R_KF(:,:,i,j)/scale;
  end
  beta_KF(i) = scale;
  gap_KF(i) = scale/gamma_KF(i); 
end
beta_KF, gap_KF

%%  MATLAB commands to generate the figures as for Example 6.2
% 
figure, mesh(distinf_KF)
colormap hsv
% colormap winter
% colormap parula
title('Norms of final residual models')
ylabel('Residual numbers')
xlabel('Model numbers')

% generate input signals for Ex. 6.2
d = diag([ 1 1 0.01 0.01 0.01 0.01 0.03 0.03]);
t = (0:0.01:10)';  ns = length(t);
usin = gensig('sin',pi,t(end),0.01)+1.5;
usquare = gensig('square',pi*2,t(end),0.01)+0.3;
u = [ usquare usin (rand(ns,mw)-0.5)]*d;  
 
%
figure
k1 = 0;
for j = 1:N,
  k1 = k1+1; 
  k = k1;
  for i=1:N, 
    subplot(N,N,k), 
    [r,t]=lsim(R_KF(:,:,j,i),u,t);
    % use a Narendra filter with (alpha,beta,gamma) = (0.9,0.1,10)
    alpha = 0.9; beta = 0.1; gamma = 10;
    theta=alpha*sqrt(r(:,1).^2+r(:,2).^2)+...
        beta*sqrt(lsim(tf(1,[1 gamma]),r(:,1).^2+r(:,2).^2,t));
    plot(t,theta), 
    if i == 1 title(['Model ',num2str(j)]), end
    if i == j, ylim([0 1]), end 
    %if abs(i-j) < 2, ylim([0 2]), end 
    %ylim([1 10])
    if j == 1, ylabel(['\theta_', num2str(i)],'FontWeight','bold'), end
    if i == N & j == 5, xlabel('Time (seconds)','FontWeight','bold'), end
    k = k+N;
  end
end


%% MATLAB commands to solve the AMDP in Example 6.2 
%  Least order synthesis (LOS)
%  Uses the Control Toolbox and Descriptor Systems Tools

% optimal second order design with scalar output filters

% Lateral aircraft model without faults
A = [-.4492 0.046 .0053 -.9926;
     0 0 1 0.0067;
     -50.8436 0 -5.2184 .722;
     16.4148 0 .0026 -.6627];
n = size(A,1); p = 2; mw = n+p;
Bu = [0.0004 0.0011; 0 0; -1.4161 .2621; -0.0633 -0.1205]; Bw = eye(n,n+p); 
C = 180/pi*eye(p,n); mu = size(Bu,2); Du = zeros(p,mu); Dw = [zeros(p,n) eye(p)];
% define the LOE faults Gamma_i
Gamma = 1 - [ 0  0 0 .5 .5 .5 1  1 1;
            0 .5 1  0 .5  1 0 .5 1 ]';
N = size(Gamma,1);
% define multiple physical fault model Gui = Gu*Gamma_i
sysuw = ss(zeros(p,mu+n+p,N,1));
for i=1:N
    sysuw(:,:,i,1) = ss(A,[Bu*diag(Gamma(i,:)) Bw],C,[Du Dw]);
end

Q01 = [eye(p) -sysuw(:,1:mu)];
% solve a minimum dynamic cover problem
% select row 1 of h*Q1i to combine it with rows 1:2
% of Q1i; the result is a least order Qi = Q2i*Q1i
h2=[ 0.7645   0.8848];
tol = 1.e-7;                     % set tolerance
Q02 = ss(zeros(1,p+mu,N,1));
for i = 1:N-1
  temp = glmcover1([h2;eye(p)]*Q01(:,:,i,1),1,tol);  
  Q02(:,:,i,1) = glcf(temp,struct('sdeg',-1)); 
end
Q02(1,1:p,N,1) = Q02(1,1:p,1,1);

% perform optimal synthesis and compute Rij 
R_LOS = ss(zeros(1,mu+mw,N,N));
Q_LOS = ss(zeros(1,p+mu,N,1));
for i = 1:N
  rwi = gir(Q02(:,1:p,i,1)*sysuw(:,mu+1:mu+mw,i,1),tol);
  [gi,go] = goifac(rwi,struct('tol',1.e-7));
  Q_LOS(:,:,i,1) = gminreal(go\Q02(:,:,i,1),tol);
  for j = 1:N
    R_LOS(:,:,i,j) = gir(Q_LOS(:,:,i,1)*[sysuw(:,:,j,1); eye(mu,mu+mw)],tol);
  end
end

% scale Qi and Rij
distinf_LOS = norm(R_LOS(:,1:mu),inf); 
beta_LOS = zeros(N,1);
for i=1:N
  scale = min(distinf_LOS(i,[1:i-1 i+1:N]));
  distinf_LOS(i,:) = distinf_LOS(i,:)/scale;
  Q_LOS(:,:,i,1) = Q_LOS(:,:,i,1)/scale;
  for j = 1:N
     R_LOS(:,:,i,j) = R_LOS(:,:,i,j)/scale;
  end
  beta_LOS(i) = scale;
end
eta_LOS = beta_LOS


%%
% MATLAB commands to generate the figures as for Example 6.2

figure,mesh(distinf_LOS)
colormap hsv
% colormap winter
% colormap parula
title('Norms of internal filters')
ylabel('Residual numbers')
xlabel('Model numbers')

% generate input signals for Ex. 6.2
d = diag([ 1 1 0.01 0.01 0.01 0.01 0.03 0.03]);
t = (0:0.01:10)';  ns = length(t);
usin = gensig('sin',pi,t(end),0.01)+1.5;
usquare = gensig('square',pi*2,t(end),0.01)+0.3;
u = [ usquare usin (rand(ns,mw)-0.5)]*d;  
 
%
figure
k1 = 0;
for j = 1:N,
  k1 = k1+1; 
  k = k1;
  for i=1:N, 
    subplot(N,N,k), 
    [r,t]=lsim(R_LOS(:,:,j,i),u,t);
    % use a Narendra filter with (alpha,beta,gamma) = (0.2,0.1,10)
    alpha = 0.2; beta = 0.1; gamma = 10;
    theta=alpha*abs(r(:,1))+...
        beta*sqrt(lsim(tf(1,[1 gamma]),abs(r(:,1)),t));
    plot(t,theta), 
    if i == 1 title(['Model ',num2str(j)]), end
    if i == j, ylim([0 1]), end 
    %if abs(i-j) < 2, ylim([0 2]), end 
    %ylim([1 10])
    if j == 1, ylabel(['\theta_', num2str(i)],'FontWeight','bold'), end
    if i == N & j == 5, xlabel('Time (seconds)','FontWeight','bold'), end
    k = k+N;
  end
end



