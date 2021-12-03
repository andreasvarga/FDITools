% MATLAB commands to solve the AMDP in Example 6.2

% Uses the Control System Toolbox and Descriptor Systems Tools (DSTOOLS)

% Lateral aircraft model without faults
A = [-.4492 0.046 .0053 -.9926;
     0 0 1 0.0067;
     -50.8436 0 -5.2184 .722;
     16.4148 0 .0026 -.6627];
Bu = [0.0004 0.0011; 0 0; -1.4161 .2621; -0.0633 -0.1205];
[n,mu] = size(Bu); p = 2; mw = n+p; m = mu+mw; 
Bw = eye(n,mw);
C = 180/pi*eye(p,n);  Du = zeros(p,mu); Dw = [zeros(p,n) eye(p)];
% define the LOE faults Gamma_i
Gamma = 1 - [ 0  0 0 .5 .5 .5 1  1 1;
            0 .5 1  0 .5  1 0 .5 1 ]';
N = size(Gamma,1);
% define multiple physical fault model Gui = Gu1*Gamma_i, Gwi = Gw1
sysuw = ss(zeros(p,m,N,1));
for i=1:N
    sysuw(:,:,i,1) = ss(A,[Bu*diag(Gamma(i,:)) Bw],C,[Du Dw]);
end

% optimal H-inf design (standard case)
Q1 = [eye(p) -sysuw(:,1:mu)];

% perform optimal synthesis
R = ss(zeros(p,mu+mw,N,N)); Q = ss(zeros(p,p+mu,N,1));
tol = 1.e-7;
for i = 1:N
  rwi = gir(Q1(:,1:p,i,1)*sysuw(:,mu+1:m,i,1),tol);
  [gi,go] = goifac(rwi,struct('tol',1.e-7));
  Q(:,:,i,1) = gminreal(go\Q1(:,:,i,1),tol);
  for j = 1:N
    R(:,:,i,j) = gir(Q(:,:,i,1)*[sysuw(:,:,j,1); eye(mu,m)],tol);
  end
end

% scale Qi and Rij; determine gap 
distinf = norm(R(:,1:mu),inf); 
beta = zeros(N,1);
for i=1:N
  scale = min(distinf(i,[1:i-1 i+1:N]));
  distinf(i,:) = distinf(i,:)/scale;
  Q(:,:,i,1) = Q(:,:,i,1)/scale;
  for j = 1:N
     R(:,:,i,j) = R(:,:,i,j)/scale;
  end
  beta(i) = scale;
end
gap = beta

%%  MATLAB commands to generate the figures for Example 6.2
% 
figure, mesh(distinf)
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
    [r,t]=lsim(R(:,:,j,i),u,t);
    % use a Narendra filter with (alpha,beta,gamma) = (0.9,0.1,10)
    alpha = 0.9; beta = 0.1; gamma = 10;
    theta=alpha*sqrt(r(:,1).^2+r(:,2).^2)+...
        beta*sqrt(lsim(tf(1,[1 gamma]),r(:,1).^2+r(:,2).^2,t));
    plot(t,theta), 
    if i == 1 title(['Model ',num2str(j)]), end
    if i == j, ylim([0 1]), end 
    if j == 1, ylabel(['\theta_', num2str(i)],'FontWeight','bold'), end
    if i == N & j == 5, xlabel('Time (seconds)','FontWeight','bold'), end
    k = k+N;
  end
end
