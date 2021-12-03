%% Example 6.1 - Solution of an EMDP using EMDSYN

% Uses the Control System Toolbox, Descriptor Systems Tools (DSTOOLS V0.71) 
% and Fault Detection and Isolation Tools (FDITOOLS V1.0)

clear variables
% Lateral aircraft model without faults
A = [-.4492 0.046 .0053 -.9926;
       0    0     1     0.0067;
   -50.8436 0   -5.2184  .722;
    16.4148 0     .0026 -.6627];
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

% setup synthesis model
sysu = mdmodset(sysu,struct('controls',1:mu));

% call of EMDSYN with the options for stability degree -1 and pole -1 for
% the filters, tolerance and a design matrix H to form a linear combination
% of the left nullspace basis vectorsH = [ 0.7645 0.8848 0.5778 0.9026 ];
H = [ 0.7645 0.8848 0.5778 0.9026 ];
opt_emdsyn = struct('sdeg',-1,'poles',-1,'HDesign',H);
[Q,R,info] = emdsyn(sysu,opt_emdsyn); info.MDperf

%% MATLAB commands to generate the figures for Example 6.1

figure, 
mesh(info.MDperf)
colormap hsv
title('Norms of final residual models')
ylabel('Residual numbers')
xlabel('Model numbers')

%
figure
k1 = 0;
for j = 1:N
  k1 = k1+1; 
  k = k1;
  for i=1:N 
    subplot(N,N,k) 
    [r,t] = step(R{j,i},4); 
    plot(t,r(:,:,1),t,r(:,:,2)), 
    if i == 1, title(['Model ',num2str(j)]), end
    if i == j, ylim([-1 1]), end 
    if j == 1, ylabel(['r^(^', num2str(i),'^)'],'FontWeight','bold'), end
    if i == N && j == 5, xlabel('Time (seconds)','FontWeight','bold'), end
    k = k+N;
  end
end

