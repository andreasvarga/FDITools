% Example 7.3 - Nullspace-based synthesis
% Uses the Control System Toolbox and Descriptor System Tools

% define the state-space realizations of Gu, Gd and Gf
n = 5; p = 3; mu = 3; md = 1; mf = 5;       % enter dimensions
% define matrices of the state space realization
A = [ 0 0 1.132 0 -1;
0 -0.0538 -0.1712 0 0.0705;
0 0 0 1 0;
0 0.0485 0 -0.8556 -1.013;
0 -0.2909 0 1.0532 -0.6859];
Bu = [0 0 0;-0.12 1 0;0 0 0;4.419 0 -1.665;1.575 0 -0.0732];
Bd = Bu(:,mu); Bf = [zeros(n,p) Bu(:,1:mu-1)];
C=eye(p,n); Du=zeros(p,mu); Dd=zeros(p,md); Df=eye(p,mf);
sys = ss(A,[Bu Bd Bf],C,[Du Dd Df]);        % define system

% compute [Q Rf], where Q is a left nullspace basis of
% [Gu Gd;I 0] and Rf = Q*[Gf;0];
[Q_Rf,info] = glnull([sys;eye(mu,mu+md+mf)],struct('m2',mf));
info.degs    % polynomial basis degrees

% determine 1-st and 2-nd order scalar output designs
Q1_Rf1 = glmcover1([0 1;eye(2)]*Q_Rf,1);    %  [Q1 Rf1]
Q2_Rf2 = glmcover1([1 1;eye(2)]*Q_Rf,1);    %  [Q2 Rf2]

% compute stable left coprime factorizations
opt_glcf = struct('sdeg',-1,'smarg',-2);
Q1_Rf1 = glcf(Q1_Rf1,opt_glcf);
Q2_Rf2 = glcf(Q2_Rf2,opt_glcf);

% compute Q1 and Rf1; check admissibility
Q1 = Q1_Rf1(:,1:p+mu); Rf1 = Q1_Rf1(:,p+mu+1:end);
g1 = abs(evalfr(Rf1,1i)), max(g1)/min(g1)
% compute Q2 and Rf2; check admissibility
Q2 = Q2_Rf2(:,1:p+mu); Rf2 = Q2_Rf2(:,p+mu+1:end);
g2 = abs(evalfr(Rf2,1i)), max(g2)/min(g2)

% display results
tf(Q1), tf(Q2)

