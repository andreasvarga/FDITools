% Example 7.3 - Nullspace-based synthesis
% Uses the Control System Toolbox and Descriptor System Tools

% define the state space realizations of Gu, Gd and Gf
A = diag([ -2 -3 -1 2 3]);
Bu = [1 1 0 0 0]'; Bd = [0 0 2 0 0]';
Bf = [0 0 0 2 2;0 0 0 0 0]';
C = [-1 0 -1 1.5 0; 0 -1 0 0 2.5];
Du = [ 1 1]'; Dd = [1 0]'; Df = [1 0; 1 1];
p = 2; mu = 1; md = 1; mf = 2;       % enter dimensions
sys = ss(A,[Bu Bd Bf],C,[Du Dd Df]); % define system

% compute [Q Rf], where Q is a left nullspace basis of
% [Gu Gd;I 0] and Rf = Q*[Gf;0];
Q_Rf = glnull([sys;eye(mu,mu+md+mf)],struct('m2',mf));

% compute a stable left coprime factorization of [Q Rf] using
% explicitly computed output injection matrix
[al,b,cl,d,el] = dssdata(Q_Rf);
k = gsfstab(al',el',cl',-3,-2).';
M = dss(al+k*cl,k,cl,1,el); 
Q_Rf = dss(al+k*cl,b+k*d,cl,d,el);    % Q_Rf <- M*Q_Rf

% alternative computation (comment out next line)
% [Q_Rf,M] = glcf(Q_Rf,struct('sdeg',-3,'smarg',-2));


% compute Q and Rf; scale to meet example
Q = sqrt(2)*Q_Rf(:,1:p+mu); Rf = sqrt(2)*Q_Rf(:,p+mu+1:end);

% display results
minreal(zpk(Q)), minreal(zpk(Rf)), minreal(zpk(M))
