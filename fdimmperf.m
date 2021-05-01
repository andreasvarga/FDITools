function gamma = fdimmperf(R,sysr,nrmflag,S)
%FDIMMPERF  Model-matching performance of FDI filters
%  GAMMA = FDIMMPERF(R,SYSR,NRMFLAG) for LTI systems R and SYSR and 
%  NRMFLAG = Inf or NRMFLAG = 2, computes GAMMA = ||R - SYSR||_NRMFLAG,  
%  NRMFLAG-norm of the model-matching error between the LTI systems R   
%  and SYSR, where R has the input-output form
%
%       y = Ru*u + Rd*d + Rf*f + Rw*w                             
%
%  and SYSR has the input-output form
%
%       yr = Mru*u + Mrd*D + Mrf*f + Mrw*w                        
%
%  with the Laplace- or Z-transformed plant outputs y and yr, control 
%  inputs u, disturbance inputs d, fault inputs f and noise inputa w, 
%  and with Ru, Rd, Rf and Rw, respectively, Mru, Mrd, Mrf and Mrw,
%  the corresponding transfer function  matrices. GAMMA is computed as
%  the H-infinity norm 
%
%     GAMMA = ||[Ru-Mru Rd-Mrd Rf-Mrf Rw-Mrw]||_NRMFLAG
%
%  The inputs u, d, f and w of R and SYSR correspond to the standard 
%  input groups named 'controls', 'disturbances', 'faults' and 'noise', 
%  respectively. Each of these groups can be void, in which case the 
%  corresponding transfer function matrix is assumed zero. 
%
%  If R and SYSR are N x 1 arrays of LTI systems, R{i}, i = 1, ..., N,  
%  and SYSR{i}, i = 1, ..., N, respectively, with the i-th system R{i} 
%  having the input-output form
%
%       y_i = Ru_i*u + Rd_i*d + Rf_i*f + Rw_i*w ,                    (1)
%
%  and the i-th system SYSR{i} having the input-output form
%
%       yr_i = Mru_i*u + Mrd_i*d + Mrf_i*f + Mrw_i*w ,                (2)
%
%  then GAMMA is an N-dimensional vector whose i-th component is
%
%   GAMMA(i) = ||[Ru_i-Mru_i Rd_i-Mrd_i Rf_i-Mrf_i Rw_i-Mrw_i]||_NRMFLAG
%
%  GAMMA = FDIMMPER(R,SYSR) is equivalent to the call 
%  GAMMA = FDIMMPER(R,SYSR,Inf). 
%
%  GAMMA = FDIMMPER(R,[],NRMFLAG), for R a LTI system, computes
%  GAMMA = ||Gw||_NRMFLAG. 
%  If R is an N x 1 array of LTI systems, R{i}, i = 1, ..., N, where 
%  the i-th system R{i} has an input-output form in (1), then GAMMA is 
%  an N-dimensional vector whose i-th component is 
%
%     GAMMA(i) = ||Gw_i||_NRMFLAG  .
%
%  GAMMA = FDIMMPER(R) is equivalent to the call
%  GAMMA = FDIMMPER(R,[],Inf).  
%
%  GAMMA = FDIMMPER(R,[],NRMFLAG,S) computes for a LTI system R and 
%  a q x mf logical structure matrix S, the performance value 
%  GAMMA = ||[Rfc Rw]||_NRMFLAG, where Rfc is formed as follows: 
%  the (i,j)-th element of Rfc is the (i,j)-th element of Rf if 
%  S(i,j) = false and zero otherwise. 
%  If R is an N x 1 array of stable LTI systems, with R{i} as in (1)
%  and S is an N x mf logical structure matrix, then GAMMA is an
%  N-dimensional vector whose i-th component is 
%
%     GAMMA(i) = ||[Rfc_i Rw_i]||_NRMFLAG, 
%
%  where Rfc_i is formed as follows: the j-th column of Rfc_i is the j-th 
%  column of Rf_i if S(i,j) = false and zero otherwise. 

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 13-08-2018.
%  Revisions: A. Varga, 16-10-2018.

narginchk(1,4)
nargoutchk(0,1)

if nargin < 2
    sysr = [];
end
if nargin < 3
   nrmflag = inf;
else
   if nrmflag ~= 2 && nrmflag ~= inf
      error('The third input argument must be set to 2 or Inf')
   end
end       

% check input system form
if isa(R,'cell')
   N = length(R);
   validateattributes(R, {'cell'},{'vector'},'','R')
   % select the indices of non-empty systems 
   syssel = zeros(1,N);
   for i = 1:N
       if ~isempty(R{i}) 
          if ~isa(R{i},'lti')
             error('R must be cell array of LTI system objects')
          else
             syssel(i) = i;
          end 
       end 
   end
   syssel = syssel(syssel > 0); 
   if isempty(syssel)
      error('R must be a non-empty cell array of LTI system objects')
   end
   j = syssel(1); 
   Ts = R{j}.Ts;
   if isfield(R{j}.InputGroup,'controls')
      % controls
      inpu = R{j}.InputGroup.controls;
      mu = length(inpu);
   else
      mu = 0; 
      inpu = []; 
   end    
   if isfield(R{j}.InputGroup,'disturbances')
      % disturbances 
      inpd = R{j}.InputGroup.disturbances;
      md = length(inpd);
   else
      md = 0; 
      inpd = []; 
   end    
   if isfield(R{j}.InputGroup,'faults')
      % faults
      inpf = R{j}.InputGroup.faults;
      mf = length(inpf);
   else
      mf = 0; 
      inpf = []; 
   end    
   if isfield(R{j}.InputGroup,'noise')
      % noise 
      inpw = R{j}.InputGroup.noise;
      mw = length(inpw);
   else
      mw = 0; 
      inpw = []; 
   end    
   for i = syssel(2:end)
       if Ts ~= R{i}.Ts
          error('All component models of R must have the same sampling time')
       end
       if isfield(R{i}.InputGroup,'controls')
          % controls
          inpui = R{i}.InputGroup.controls;
          mui = length(inpu);
       else
          mui = 0; 
          inpui = []; 
       end    
       if mu ~= mui || ~isequal(inpu,inpui)
          error('All component models of R must have the same number of control inputs')
       end
       if isfield(R{j}.InputGroup,'disturbances')
          % disturbances 
          inpdi = R{i}.InputGroup.disturbances;
          mdi = length(inpd);
       else
          mdi = 0; 
          inpdi = []; 
       end    
       if md ~= mdi || ~isequal(inpd,inpdi)
          error('All component models of R must have the same number of disturbance inputs')
       end
       if isfield(R{i}.InputGroup,'faults')
          % fault inputs
          inpfi = R{i}.InputGroup.faults;
          mfi = length(inpfi);
       else
          mfi = 0; 
          inpfi = []; 
       end    
       if mf ~= mfi || ~isequal(inpf,inpfi)
          error('All component models of R must have the same number of fault inputs')
       end
       if isfield(R{i}.InputGroup,'noise')
          % noise inputs
          inpwi = R{i}.InputGroup.noise;
          mwi = length(inpwi);
       else
          mwi = 0; 
          inpwi = []; 
       end    
       if mw ~= mwi || ~isequal(inpw,inpwi)
          error('All component models of R must have the same number of noise inputs')
       end
   end
   if ~isempty(sysr)
      validateattributes(sysr, {'cell'},{'vector'},'','SYSR')
      if N ~= length(sysr)
         error('R and SYSR must be cell arrays of the same lenght')
      end
      % select the indices of non-empty systems 
      sysrsel = zeros(1,N);
      for i = 1:N
          if ~isempty(sysr{i}) 
             if ~isa(sysr{i},'lti')
                error('SYSR must be a cell array of LTI system objects')
             else
                sysrsel(i) = i;
             end 
          end 
      end
      sysrsel = sysrsel(sysrsel > 0); 
      if ~isequal(syssel,sysrsel)
         error('R and SYSR must contain the same empty cells')
      end
      j = syssel(1); 
      Trs = sysr{j}.Ts;
      if isfield(sysr{j}.InputGroup,'controls')
         % controls
         inpru = sysr{j}.InputGroup.controls;
         mru = length(inpru);
      else
         mru = 0; 
         inpru = []; 
      end    
      if isfield(sysr{j}.InputGroup,'disturbances')
         % disturbances 
         inprd = sysr{j}.InputGroup.disturbances;
         mrd = length(inprd);
      else
         mrd = 0; 
         inprd = []; 
      end    
      if isfield(sysr{j}.InputGroup,'faults')
         % faults
         inprf = sysr{j}.InputGroup.faults;
         mrf = length(inprf);
      else
         mrf = 0; 
         inprf = []; 
      end    
      if isfield(sysr{j}.InputGroup,'noise')
         % noise
         inprw = sysr{j}.InputGroup.noise;
         mrw = length(inpw);
      else
         mrw = 0; 
         inprw = []; 
      end    
      for i = syssel(2:end)
          if Trs ~= sysr{i}.Ts
             error('All component models of SYSR must have the same sampling time')
          end
          if isfield(sysr{i}.InputGroup,'controls')
             % controls
             inprui = sysr{i}.InputGroup.controls;
             mrui = length(inprui);
          else
             mrui = 0; 
             inprui = []; 
          end    
          if mru ~= mrui || ~isequal(inpru,inprui)
             error('All component models of SYSR must have the same number of control inputs')
          end
          if isfield(sysr{i}.InputGroup,'disturbances')
             % disturbances 
             inprdi = sysr{i}.InputGroup.disturbances;
             mrdi = length(inprdi);
          else
             mrdi = 0; 
             inprdi = []; 
          end    
          if mrd ~= mrdi || ~isequal(inprd,inprdi)
             error('All component models of SYSR must have the same number of disturbance inputs')
          end
          if isfield(sysr{i}.InputGroup,'faults')
             % fault inputs
             inprfi = sysr{i}.InputGroup.faults;
             mrfi = length(inprfi);
          else
             mrfi = 0; 
             inprfi = []; 
          end    
          if mrf ~= mrfi || ~isequal(inprf,inprfi)
             error('All component models of SYSR must have the same number of fault inputs')
          end
          if isfield(sysr{i}.InputGroup,'noise')
             % noise inputs
             inprwi = sysr{i}.InputGroup.noise;
             mrwi = length(inprwi);
          else
             mrwi = 0; 
             inprwi = []; 
          end    
          if mrw ~= mrwi || ~isequal(inprw,inprwi)
             error('All component models of SYSR must have the same number of noise inputs')
          end
      end
      if Trs ~= Ts
          error('R and SYSR must have the same sampling time')
      end
      if mru && mu ~= mru 
         error('Incompatible control input groups between R and SYSR')
      end
      if mrd && md ~= mrd
         error('Incompatible disturbance input groups between R and SYSR')
      end
      if mrf && mf ~= mrf  
         error('Incompatible fault input groups between R and SYSR')
      end
      if mrw && mw ~= mrw  
         error('Incompatible noise input groups between R and SYSR')
      end
   end
elseif isa(R,'lti')
   [~,~,N] = size(R);
   if N ~= 1
       error('No multiple models supported')
   end
   if isfield(R.InputGroup,'controls')
      % controls
      inpu = R.InputGroup.controls;
      mu = length(inpu);  
   else
      inpu = []; mu = 0;
   end
   if isfield(R.InputGroup,'disturbances')
      % disturbances
      inpd = R.InputGroup.disturbances;
      md = length(inpd);  
   else
      inpd = []; md = 0;
   end
   if isfield(R.InputGroup,'faults')
      % faults
      inpf = R.InputGroup.faults;
      mf = length(inpf);
   else
      mf = 0;
      inpf = []; 
   end    
   if isfield(R.InputGroup,'noise')
      % faults
      inpw = R.InputGroup.noise;
      mw = length(inpw);
   else
      mw = 0;
      inpw = []; 
   end    
   if ~isempty(sysr)
      if sysr.Ts ~= R.Ts
          error('R and SYSR must have the same sampling time')
      end
      % decode input information for SYSR 
      if isfield(sysr.InputGroup,'controls')
         % controls
         inpru = sysr.InputGroup.controls;
         mru = length(inpru);
         if mu ~= mru 
            error('Incompatible control input groups between R and SYSR')
         end
      else
         inpru = []; 
      end
      if isfield(sysr.InputGroup,'disturbances')
         % disturbances
         inprd = sysr.InputGroup.disturbances;
         mrd = length(inprd);
         if md ~= mrd
            error('Incompatible disturbance input groups between R and SYSR')
         end
      else
         inprd = [];
      end
      if isfield(sysr.InputGroup,'faults')
         % faults
         inprf = sysr.InputGroup.faults;
         mrf = length(inprf);
         if mf ~= mrf  
            error('Incompatible fault input groups between R and SYSR')
         end
      else
         inprf = []; 
      end
      if isfield(sysr.InputGroup,'noise')
         % noise
         inprw = sysr.InputGroup.noise;
         mrw = length(inprw);
         if mw ~= mrw  
            error('Incompatible noise input groups between R and SYSR')
         end
      else
         inprw = []; 
      end
   end
else
   error('R and SYSR must be either cell arrays or LTI system objects')
end

if nargin < 4
   S = [];
else
   if ~isempty(S)
       validateattributes(S, {'logical'},{'binary'},'','S',4)
       if mf ~= size(S,2)
          error('Structure matrix incompatible with the number of faults')
       end
       nb = size(S,1);
       if N == 1
          if size(R,1) ~= nb 
             error('Structure matrix incompatible with the number of outputs')
          end
       else
          if N ~= nb 
             error('Structure matrix incompatible with the number of systems')
          end
       end
   end
end

m = mu+md+mf+mw;
if m == 0
   gamma = zeros(N,1);
   return
end


if isempty(sysr)
   if isa(R,'cell')
      gamma = zeros(N,1);  
      if isempty(S) 
         for i = syssel
             gamma(i) = norm(R{i}(:,inpw),nrmflag); 
         end
      else
         for i = syssel
             gamma(i) = norm(R{i}(:,[inpf(~S(i,:)) inpw]),nrmflag); 
         end
     end
   else
      if isempty(S) 
         gamma = norm(R(:,inpw),nrmflag);
      else
         sysref = ss(zeros(nb,mf));  sysref.Ts = R.Ts; 
         for i = 1:nb
             inpwi = inpf(~S(i,:));
             sysref(i,inpwi) = R(i,inpwi);
         end
         gamma = norm([sysref R(:,inpw)],nrmflag);     
      end
   end
   return
else
   if ~isempty(S)
       warning('Fourth input argument S ignored')
   end
end



rinp = zeros(0,m);
if ~isempty(inpru) 
   rinp = [rinp;eye(mu,m)];
end
if ~isempty(inprd) 
   rinp = [rinp; zeros(md,mu) eye(md,m-mu)];
end
if ~isempty(inprf) 
   rinp = [rinp; zeros(mf,mu+md) eye(mf,m-mu-md)];
end
if ~isempty(inprw) 
   rinp = [rinp; zeros(mw,m-mw) eye(mw)];
end 
if isa(R,'cell')
   gamma = zeros(N,1);  
   for i = syssel
       % substract explicitly [ Gru_i Grd_i Grf_i Grw_i ] 
       syst = R{i}(:,[inpu inpd inpf inpw])-sysr{i}(:,[inpru inprd inprf inprw])*rinp;
       gamma(i) = norm(syst,inf);
       if nrmflag == 2 
          if Ts == 0 
             % set negligible elements in the feedthrough matrix to zero 
             syst.d(abs(syst.d) <= sqrt(eps)*gamma(i)) = 0; 
          end
          gamma(i) = norm(prescale(syst),2);
       end
   end
else
   % substract explicitly [ Gru Grd Grf Grw ] 
   syst = R(:,[inpu inpd inpf inpw])-sysr(:,[inpru inprd inprf inprw])*rinp;
   gamma = norm(syst,inf);
   if nrmflag == 2 
      if R.Ts == 0 
         % set negligible elements in the feedthrough matrix to zero 
         syst.d(abs(syst.d) <= sqrt(eps)*gamma) = 0; 
      end
      gamma = norm(prescale(syst),2);
   end
end

% end FDIMMPERF
end
