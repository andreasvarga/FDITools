function [beta,ind] = hinfminus(sys,freq)
%HINFMINUS  H-(infinity-) index of a stable transfer function matrix 
%  [BETA,IND] = HINFMINUS(SYS,FREQ) computes the H-(infinity-) index 
%  of the transfer function matrix G of the stable LTI system SYS.
%  If FREQ is not specified or is an empty array, then BETA
%  is the minimum H-infinity norm of the columns of G. If FREQ is a real
%  vector of frequency values, then BETA is the minimum of the norms of 
%  the columns of the frequency responses of G evaluated for all values 
%  contained in FREQ. IND is the index of column for which the minimum is
%  achieved. 
%
%  Note: The stability of SYS is not checked.

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 19-04-2018.
%  Revisions: 

narginchk(1,2)
nargoutchk(0,2)

% check input system form
if ~isa(sys,'lti')
   error('The input system SYS must be an LTI system object')
end
m = size(sys,2);

if m == 0
   beta = []; ind = [];
   return
end

if nargin == 1
    freq = [];
else
    if ~isempty(freq)
        validateattributes(freq, {'double'},{'real','vector','>=',0},'','FREQ',2)
    end
end

beta = inf; ind = 1;
if ~isempty(freq) 
   % evaluate BETA as the minimum of the norms of columns of the frequency 
   % responses of G evaluated over all frequencies contained in FREQ 
   if isequal(freq,0) 
      % use DCGAIN if only frequency 0 is present
      frgains = dcgain(sys); 
      for j = 1:m
          temp = norm(frgains(:,j));
          if beta > temp
             beta = temp; ind = j;
          end             
      end
   else
      for i = 1:length(freq)
          frgains = freqresp(sys,freq(i)); 
          for j = 1:m
              temp = norm(frgains(:,j));
              if beta > temp
                 beta = temp; ind = j;
              end             
          end
      end  
   end
else 
   % evaluate BETA as the minimum of H-infinity norms of the columns of G
   for j = 1:m
       temp = norm(sys(:,j),inf);
       if beta > temp
          beta = temp; ind = j;
       end             
   end
end

% end HINFMINUS
end
