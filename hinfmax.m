function [gamma,ind] = hinfmax(sys,freq)
%HINFMAX  Maximum of H-inf norms of columns of a transfer function matrix 
%  [GAMMA,IND] = HINFMAX(SYS) computes GAMMA, the maximum of the H-infinity
%  norms of the columns of the transfer function matrix G of the stable LTI 
%  system SYS. IND is the index of the column for which the maximum is
%  achieved. 
%
%  [GAMMA,IND] = HINFMAX(SYS,FREQ) computes for a real vector of frequency 
%  values FREQ, GAMMA, the maximum of the 2-norms of the columns of the 
%  frequency responses of G evaluated for all values contained in FREQ. 
%  IND is the index of the column for which the maximum is achieved. 
%
%  Note: The stability of SYS is not checked.

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 11-06-2018.
%  Revisions: 

narginchk(1,2)
nargoutchk(0,2)

% check input system form
if ~isa(sys,'lti')
   error('The input system SYS must be an LTI system object')
end
m = size(sys,2);

if m == 0
   gamma = []; ind = [];
   return
end

if nargin == 1
    freq = [];
else
    if ~isempty(freq)
        validateattributes(freq, {'double'},{'real','vector','>=',0},'','FREQ',2)
    end
end

gamma = 0; ind = 1;
if ~isempty(freq) 
   % evaluate GAMMA as the maximum of the norms of columns of the frequency 
   % responses of G evaluated over all frequencies contained in FREQ 
   if isequal(freq,0) 
      % use DCGAIN if only frequency 0 is present
      frgains = dcgain(sys); 
      for j = 1:m
          temp = norm(frgains(:,j));
          if gamma < temp
             gamma = temp; ind = j;
          end             
      end
   else
      for i = 1:length(freq)
          frgains = freqresp(sys,freq(i)); 
          for j = 1:m
              temp = norm(frgains(:,j));
              if gamma < temp
                 gamma = temp; ind = j;
              end             
          end
      end  
   end
else 
   % evaluate BETA as the minimum of H-infinity norms of the columns of G
   for j = 1:m
       temp = norm(sys(:,j),inf);
       if gamma < temp
          gamma = temp; ind = j;
       end             
   end
end

% end HINFMAX
end
