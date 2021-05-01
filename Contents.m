% FDITOOLS - Fault detection and isolation filter synthesis tools.
% Version 1.0.6           May 1, 2021
% Copyright 2016-2021 A. Varga
%
% Demonstration.
%   FDIToolsdemo - Demonstration of FDITOOLS.
%
% Setup of synthesis models.
%   fdimodset  - Setup of models for solving FDI synthesis problems.
%   mdmodset   - Setup of models for solving model detection synthesis problems.
%
% FDI related analysis.
%   fdigenspec - Generation of achievable FDI specifications.
%   fdichkspec - Feasibility analysis of a set of FDI specifications.
%
% Model detection related analysis.
%   mddist     - Computation of distances between component models. 
%   mddist2c   - Computation of distances to a set of component models. 
%
% Performance evaluation of FDI filters
%   fditspec   - Computation of the weak or strong structure matrix.
%   fdisspec   - Computation of the strong structure matrix.
%   fdifscond  - Fault sensitivity condition of FDI filters. 
%   fdif2ngap  - Fault-to-noise gap of FDI filters.
%   fdimmperf  - Model-matching performance of FDI filters.
%
% Performance evaluation of model detection filters
%   mdperf     - Distance mapping performance of model detection filters.
%   mdmatch    - Distance matching performance of model detection filters.
%   mdgap      - Noise gaps of model detection filters.
%   
% Synthesis of FDI filters.
%   efdsyn     - Exact synthesis of fault detection filters.
%   afdsyn     - Approximate synthesis of fault detection filters.
%   efdisyn    - Exact synthesis of fault detection and isolation filters.
%   afdisyn    - Approximate synthesis of fault detection and isolation filters.
%   emmsyn     - Exact model matching based synthesis of FDI filters.
%   ammsyn     - Approximate model matching based synthesis of FDI filters.
%
% Synthesis of model detection filters.
%   emdsyn     - Exact synthesis of model detection filters.
%   amdsyn     - Approximate synthesis of model detection filters.
%
% Miscellaneous. 
%   hinfminus  - H-(infinity-) index of a stable transfer function matrix. 
%   hinfmax    - Maximum of H-inf norms of columns of a transfer function matrix. 
%   efdbasesel - Selection of admissible basis vectors to solve the EFDP. 
%   afdbasesel - Selection of admissible basis vectors to solve the AFDP. 
%   emmbasesel - Selection of admissible basis vectors to solve the strong EFDIP. 
%   ammbasesel - Selection of admissible basis vectors to solve the strong AFDIP. 
%   emdbasesel - Selection of admissible basis vectors to solve the EMDP. 
