% This function takes vessel locations, details of each cycle and
% an estimated spin velocity, v (cm/s), and determines the encoding matrix A.
%
% Vessel locations is an Nx2 array [LRPos APPos] R and P are negative and
% distances are in mm.
%
% Cycle details is a Mx3 array (TagAngle, TagMode, vA, vB; ...), where
% the tagging angle is specified in degrees (0 is AP, 90 is LR), TagMode
% is an integer (0 is tag all, 1 control all, 2 tag A, 3 tag B), and vA
% and vB specify the tag/control positions in mm from isocentre.
% PerfectInv determines whether perfect inversion is assumed (default true)
% and Laminar determines whether laminar rather than plug flow is assumed
% (default true).
%
% function A = Enc_Mtx(VesLoc, CycDets, v, PerfectInv, Laminar)

function A = hw_Enc_Mtx_ProjCoW(VesLoc, CycDets, v, PerfectInv, Laminar, Unipolar,BipolarIEPath)
  
  % Deal with optional arguments
  if nargin < 3; v = 30; end
  if nargin < 4; PerfectInv = true; end
  if nargin < 5; Laminar = true; end
  if nargin < 6; Unipolar = false; end
  
  % Convert velocity to m/s
  v = v / 100;
  
  % Load the data describing the inversion profile
  if Laminar == true && Unipolar == false
	  IE = load(BipolarIEPath).hw_PcaslParas_vs_pos_and_v_gaussian_dsv;
  elseif Laminar == false && Unipolar == false
	  IE = load(['///']);
  elseif Laminar == true && Unipolar == true % Unipolar
      IE = load(['///']);
  else
      error('No unipolar non-laminar simulation available');
  end
  
  % Record some parameters used below and convert m-> mm
  Pa = IE.Pa(1)*1000;  Pb = IE.Pb(1)*1000;  XPrime = IE.Ps(1,:)*1000; 
  MinXPrime = min(XPrime);  MaxXPrime = max(XPrime);

  % Determine the tagging and control efficiency at the exact tag and
  % control locations
  % NB. At exactly on the vessel, tag is identical for selective or
  % non-selective cycles (in simulations at least)
  TagPosIdx  = find(IE.Ps(1,:) == IE.Pa(1));
  CntlPosIdx = find(IE.Ps(1,:) == IE.Pb(1));
  IETag  = interp1(IE.v, IE.IE(TagPosIdx, :), v);
  IECntl = interp1(IE.v, IE.IE(CntlPosIdx,:), v);
  
  % Determine the encoding matrix for each cycle
  M = size(CycDets,1); N = size(VesLoc,1);
  A = zeros(M, N);
  
  for CycNo = 1:M
  
    switch CycDets(CycNo,2) % Tag mode
     case 0 % Tag All
      A(CycNo, :) = IETag;
      SelTag = false;
      
     case 1 % Control all
      A(CycNo, :) = IECntl;
      SelTag = false;
      
     case 2 % For selective tagging, set the offset used to
            % calculate the vessel position in IE space
      OffSet = 0;  
      SelTag = true;
      
     case 3
      OffSet = Pb - Pa;
      SelTag = true;
      
     otherwise
      error(['Tag mode not recognised for cycle ' num2str(CycNo)])
      
    end

    % Work out matrix entries for selective tagging
    if SelTag == true
      
      % Extract vessel locations in the relevant direction
      
      % Convert degrees to radians
      theta = todeg2rad(CycDets(CycNo,1));
      
      % Determine a unit vector in the direction of the gradient
      G = [ sin(theta), cos(theta) ];
      
      % The distance from isocentre along the gradient direction is given
      % by the dot product of the gradient unit vector with the vessel location
      TempVesLocs = sum( VesLoc .* repmat( G, [size(VesLoc,1) 1] ) , 2);

      % Convert vessel locations into the space in which they were
      % simulated
      vA = CycDets(CycNo, 3);  vB = CycDets(CycNo, 4);
      
      TempVesLocsPrime = Pa + OffSet + ...
          (TempVesLocs - vA) / (vB - vA) * (Pb - Pa);
      
      % Wrap locations into the simulated space since the profile is
      % periodic
      TempVesLocsPrime = mod(TempVesLocsPrime-MinXPrime, MaxXPrime-MinXPrime) ...
          + MinXPrime;
      
      % Interpolate the simulated data to fill in the encoding matrix
      InvEffs = interp2(IE.v, XPrime, IE.IE, v, TempVesLocsPrime);
      A(CycNo,:) = InvEffs;
      
    end
    
    
  end
  
  % If assuming perfect inversion and control efficiency, correct for
  % this here
  if PerfectInv == true
    A = (A - IECntl) / (IETag - IECntl);
  end
  
  % Convert from IE to Mz
  A = -2 * A + 1;
  
  % Add a column for static tissue
  A(:,N+1) = 1;
  
       
