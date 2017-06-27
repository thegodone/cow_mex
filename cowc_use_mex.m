function [Warping,XWarped,Diagnos,Table] = cowc_use_mex(T,X,Seg,Slack)
% function [Warping,XWarped,Diagnos] = cow(T,X,Seg,Slack,Options);
% Correlation Optimized Warping function with linear interpolation
% Giorgio Tomasi / Frans van den Berg 070821 (GT)
%
% Thomas Skov 061219 - line 277 changed to work in MATLAB version 6.5
%
% speed improvment:
%
% Guillaume Godin  27-05-2017 - line 172 call nakeinterp1 
% Guillaume Godin / Addisalem SHIFERAW 27-05-2017 - line 146 call mex
% version sub routine to improve the speed
%
% in:  T (1 x nt) target vector
%      X (mP x nP) matrix with data for mP row vectors of length nP to be warped/corrected
%      Seg (1 x 1) segment length; number of segments N = floor(nP/m)
%       or (2 x N+1) matrix with segment (pre-determined) boundary-points
%                    first row = index in "xt", must start with 1 and end with "nt"
%                    second row = index in "xP", must start with 1 and end with "nP"
%      Slack (1 x 1) 'slack' - maximum range or degree of warping in segment length "m"
%      Options (1 x 5) 1 : triggers plot and progress-text (note: only last row/object in "xP" is plotted)
%                      2 : correlation power (minimum 1th power, maximum is 4th power)
%                      3 : force equal segment lengths in "xt" and "xP" instead of filling up "xt" with N boundary-points
%                          (notice that different number of boundaries in "xt" and "xP" will generate an error)
%                      4 : fix maximum correction to + or - options(4) points from the diagonal
%                      5 : save in "diagnos" the table with the optimal values of loss function and predecessor (memory
%                          consuming for large problems - on how to read the tables are in the m-file
%              default [0 1 0 0 0] (no plot; power 1; no forced equal segment lengths; no band constraints; no Table in "diagnos")
%
% out: Warping (mP x N x 2) interpolation segment starting points (in "nP"
%          units) after warping (first slab) and before warping (second slab)
%          (difference of the two = alignment by repositioning segment
%          boundaries; useful for comparing correction in different/new objects/samples)
%      XWarped (mP x nt) corrected vectors (from "xP" warped to mach "xt")
%      Diagnos (struct) warping diagnostics: options, segment, slack,
%          index in target ("xt", "warping" is shift compared to this) and sample ("xP"), search range in "xP", computation time
%          (note: diagnostics are only saved for one - the last - signal in "xP")
%
% based on: Niels-Peter Vest Nielsen, Jens Micheal Carstensen and Jørn Smedegaard 'Aligning of singel and multiple
%           wavelength chromatographic profiles for chemometric data analysis using correlation optimised warping'
%           J. Chrom. A 805(1998)17-35
%
% Reference: Correlation optimized warping and dynamic time warping as preprocessing methods for chromatographic Data
%            Giorgio Tomasi, Frans van den Berg and Claus Andersson, Journal of Chemometrics 18(2004)231-241
%
% Authors:
% Giorgio Tomasi / Frans van den Berg
% Royal Agricultural and Veterinary University - Department of Food Science
% Quality and Technology - Spectroscopy and Chemometrics group - Denmark
% email: gt@kvl.dk / fb@kvl.dk - www.models.kvl.dk
% this is the basic code without options setting as default for c++ code
% conversion
Options = [0 1 0 0 0];

%% Initialise
[nX,pX] = size(X);         % nX     : number of signals that are to be aligned
% pX     : number of data points in each signal
pT      = size(T,2);       % pT     : number of data points in the target
XWarped = zeros(nX,pT);    % XWarped: initialise matrix of warped signals
Time    = zeros(1,1);      % Time   : processing time
%% Initialise segments
Seg        = round(Seg);      % Only integers are currently allowed as segment boundaries

if Seg > min(pX,pT)
    error('Segment length is larger than length of the signal');
end

nSeg               = floor((pT - 1) / (Seg - 1));
LenSeg=zeros(2,nSeg);
LenSeg(1:2,1:nSeg) = Seg - 1;
if floor((pX - 1) / (Seg - 1)) ~= nSeg
    error('For non-fixed segment lengths the target and the signal do not have the same number of segments (try Options(3))');
end

temp = rem(pT - 1,LenSeg(1,1)); % The remainders are attached to the last segment in the target and in the reference
if (temp > 0)
    LenSeg(1,nSeg) = LenSeg(1,nSeg) + temp;
end

temp = rem(pX - 1,LenSeg(2,1));

if temp > 0
    LenSeg(2,nSeg) = LenSeg(2,nSeg) + temp;
end


if any(LenSeg(:) <= Slack + 2) % Two points are the minimum required for linear interpolation
    error('The slack cannot be larger than the length of the segments');
end

bT      = cumsum([1,LenSeg(1,:)]);
bP      = cumsum([1,LenSeg(2,:)]);
Warping = zeros(nX,nSeg + 1);

%% Check slack
if length(Slack) > 1 % Different slacks for the segment boundaries will be implemented
    if size(Slack,2) <= nSeg
        error('The number of slack parameters is not equal to the number of optimised segments');
    end
    fprintf('\n Multiple slacks have not been implemented yet')
    return
end
Slacks_vec = zeros(2*Slack+1,1);
Slacks_vec(1:2*Slack+1,1) = -Slack:Slack;                     % All possible slacks for a segment boundary
%% Set feasible points for boundaries
Bounds      = ones(2,nSeg + 1);
% Slope Constraints
offs        = (Slack * [-1,1]') * (0:nSeg);
Bounds_a    = bP(ones(2,1),1:nSeg + 1) + offs;
Bounds_b    = bP(ones(2,1),1:nSeg + 1) + offs(:,nSeg + 1:-1:1);
Bounds(1,:) = max(Bounds_a(1,:),Bounds_b(1,:));
Bounds(2,:) = min(Bounds_a(2,:),Bounds_b(2,:));


%% Calculate first derivatives for interpolation
% Xdiff = diff(X,1,2);
% 
% %% Calculate coefficients and indexes for interpolation
% Int_Coeff = cell(nSeg,1);
% Int_Index = Int_Coeff;
% [A,B]                             = InterpCoeff(LenSeg(1,1) + 1,LenSeg(2,1) + Slacks_vec + 1,Slacks_vec);
% [Int_Coeff{1:nSeg - 1}]           = deal(A);
% [Int_Index{1:nSeg - 1}]           = deal(B);
% [A_last, B_last] =InterpCoeff(LenSeg(1,nSeg) + 1,LenSeg(2,nSeg) + Slacks_vec + 1,Slacks_vec);
% Int_Coeff{nSeg} = A_last;
% Int_Index{nSeg} = B_last;
% [Int_Coeff{nSeg},Int_Index{nSeg}] = InterpCoeff(LenSeg(1,nSeg) + 1,LenSeg(2,nSeg) + Slacks_vec + 1,Slacks_vec);
% 

%% Dynamic Programming Section
Table_Index    = cumsum([0,diff(Bounds) + 1]);       % Indexes for the first node (boundary point) of each segment in Table
Table          = zeros(3,Table_Index(nSeg + 2),nX);  % Table: each column refer to a node
%        (1,i) position of the boundary point in the signal
%        (2,i) optimal
%        value of the loss function up to node (i)
%        (3,i) pointer to optimal preceding node (in Table)
Table(2,2:end,1:nX) = -Inf;                          % All loss function values apart from node (1) are set to -Inf
for i_seg = 1:nSeg + 1                               % Initialise Table
    v       = (Bounds(1,i_seg):Bounds(2,i_seg))';
    Table(1,Table_Index(i_seg) + 1:Table_Index(i_seg + 1),:) = v(:,ones(nX,1));
end

%INPUT
% Slacks_vec - column array 
% LenSeg - column array (2 dim)
% Table_Index - row array
%
% Forward phase

%[ind_Array, pos_value_array] =cowc_mex(nSeg,Slacks_vec, LenSeg, Table_Index, A,B,A_last, B_last, T, X, Table,Bounds, bT, Xdiff);
[ind_Array, pos_value_array] =cowc_mex(T,X, Seg,Slack);
Table(2,:,:) = ind_Array;
Table(3,:,:) = pos_value_array;

for i_sam = 1:nX                                  % Loop over samples/signals
    % Backward phase
    Pointer                 = size(Table,2);           % Backtrace optimal boundaries using the pointers in Table
    Warping(i_sam,nSeg + 1) = pX;
    for i_bound = nSeg:-1:1
        Pointer                =Table(3,Pointer,i_sam);
        Warping(i_sam,i_bound) = Table(1,Pointer,i_sam);
    end
    
end
Warping(:,:,2) = bT(ones(nX,1),:);

%% Output
if (nargout > 1)  % Reconstruct aligned signals
    for i_seg = 1:nSeg
        indT = bT(i_seg):bT(i_seg + 1);
        lenT = bT(i_seg + 1) - bT(i_seg);
        for i_sam = 1:nX
            indX                = Warping(i_sam,i_seg):Warping(i_sam,i_seg + 1);
            lenX                = Warping(i_sam,i_seg + 1) - Warping(i_sam,i_seg);
            % NB the right handside expression must be transposed to fit MATLAB version 6.5
            XWarped(i_sam,indT) = nakeinterp1(indX' - Warping(i_sam,i_seg) + 1,X(i_sam,indX)',(0:lenT)'/lenT * lenX + 1)';
        end
    end
end
if (nargout > 2)    % Save some diagnostics if requested
    Diagnos = struct('indexP',bP,'indexT',bT,'Nsegments',nSeg,'options',Options,'rangeP',Bounds',...
        'segment_length',LenSeg,'slack',Slack,'table',[],'time',Time);
    
end