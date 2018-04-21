%% [Results] = ESS(filename)
% This function segments the sclera region in an eye image. 
% Input: 
% -infilename: the name of the file storing the image.
% -outfilename: the name of the output file image. if outfilename is [], the function compute the name of the output file automatically.
% 
% Output: 
% -Results: a structure with the following fields
%    + NormImage (double) : the gray level image obtained by applying color correction. 
%    + Binary (logical)   : bynary image after gray level clustering.
%    + Candidates (double): foreground regions that are candidates to be selected as sclera or part of it.
%    + Sclera (double)    : the binary image composed by pixels selected as  belonging to the sclera regions.

function [Results] = ESS(infilename, outfilename)
if nargin<2, outfilename=[]; end;
% we set the output variable
Results = [];

% we try to load the input image
try
    A = imread(infilename);
catch
    sprintf('H_GLC (Error): file %s not found', infilename)
    Results = [];
    return;
end

% we check if the image format is compatible with next processing steps

% gray leve image
if(size(A,3)==3)
    sprintf('do nothing');
else
    sprintf('H_GLC (Error): the input must be a color image in RGB format')
    return;
end

% we resize the image to a fixed width/legnth of 500 pixels, according to which is the maximum dimension in the input image 
s = size(A);
s = s(1:2);
q=500/max(size(A));
A = imresize(A, q);

% we convert the input image in double format to ensure that arithmetic
% operations on pixels are performed by considering decimals.
A=double(A);

% we separate the three colo channels
R=A(:,:,1);
G=A(:,:,2);
B=A(:,:,3);

% we select the best parameter for the pixel mapping in the range [0,1]
mn = min([mean(R(:))+std2(R(:))/2,max(G(:))+std2(G(:))/2, max(B(:))+std2(B(:))/2]);

R = hsigmoidal(R, mn);
G = hsigmoidal(G, mn);
B = hsigmoidal(B, mn);

% we fuse the three color channels in a single gray scale image
Qsc = B+G-R; 

% we remap pixel values in the range [0,255]
QscIm = 255*(Qsc-min(Qsc(:)))/(max(Qsc(:))-min(Qsc(:)));

% we clusterize gray levels in order to separate the foreground from the
% background
[Mask, Q] = ClusterizeGrayLevels(QscIm, 0);
    
% we check if the outfilename is empty
nomeshort = outfilename;
if(isempty(nomeshort))
    [d_,n_,e_]=fileparts(infilename);
    nomeshort = [n_,'_res',e_];
    %{
    try
        slp = findstr(infilename, '\');
        slp = slp(length(slp));
        if(slp>0)
            nomeshort = infilename(slp+1:length(infilename));
        else
            nomeshort = infilename;
        end
    end
    %}
end

% we select which regions represents the sclera among the best candidates 
[C, Can] = SelectCompactComponents(Mask, Qsc, 0.03);

% we set the fields of the output structure
Results.NormImage = QscIm;
Results.Binary = Mask;
Results.Candidates = Can;
Results.Sclera = imfill(C, 'holes');

% we write the output image in tif format (non compressed)
Mresz =imresize(uint8(255*(Results.Sclera>0)), s);
nome = sprintf('.\\%s', nomeshort);
nome=strrep(nome, 'E', 'M');
nome=strrep(nome, 'jpg', 'tif');

imwrite(Mresz, nome);



        
%% Side Functions



% CLUSTERIZEGRAYLEVELS
function [Mask, Q] = ClusterizeGrayLevels(A, thbin)

th = 1.5;

B = double(A); %double(imresize(uint8(A), 2, 'bicubic'));

C=uint8(B);
% C = adapthisteq(uint8(B));
C = medfilt2(uint8(C), [8 8]);

% Calcoliamo la probabilità che ciascun pixel sia presente nell'immagine
h = histc(C(:), 0:255);
P = h/sum(h);

% clacoliamo la sparsità di ciascun tono di grigio
[SP] = ComputeSparseness(C);


T=zeros(size(B));
for i=1:size(B,1)
    for j=1:size(B,2)
        T(i,j)=SP(1+C(i,j));
    end
end

% estraiamo un vettore con soli i valori [0, 255], che realmente sono
% presenti nell'immagine
Pixels = double(unique(C(:)));

% contiamo quanti valori diversi sono presenti nell'immagine
n = length(Pixels);

delta = 1; 
% costruiamo una matrice delle distanze. La posizione (i,j) nella matrice
% indica quanto il pixel con tonalità i dista dal pixel con tonalità j
if(n > 1)
    M = zeros(n, n);
    for i=1:n
        for j=1:n
            M(i,j)=abs(P(Pixels(i)+1)*Pixels(i)-P(Pixels(j)+1)*Pixels(j))+delta*log(abs(Pixels(i)-Pixels(j))+1) + min(SP(i), SP(j))/(max(SP(i),SP(j))+1);
        end
    end
end

% clusterizziamo la matrice delle distanze fra toni di grigio
Clusters = zeros(1,n);
Assigned = zeros(1,n);
dist     = zeros(1,n);
nclusters = 0;
for i=1:n
    dist(:)=1000000;
    for j=1:nclusters
        k=1;
        while(Clusters(j,k)>0)
          dist(j) = min(dist(j), M(i,Clusters(j,k)));
          k = k+1;
        end
    end
    dminpos = find(dist==min(dist),1);
    dmin = min(dist);
    
    if(dmin < th)
        pos = find(Clusters(dminpos,:)==0,1);
        Clusters(dminpos, pos) = i;
        Assigned(i)=1;
    else
        nclusters = nclusters + 1;
        Clusters(nclusters, 1) = i;
        Assigned(i) = 1;
    end
end

Q = zeros(size(B));
for i=1:n
    c = (C==(Pixels(i)));
    [x,y] = find(Clusters==i);
    Q = Q + c*x;
end

S=Q(:);
if(thbin==0)
    thbin = mean(S(:))+0.75*std(S(:));
end


THfinal = 0.4;

Mask = Q>thbin; 




% COMPUTESPARSENESS
function [SP] = ComputeSparseness(A)

SP = zeros(1,256);

for i=0:255
    [x, y] = find(A==i);
    
    Xc = mean(x);
    Yc = mean(y);
    
    t  = sqrt(mean((x-Xc).^2 + (y-Yc).^2));
    if(isnan(t)==0)
        SP(i+1) = t;
    end
end


function val = hsigmoidal(x, xmax)

a = 2 + sqrt(3);
b = 7 - 4*sqrt(3);

val = (1 - b.^(x./xmax))./(a*b.^(x./xmax) + 1);


function [C, Can] = SelectCompactComponents(A, Ag, prc)

M = sum(sum(A));

[L, n] = bwlabel(A>0);

Can = zeros(size(A));
C = zeros(size(A));
d = zeros(1,n);
cm = zeros(1,n);

diago = sqrt(size(A,1)^2+size(A,2)^2);
for i=1:n
    c = (L==i);
    
    if(sum(c(:))>prc*M)
        [X,Y] = find(c>0);
        xm = mean(X);
        ym = mean(Y);
        dst = 0.5*sqrt((xm-size(A,1)/2)^2+(ym-size(A,2)/2)^2)/diago;
        cont = bwmorph(c, 'remove');
        cont = Ag.*cont;
        cont = mean(cont(cont>0));
        cf = imfill(c, 'holes').*Ag/(cont+dst);
        a = bwconvhull(c);
        cm(i) = sum(c(:))/sum(a(:));
        d(i) = cm(i)*mean(cf(cf>0));
        Can = Can + c;
    end
end

[ds, I] = sort(d, 'descend');

try
    C = C + (L==I(1));
    
    [X,Y] = find(C==1);
    xmC = round(mean(X));
    i=2;
    while(ds(i) >0)
        c = (L==I(i));
        [X,Y] = find(c==1);
        xm = min(abs(X-xmC));

        if(xm==0 && cm(I(i))>0.7*cm(I(1)))
            C = C + (L==I(i));
        end
        i=i+1;
    end
end

function convex_hull = bwconvhull(varargin)
%BWCONVhull Generate convex hull image from binary image.
%   CH = BWCONVHULL(BW) computes the convex hull of all objects in BW and
%   returns CH, binary convex hull image.  BW is a logical 2D image and CH
%   is a logical convex hull image, containing the binary mask of the
%   convex hull of all foreground objects in BW.
%
%   CH = BWCONVHULL(BW,METHOD) specifies the desired method for computing
%   the convex hull image.  METHOD is a string and may have the following
%   values:
%
%      'union'   : Compute convex hull of all foreground objects, treating
%                  them as a single object.  This is the default method.
%      'objects' : Compute the convex hull of each connected component of
%                  BW individually.  CH will contain the convex hulls of
%                  each connected component.
%
%   CH = BWCONVHULL(BW,'objects',CONN) specifies the desired connectivity
%   used when defining individual foreground objects.  The CONN parameter
%   is only valid when the METHOD is 'objects'.  CONN may have the
%   following scalar values:
%
%      4 : two-dimensional four-connected neighborhood
%      8 : two-dimensional eight-connected neighborhood {default}
%
%   Additionally, CONN may be defined in a more general way, using a 3-by-3
%   matrix of 0s and 1s.  The 1-valued elements define neighborhood
%   locations relative to the center element of CONN.  CONN must be
%   symmetric about its center element.
%
%   Example
%   -------
%   subplot(2,2,1);
%   I = imread('coins.png');
%   imshow(I);
%   title('Original');
%   
%   subplot(2,2,2);
%   BW = I > 100;
%   imshow(BW);
%   title('Binary');
%   
%   subplot(2,2,3);
%   CH = bwconvhull(BW);
%   imshow(CH);
%   title('Union Convex Hull');
%   
%   subplot(2,2,4);
%   CH_objects = bwconvhull(BW,'objects');
%   imshow(CH_objects);
%   title('Objects Convex Hull');
%
%   See also BWCONNCOMP, BWLABEL, LABELMATRIX, REGIONPROPS.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2011/08/09 17:49:02 $

[BW, method, conn] = parseInputs(varargin{:});

% Label the image
if strcmpi(method,'union')
    % 'union' : label all 'true' pixels as a single region
    labeled_image = uint8(BW);
else
    % 'objects' : label as normal
    labeled_image = bwconncomp(BW,conn);
end

% Call regionprops
blob_props = regionprops(labeled_image,'BoundingBox','ConvexImage');
num_blobs = length(blob_props);
[rows columns] = size(BW);

% Loop over all blobs getting the CH for each blob one at a time, then add
% it to the cumulative CH image.
convex_hull = false(rows, columns);
for i = 1 : num_blobs
    m = blob_props(i).BoundingBox(4);
    n = blob_props(i).BoundingBox(3);
    r1 = blob_props(i).BoundingBox(2) + 0.5;
    c1 = blob_props(i).BoundingBox(1) + 0.5;
    rows = (1:m) + r1 - 1;
    cols = (1:n) + c1 - 1;
    convex_hull(rows,cols) = convex_hull(rows,cols) | blob_props(i).ConvexImage;
end


%------------------------------------------------
function [BW,method,conn] = parseInputs(varargin)

%narginchk(1,3);

BW = varargin{1};
validateattributes(BW, {'logical' 'numeric'}, {'2d', 'real', 'nonsparse'}, ...
    mfilename, 'BW', 1);

if ~islogical(BW)
    BW = BW ~= 0;
end

if nargin == 1
    % BWCONVHULL(BW)
    method = 'union';
    conn = 8;
    
elseif nargin == 2
    % BWCONVHULL(BW,METHOD)
    method = varargin{2};
    conn = 8;
    
else
    % BWCONVHULL(BW,METHOD,CONN)
    method = varargin{2};
    conn = varargin{3};
    
    % special case so that we go through the 2D code path for 4 or 8
    % connectivity
    if isequal(conn, [0 1 0;1 1 1;0 1 0])
        conn = 4;
    end
    if isequal(conn, ones(3))
        conn = 8;
    end
    
end

% validate inputs (accepts partial string matches)
method = validatestring(method,{'union','objects'},mfilename,'METHOD',2);

% validate connectivity
is_valid_scalar = isscalar(conn) && (conn == 4 || conn == 8);
if is_valid_scalar
    return
end

% else, validate 3x3 connectivity matrix

% 3x3 matrix...
is_valid_matrix = isnumeric(conn) && isequal(size(conn),[3 3]);
% with all 1's and 0's...
is_valid_matrix = is_valid_matrix && all((conn(:) == 1) | (conn(:) == 0));
% whos center value is non-zero
is_valid_matrix = is_valid_matrix && conn((end+1)/2) ~= 0;
% and which is symmetrix
is_valid_matrix = is_valid_matrix && isequal(conn(1:end), conn(end:-1:1));

if ~is_valid_matrix
    error(message('images:bwconvhull:invalidConnectivity'))
end