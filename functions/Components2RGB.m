% This function takes an image with Nvc components and converts it into
% an RGB map, with each component taking a different colour.  The input
% array must have its fourth dimension as vascular components.
%
% function out = Components2RGB(in)

function [out, Col] = Components2RGB(in,UseDCs,PermOrder,CustomCols)
  
  if nargin < 2; UseDCs = false; end
  if nargin < 3; PermOrder = 1:size(in,4); end
  if nargin < 4; CustomCols = []; end
  
  % Find the number of components
  s = size(in); s((ndims(in)+1):5) = 1;
  Nvc = s(4);
  Ncoils = s(5);
  
  % Decompose into a magnitude and vector image
  Mag = sqrt(sum(abs(in).^2,4));
  Vec = in./repmat(Mag,[1 1 1 Nvc]);

  % Permute so vascular components are the last dimension
  Vec = permute(Vec,[1 2 3 5 4]);

  % Reshape into a matrix form: voxels by components
  Vec_mat = reshape(Vec,[],Nvc);

  % Define the colours to transform components with
  if isempty(CustomCols)
    Col = VTColours(Nvc,UseDCs);
  else
    Col = CustomCols;  
    Col = Col./repmat(sqrt(sum(Col.^2,2)),[1 3]); % Normalise
  end
  
  % Permute
  Col = Col(PermOrder,:);
  
  % Multiply the matrices to transform components to colour in each
  % voxel. NB. out_mat is voxels by RGB.
  out_mat = Vec_mat * Col;
  
  % Normalise these vectors
  out_mat = out_mat ./ repmat(sqrt(sum(out_mat.^2,2)),[1 3]);
  
  % Transform back to the original dimensions, with RGB as the last
  % dimension. 
  out = reshape(out_mat, [s(1:3) Ncoils 3]);
  
  % Permute so RGB values are the fourth dimension.
  out = permute(out,[1 2 3 5 4]);
  
  % Multiple the output by the magnitude image
  out = out.* repmat(Mag,[1 1 1 size(out,4)]);

  % Zero values where the magnitude is zero
  out(repmat(Mag,[1 1 1 size(out,4)]) == 0) = 0;

  % Rounding errors can cause RGB values to exceed unity, so correct for
  % this here
  out(out > 1) = 1;