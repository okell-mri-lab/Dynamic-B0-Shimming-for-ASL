% Returns colours for mapping vascular components (N x 3 matrix)
%
% cols = VTColours(Nvc)
function cols = VTColours(Nvc,UseDistinguishableColours)

if nargin < 2; UseDistinguishableColours = false; end

  if Nvc <=3
    cols = eye(3);
  elseif (Nvc <= 13) && ~(UseDistinguishableColours)
     %        R    G    B
    cols = [  1    0    0 ; ... % Red
              0    1    0 ; ... % Green
              0  0.5    1 ; ... % Light blue
              1    0  0.5 ; ... % Magenta
              1    1    0 ; ... % Yellow
              0    1  0.5 ; ... % Turquoise
              1   0.5   0 ; ... % Orange
              0.5  0    1 ; ... % Purple
              1    1    1 ; ... % Grey
              1   0.7 0.2 ; ... % ??
              1   0.2 0.7 ; ... % ??            
              0.7  1  0.2 ; ... % ??
              0.7 0.2   1 ; ... % ??              
              ];   
  else
    cols = distinguishable_colors(Nvc,[0 0 0]);
  end
  
  % Ensure all colours have equal brightness
  for r = 1:size(cols,1); 
    cols(r,:) = cols(r,:)/sqrt(sum(cols(r,:).^2));
  end
  
  % Only return the requested colsours
  cols = cols(1:Nvc,:);