% Returns the pathname for a file selected by the user with the
% filter criteria set to anything (rather than the default matlab file formats)
%
% function out = togetfile(Comment);

function out = togetfile(Comment);
  
  if nargin == 0
    Comment = '';
  end
  
  [FileName, PathName] = uigetfile('*.*',Comment);
  
  out = [PathName '/' FileName];