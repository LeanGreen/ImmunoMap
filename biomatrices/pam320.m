function [matrix, matrixInfo] = pam320
%PAM320 substitution matrix in 1/4 bit units,
%   Expected score = -0.741, Entropy = 0.224 bits
%   Lowest score = -8, Highest score = 22
%
%   Order:
%   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
%
%   [MATRIX,MATRIXINFO] = PAM320 returns a structure of information about
%   the matrix with fields Name, Scale, Entropy, ExpectedScore, LowestScore,
%   HighestScore and Order.

% Source:  ftp://ftp.ncbi.nih.gov/blast/matrices/


matrix = [...
   1 -1  0  1 -2  0  0  1 -1  0 -2 -1 -1 -4  1  1  1 -6 -4  0  0  0  0 -8;...
  -1  6  0 -1 -4  1 -1 -2  2 -2 -3  4  0 -5  0  0 -1  3 -4 -2  0  0 -1 -8;...
   0  0  2  2 -4  1  2  1  2 -2 -3  1 -2 -4  0  1  0 -4 -2 -2  2  1  0 -8;...
   1 -1  2  4 -5  2  3  1  1 -2 -4  0 -2 -6  0  0  0 -7 -5 -2  3  3 -1 -8;...
  -2 -4 -4 -5 15 -5 -5 -3 -4 -2 -6 -6 -5 -4 -3  0 -2 -8  1 -2 -4 -5 -3 -8;...
   0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0  0  0 -5 -4 -2  2  3  0 -8;...
   0 -1  2  3 -5  2  4  1  1 -2 -3  0 -2 -5  0  0  0 -7 -5 -2  3  3  0 -8;...
   1 -2  1  1 -3 -1  1  5 -2 -2 -4 -1 -3 -5  0  1  0 -7 -5 -1  1  0 -1 -8;...
  -1  2  2  1 -4  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2  1  2  0 -8;...
   0 -2 -2 -2 -2 -2 -2 -2 -2  4  3 -2  2  1 -2 -1  0 -5 -1  4 -2 -2  0 -8;...
  -2 -3 -3 -4 -6 -2 -3 -4 -2  3  7 -3  4  3 -2 -3 -2 -2  0  2 -3 -3 -1 -8;...
  -1  4  1  0 -6  1  0 -1  0 -2 -3  5  0 -5 -1  0  0 -3 -5 -2  1  1 -1 -8;...
  -1  0 -2 -2 -5 -1 -2 -3 -2  2  4  0  6  1 -2 -1 -1 -4 -2  2 -2 -2  0 -8;...
  -4 -5 -4 -6 -4 -5 -5 -5 -2  1  3 -5  1 11 -5 -3 -3  1  9 -1 -5 -5 -2 -8;...
   1  0  0  0 -3  0  0  0  0 -2 -2 -1 -2 -5  6  1  1 -6 -5 -1  0  0  0 -8;...
   1  0  1  0  0  0  0  1 -1 -1 -3  0 -1 -3  1  1  1 -3 -3 -1  1  0  0 -8;...
   1 -1  0  0 -2  0  0  0 -1  0 -2  0 -1 -3  1  1  2 -5 -3  0  0  0  0 -8;...
  -6  3 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  1 -6 -3 -5 22  1 -6 -6 -6 -4 -8;...
  -4 -4 -2 -5  1 -4 -5 -5  0 -1  0 -5 -2  9 -5 -3 -3  1 12 -2 -3 -4 -2 -8;...
   0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4 -2 -2  0 -8;...
   0  0  2  3 -4  2  3  1  1 -2 -3  1 -2 -5  0  1  0 -6 -3 -2  2  2  0 -8;...
   0  0  1  3 -5  3  3  0  2 -2 -3  1 -2 -5  0  0  0 -6 -4 -2  2  3  0 -8;...
   0 -1  0 -1 -3  0  0 -1  0  0 -1 -1  0 -2  0  0  0 -4 -2  0  0  0 -1 -8;...
  -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  1;...
    ];

if nargout >1
    matrixInfo.Name = 'PAM320';
    matrixInfo.Scale = 1/4;
    matrixInfo.Entropy = 0.224 ;
    matrixInfo.ExpectedScore = -0.741;
    matrixInfo.LowestScore = -8;
    matrixInfo.HighestScore = 22;
    matrixInfo.Order = 'ARNDCQEGHILKMFPSTWYVBZX*';
end

