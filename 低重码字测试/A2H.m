%A(i,j)=0 for zero-entry (z*z zero-matrix), in [1,z] for non-zero entry (shifted version of z*z identity matrix)
%H is sparse logic matrix
function H = A2H(A,z)

  [mrow, ncol] = size(A);

  A = A - 1;%?
  nzmax = sum(sum(A~=-1))*z; %total degree or weight (number of non-zero entries)
  irow(1:nzmax) = 0; jcol(1:nzmax) = 0; %coordinates (x,y) for each non-zero entry in H

  k=1:z;
  for i1=1:mrow, for j1=1:ncol
    if A(i1,j1)~=-1
      irow(k) = i1*z-z+1 + (0:z-1);
      jcol(k) = j1*z-z+1 + mod(A(i1,j1)+(0:z-1),z);
      k = k+z;
    end
  end, end
 H = logical( sparse(irow',jcol',ones(1,nzmax),mrow*z,ncol*z,nzmax) );

return