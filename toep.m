function mat=toep(column,row)

%	TOEP Toeplitz matrix
%
%	mat=toep(column,row)
%
%	The length(column) x length(row) element toeplitz matrix is generated from the
%	column and row vectors.  If unspecified, the column vector is entered after the
%	prompt from the keyboard, and the default assignment row=column is used.  For
%	column and row matrices, the initial column and row vectors, respectively, are
%	used.  For unequal column and row vector initial elements, the column vector
%	initial element is used for the toeplitz matrix.

%	Implemented using MATLAB 5.3.1
%
%	Example:
%
%	» m=toep((1:3)'*(1+i),[1+i (4:6)*(1-i)])
%
%	m =
%
%	   1.0000 + 1.0000i   4.0000 - 4.0000i   5.0000 - 5.0000i   6.0000 - 6.0000i
%	   2.0000 + 2.0000i   1.0000 + 1.0000i   4.0000 - 4.0000i   5.0000 - 5.0000i
%	   3.0000 + 3.0000i   2.0000 + 2.0000i   1.0000 + 1.0000i   4.0000 - 4.0000i
%
%	Copyright (c) 2000
%	Tom McMurray
%	mcmurray@teamcmi.com

%	assign default input column vector

if ~nargin
   column=input('enter column vector or return for 0 output\n');
   if isempty(column)
      mat=0;
      return
   end
end

%	while column is unsupported, enter supported column or return for 0 output

while isempty(column)|~isnumeric(column)|ndims(column)>2
   column=input(['column vector is empty, nonnumeric, or > 2 dimensional:\nenter '...
         'column vector or return for 0 output\n']);
   if isempty(column)
      mat=0;
      return
   end
end

%	if column is a matrix, modify to the initial column vector

[sizecol1,sizecol2]=size(column);
if min(sizecol1,sizecol2)~=1
   column=column(:,1);
   warning('using initial column vector of column matrix')
   
%	else if column is a row vector, modify to a column vector
   
elseif sizecol2~=1
   column=column(:);
   sizecol1=sizecol2;
end

%	assign default input row vector

if nargin<2
   row=column;
end

%	while row is unsupported, enter supported row or return for row=column

while isempty(row)|~isnumeric(row)|ndims(row)>2
   row=input(['row vector is empty, nonnumeric, or > 2 dimensional:\nenter row '...
         'vector or return for row=column\n']);
   if isempty(row)
      row=column;
   end
end

%	if row is a matrix, modify the initial row vector to a column vector

[sizerow1,sizerow2]=size(row);
if min(sizerow1,sizerow2)~=1
   row=row(1,:).';
   warning('using initial row vector of row matrix')
   
%	else if row is a column vector, sizerow2=sizerow1
   
elseif sizerow1~=1
   sizerow2=sizerow1;
   
%	else row is a row vector, modify to a column vector
   
else
   row=row(:);
end

%	if column(1)~=row(1), use column vector initial element

if column(1)~=row(1)
   warning(['modifying row initial element = ' num2str(row(1))...
         ' to column initial element = ' num2str(column(1))])
end

%	generate column vector from input vectors

vec=[row(sizerow2:-1:2);column];
lenvec=sizecol1+sizerow2-1;

%	generate toeplitz matrix

mat=zeros(sizecol1,sizerow2);
for k=1:sizerow2
   k1=k-1;
   mat(:,k)=vec(sizerow2-k1:lenvec-k1);
end