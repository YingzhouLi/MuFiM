function Btree = BiTreePartition(A,cutoff)
% BITREEPARTITION Generates binary tree partition
%   Btree = BITREEPARTITION(A) generates the binary tree partition for the
%   graph of sparse matrix A. The graph of A is recursively bipartitioned
%   via Metis. The default cutoff for the leaf level is 64 points.
%
%   Btree = BITREEPARTITION(A,cutoff) generates the binary tree partition
%   for the graph of sparse matrix A. The leaf node contains less than
%   cutoff points.
%
%   See also MULTIFRONTAL, SYMBOLMF.

%   Copyright 2016 Yingzhou Li, Stanford University

if nargin < 2
    cutoff = 128;
end

N = size(A,1);

Btree = BiTreePartitionRecursion(1:N,[],cutoff);

    function Btree = BiTreePartitionRecursion(gidx,actidx,cutoff)
        
        if length(gidx) <= cutoff
            Btree.type = 'leaf';
            Btree.idx = sort(gidx);
            Btree.actidx = sort(FindNonZero(gidx,actidx));
            return;
        end
        
        [lidx,ridx,sepidx] = METIS_SepPartition(A(gidx,gidx));
        lidx = gidx(lidx);
        ridx = gidx(ridx);
        sepidx = gidx(sepidx);
        
        Btree.type = 'node';
        Btree.idx = sort(sepidx);
        Btree.actidx = sort(FindNonZero(gidx,actidx));
        
        subidx = [ FindNonZero(lidx,actidx) sepidx];
        Btree.ltree = BiTreePartitionRecursion(lidx,subidx,cutoff);
        
        subidx = [ FindNonZero(ridx,actidx) sepidx];
        Btree.rtree = BiTreePartitionRecursion(ridx,subidx,cutoff);
        
    end

    function subidx = FindNonZero(idx,actidx)
        
        [I,J] = find(A(:,actidx));
        iI = ismembc(I,idx);
        fullidx = J(iI);
        subidx = actidx(unique(fullidx));
        
    end

end