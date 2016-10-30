function Btree = BiTreePartition(A,cutoff)

if nargin < 2
    cutoff = 64;
end

N = size(A,1);

Btree = BiTreePartitionRecursion(1:N,[],cutoff);


    function Btree = BiTreePartitionRecursion(gidx,actidx,cutoff)
        
        if length(gidx) <= cutoff
            Btree.type = 'leaf';
            Btree.idx = gidx;
            Btree.actidx = FindNonZero(gidx,actidx);
            return;
        end
        
        [lidx,ridx,sepidx] = METIS_SepPartition(A(gidx,gidx));
        lidx = gidx(lidx);
        ridx = gidx(ridx);
        sepidx = gidx(sepidx);
        
        Btree.type = 'node';
        Btree.idx = sepidx;
        Btree.actidx = FindNonZero(gidx,actidx);
        
        subidx = [ FindNonZero(lidx,actidx) sepidx];
        Btree.ltree = BiTreePartitionRecursion(lidx,subidx,cutoff);
        
        subidx = [ FindNonZero(ridx,actidx) sepidx];
        Btree.rtree = BiTreePartitionRecursion(ridx,subidx,cutoff);
        
    end

    function subidx = FindNonZero(idx,actidx)
        
        [~,fullidx] = find(A(idx,actidx));
        subidx = actidx(unique(fullidx));
        
    end

end