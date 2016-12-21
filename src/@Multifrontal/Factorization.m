function MF = Factorization(MF,A)
% FACTORIZATION Multifrontal factorization
%   MF = FACTORIZATION(MF,A) factorizes the matrix A as a data-sparse
%   multiplication of lower or upper trangular matrices. If A is a
%   numerically symmetric matrix, A is factorized as
%       L_1 L_2 ... L_k D L_k^T ... L_2^T L_1^T,
%   where L_i is lower trangular matrix. If A is a pattern symmetric
%   matrix, A is factorized as
%       A = L_1 L_2 ... L_k U_k ... U_2 U_1,
%   where L_i is lower trangular matrix and U_j is upper trangular matrix.
%
%   See also SYMBOLMF, MULTIFRONTAL.

%   Copyright 2016 Yingzhou Li, Stanford University

if MF.symm == 1
    [MF.Ltree,MF.Dtree] = ...
        SymmFactorizationRecursion(MF.symboltree);
elseif MF.symm == 2
    [MF.Ltree,MF.Utree] = ...
        PatSymmFactorizationRecursion(MF.symboltree);
end

%====================================================================
    function [Ltree,Dtree] = SymmFactorizationRecursion(Stree)
        
        if strcmpi(Stree.type,'node')
            [Ltree.ltree,Dtree.ltree] = ...
                SymmFactorizationRecursion(Stree.ltree);
            
            [Ltree.rtree,Dtree.rtree] = ...
                SymmFactorizationRecursion(Stree.rtree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        [L,D] = ldl(full(A(idx,idx)));
        
        ALDinv = A(actidx,idx)/L'/D;
        
        A(actidx,actidx) = A(actidx,actidx) ...
            - sparse(A(actidx,idx)/A(idx,idx)*A(actidx,idx)');
        
        Ltree.Mat = L;
        Ltree.AMatinv = ALDinv;
        Dtree.Mat = D;
        
    end

%====================================================================
    function [Ltree,Utree] = PatSymmFactorizationRecursion(Stree)
        
        if strcmpi(Stree.type,'node')
            [Ltree.ltree,Utree.ltree] = ...
                PatSymmFactorizationRecursion(Stree.ltree);
            
            [Ltree.rtree,Utree.rtree] = ...
                PatSymmFactorizationRecursion(Stree.rtree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        [L,U] = lu(full(A(idx,idx)));
        
        AUinv = A(actidx,idx)/U;
        ALinv = (L\A(idx,actidx))';
        
        A(actidx,actidx) = A(actidx,actidx) ...
            - sparse(A(actidx,idx)/A(idx,idx)*A(idx,actidx));
        
        Ltree.Mat = L;
        Ltree.AMatinv = AUinv;
        Utree.Mat = U';
        Utree.AMatinv = ALinv;
        
    end

end