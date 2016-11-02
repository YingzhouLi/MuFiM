function C = mtimes(A,B)
% MTIMES C = A B
%   C = A*B multiplies matrix A and B together when either of them is
%   factorized.
%
%   See also SYMBOLMF, MULTIFRONTAL.

%   Copyright 2016 Yingzhou Li, Stanford University

if isa(A,'Multifrontal')
    
    MF = A;
    C = B;
    
    if MF.symm == 1
        LeftMulSymmRecursionUp  (MF.symboltree,MF.Ltree);
        LeftMulSymmRecursionDiag(MF.symboltree,MF.Dtree);
        LeftMulSymmRecursionDown(MF.symboltree,MF.Ltree);
    elseif MF.symm == 2
        LeftMulSymmRecursionUp  (MF.symboltree,MF.Utree);
        LeftMulSymmRecursionDown(MF.symboltree,MF.Ltree);
    end
    
elseif isa(B,'Multifrontal')
    
    MF = B;
    C = A;
    
    if MF.symm == 1
        RightMulSymmRecursionUp  (MF.symboltree,MF.Ltree);
        RightMulSymmRecursionDiag(MF.symboltree,MF.Dtree);
        RightMulSymmRecursionDown(MF.symboltree,MF.Ltree);
    elseif MF.symm == 2
        RightMulSymmRecursionUp  (MF.symboltree,MF.Ltree);
        RightMulSymmRecursionDown(MF.symboltree,MF.Utree);
    end
    
end

%=====================================================================
    function LeftMulSymmRecursionUp(Stree,Utree)
        
        if strcmpi(Stree.type,'node')
            LeftMulSymmRecursionUp(Stree.ltree,Utree.ltree);
            LeftMulSymmRecursionUp(Stree.rtree,Utree.rtree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(idx,:) = Utree.Mat'*C(idx,:);
        C(idx,:) = C(idx,:) + Utree.AMatinv'*C(actidx,:);
        
    end

    function LeftMulSymmRecursionDiag(Stree,Dtree)
        
        if strcmpi(Stree.type,'node')
            LeftMulSymmRecursionDiag(Stree.ltree,Dtree.ltree);
            LeftMulSymmRecursionDiag(Stree.rtree,Dtree.rtree);
        end
        
        idx = Stree.idx;
        C(idx,:) = Dtree.Mat*C(idx,:);
        
    end

    function LeftMulSymmRecursionDown(Stree,Ltree)
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(actidx,:) = C(actidx,:) + Ltree.AMatinv*C(idx,:);
        C(idx,:) = Ltree.Mat*C(idx,:);
        
        if strcmpi(Stree.type,'node')
            LeftMulSymmRecursionDown(Stree.rtree,Ltree.rtree);
            LeftMulSymmRecursionDown(Stree.ltree,Ltree.ltree);
        end
        
    end

%=====================================================================
    function RightMulSymmRecursionUp(Stree,Ltree)
        
        if strcmpi(Stree.type,'node')
            RightMulSymmRecursionUp(Stree.rtree,Ltree.rtree);
            RightMulSymmRecursionUp(Stree.ltree,Ltree.ltree);
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(:,idx) = C(:,idx)*Ltree.Mat;
        C(:,idx) = C(:,idx) + C(:,actidx)*Ltree.AMatinv;
        
    end

    function RightMulSymmRecursionDiag(Stree,Dtree)
        
        if strcmpi(Stree.type,'node')
            RightMulSymmRecursionDiag(Stree.ltree,Dtree.ltree);
            RightMulSymmRecursionDiag(Stree.rtree,Dtree.rtree);
        end
        
        idx = Stree.idx;
        C(:,idx) = C(:,idx)*Dtree.Mat;
        
    end

    function RightMulSymmRecursionDown(Stree,Utree)
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        C(:,actidx) = C(:,actidx) + C(:,idx)*Utree.AMatinv';
        C(:,idx) = C(:,idx)*Utree.Mat';
        
        if strcmpi(Stree.type,'node')
            RightMulSymmRecursionDown(Stree.ltree,Utree.ltree);
            RightMulSymmRecursionDown(Stree.rtree,Utree.rtree);
        end
        
    end

end