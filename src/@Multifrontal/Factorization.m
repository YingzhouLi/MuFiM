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

P = zeros(size(A,1),1);

if MF.symm == 1
    [MF.Ltree,MF.Dtree] = ...
        SymmFactorizationRecursion(MF.symboltree);
elseif MF.symm == 2
    [MF.Ltree,MF.Utree] = ...
        PatSymmFactorizationRecursion(MF.symboltree);
end

%====================================================================
    function [Ltree,Dtree,extidx,Aupdate] = ...
            SymmFactorizationRecursion(Stree)
        
        if strcmpi(Stree.type,'node')
            [Ltree.ltree,Dtree.ltree,lidx,lA] = ...
                SymmFactorizationRecursion(Stree.ltree);
            
            [Ltree.rtree,Dtree.rtree,ridx,rA] = ...
                SymmFactorizationRecursion(Stree.rtree);
            
            [cidx,cA] = MergeUpdate(lidx,lA,ridx,rA);
        else
            cidx = [];
            cA = [];
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        
        [Aidx,Aact,Aupdate] = ExtUpdate(idx,actidx,cidx,cA);
        
        [L,D] = ldl(full(Aidx));
        
        ALDinv = Aact/L'/D;
        
        extidx = actidx;
        Aupdate = Aupdate - ALDinv/L*Aact';
        
        Ltree.Mat = L;
        Ltree.AMatinv = ALDinv;
        Dtree.Mat = D;
        
    end

%====================================================================
    function [Ltree,Utree,extidx,Aupdate] = ...
            PatSymmFactorizationRecursion(Stree)
        
        if strcmpi(Stree.type,'node')
            [Ltree.ltree,Utree.ltree,lidx,lA] = ...
                PatSymmFactorizationRecursion(Stree.ltree);
            
            [Ltree.rtree,Utree.rtree,ridx,rA] = ...
                PatSymmFactorizationRecursion(Stree.rtree);
            
            [cidx,cA] = MergeUpdate(lidx,lA,ridx,rA);
        else
            cidx = [];
            cA = [];
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        
        [Aidx,Aactidx,Aidxact,Aupdate] = PatExtUpdate(idx,actidx,cidx,cA);
        
        [L,U] = lu(Aidx);
        
        AUinv = Aactidx/U;
        ALinv = (L\Aidxact)';
        
        extidx = actidx;
        Aupdate = Aupdate - AUinv*ALinv';
        
        Ltree.Mat = L;
        Ltree.AMatinv = AUinv;
        Utree.Mat = U';
        Utree.AMatinv = ALinv;
        
    end

%====================================================================
    function [idx,A] = MergeUpdate(lidx,lA,ridx,rA)
        
        llen = length(lidx);
        rlen = length(ridx);
        [idx,~,IC] = unique([lidx,ridx],'sorted');
        lI = IC(1:llen);
        rI = IC(llen+1:llen+rlen);
        A = zeros(length(idx));
        A(lI,lI) = A(lI,lI) + lA;
        A(rI,rI) = A(rI,rI) + rA;
        
    end

    function [Aidx,Aact,Aupdate] = ExtUpdate(idx,actidx,cidx,cA)
        
        ilen = length(idx);
        alen = length(actidx);
        [I,J,V] = find(A(:,idx));
        
        P(idx) = 1:ilen;
        Aidx = zeros(ilen,ilen);
        iI = ismembc(I,idx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aidx(P(Isub) + (Jsub - 1)*ilen) = Vsub;
        
        P(actidx) = 1:alen;
        Aact = zeros(alen,ilen);
        iI = ismembc(I,actidx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aact(P(Isub) + (Jsub - 1)*alen) = Vsub;
        
        Aupdate = zeros(length(actidx));
        
        iicI = ismembc(idx,cidx);
        cicI = ismembc(cidx,idx);
        aacI = ismembc(actidx,cidx);
        cacI = ismembc(cidx,actidx);
        Aidx(iicI,iicI) = Aidx(iicI,iicI) + cA(cicI,cicI);
        Aact(aacI,iicI) = Aact(aacI,iicI) + cA(cacI,cicI);
        Aupdate(aacI,aacI) = cA(cacI,cacI);
        
    end

    function [Aidx,Aactidx,Aidxact,Aupdate] = ...
            PatExtUpdate(idx,actidx,cidx,cA)
        
        ilen = length(idx);
        alen = length(actidx);
        [I,J,V] = find(A(:,idx));
        
        P(idx) = 1:ilen;
        Aidx = zeros(ilen,ilen);
        iI = ismembc(I,idx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aidx(P(Isub) + (Jsub - 1)*ilen) = Vsub;
        
        P(actidx) = 1:alen;
        Aactidx = zeros(alen,ilen);
        iI = ismembc(I,actidx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aactidx(P(Isub) + (Jsub - 1)*alen) = Vsub;
        
        P(idx) = 1:ilen;
        [I,J,V] = find(A(:,actidx));
        Aidxact = zeros(ilen,alen);
        iI = ismembc(I,idx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aidxact(P(Isub) + (Jsub - 1)*ilen) = Vsub;
        
        Aupdate = zeros(length(actidx));
        
        iicI = ismembc(idx,cidx);
        cicI = ismembc(cidx,idx);
        aacI = ismembc(actidx,cidx);
        cacI = ismembc(cidx,actidx);
        Aidx(iicI,iicI) = Aidx(iicI,iicI) + cA(cicI,cicI);
        Aactidx(aacI,iicI) = Aactidx(aacI,iicI) + cA(cacI,cicI);
        Aidxact(iicI,aacI) = Aidxact(iicI,aacI) + cA(cicI,cacI);
        Aupdate(aacI,aacI) = cA(cacI,cacI);
        
    end

end