function [gex, dgex] = Gexs(Par, T, npa, y, dy) 
    Npa = Par.Npa;
    gex = 0;
    dgex = zeros(1,Npa);
    
    for b = 1:Par.topo.nBasePairs
        ids = Par.idx.npa_ids_of_base{b};
        
        if Par.topo.is_homo(ids(1))
            continue;
        end
        
        for k = ids
            [gb, dgb] = OTFBEF(Par, T, npa, idp, y, dy);
            gex = gex + npa(k) * gb;
            dgex = dgex + npa(k) * dgb; 
        end
    end
      
end