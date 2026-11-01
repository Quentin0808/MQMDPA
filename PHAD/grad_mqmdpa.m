function dG = grad_mqmdpa(Par, T, npa)
    [~, dG] = G_mqmdpa(Par, T, npa);
    dG = dG(:);
end

