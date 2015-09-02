function qp = qp_extract_axis(qp, j)
    Nt = length(qp.H) - 1;
    Nq = size(qp.H{end}, 1) / 2;
    
    z_ind = [j, Nq + j, 2 * Nq + j];
    x_ind = [j, Nq + j];
    
    for i = 1 : Nt
        qp.H{i} = qp.H{i}(z_ind, z_ind);
        qp.g{i} = qp.g{i}(z_ind);
        qp.C{i} = qp.C{i}(x_ind, z_ind);
        qp.c{i} = qp.c{i}(x_ind);
        qp.zMin{i} = qp.zMin{i}(z_ind);
        qp.zMax{i} = qp.zMax{i}(z_ind);
    end
    
    qp.H{end} = qp.H{end}(x_ind, x_ind);
    qp.g{end} = qp.g{end}(x_ind);
    qp.zMin{end} = qp.zMin{end}(x_ind);
    qp.zMax{end} = qp.zMax{end}(x_ind);
end

