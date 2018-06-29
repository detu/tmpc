function qp = json2qp(json)
    % Build a qp object.
    
    n_nodes = length(json.nodes);
    
    for k = n_nodes : -1 : 1
        src_node = json.nodes(k);
        nx = length(src_node.q);
        nu = length(src_node.r);
        nc = length(src_node.ld);
        
        dst_node = QpNode(nx, nu, nc);
        
        for field = {'Q', 'R', 'S', 'q', 'r', 'C', 'D', 'ld', 'ud', 'lx', 'ux', 'lu', 'uu'}
            if ~(isempty(src_node.(field{:})) && isempty(dst_node.(field{:})))
                dst_node.(field{:}) = src_node.(field{:});
            end
        end
        
        qp.nodes(k) = dst_node;        
    end
    
    n_edges = length(json.edges);
    
    for k = n_edges : -1 : 1
        % Convert 0-based node indices to 1-based.
        src_edge = json.edges(k);
        dst_edge.from = src_edge.from + 1;
        dst_edge.to = src_edge.to + 1;
        
        qp.edges(k) = dst_edge;
    end
end