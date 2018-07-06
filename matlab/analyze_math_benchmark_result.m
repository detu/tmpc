function analyze_math_benchmark_result(filename)
    j = jsondecode(fileread(filename));
    
    n = [2, 3, 5, 10, 20, 30];
    method = {'CasADi MX', 'CasADi SX', 'Eigen', 'Blaze', 'blasfeo'};
    cpu_time = zeros(length(n), length(method));
    
    for i_n = 1 : length(n)
        for i_m = 1 : length(method)
            name = benchmark_name(method{i_m}, n(i_n));
            data = j.benchmarks(strcmp({j.benchmarks.name}, name));
            cpu_time(i_n, i_m) = data.cpu_time;
        end
    end
    
    flops = repmat((n.^3).', 1, length(method)) ./ cpu_time;
    
    bar(categorical(n), flops);
    grid('on');
%     set(gca, 'YScale', 'log');
    legend(method, 'Location', 'NorthWest');
    xlabel('Matrix size');
    ylabel('Performance [Gflops]');
    title('Benchmark C=A*B');
end


function name = benchmark_name(method, n)
    switch (method)
        case 'CasADi MX'
            name = sprintf('BM_MTimes/%dx%dx%d_MX_median', n, n, n);
        case 'CasADi SX'
            name = sprintf('BM_MTimes/%dx%dx%d_SX_median', n, n, n);
        case 'Eigen'
            name = sprintf('BM_MTimes<EigenKernel<double>, %d, %d, %d>_median', n, n, n);
        case 'Blaze'
            name = sprintf('BM_MTimes<BlazeKernel<double>, %d, %d, %d>_median', n, n, n);
        case 'blasfeo'
            name = sprintf('BM_MTimes/%dx%dx%d_blasfeo_median', n, n, n);
        otherwise
            error('Unknown method %s', method);
    end
end