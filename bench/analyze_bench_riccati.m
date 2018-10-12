function analyze_bench_riccati(filename)
% Analyze the results from Riccati solvers benchmark.

    data = jsondecode(fileread(filename));
    c = cell2mat(arrayfun(@parse_benchmark, data.benchmarks, 'UniformOutput', false));
    
    algs = unique({c.alg});

    for nu = unique([c.nu])
        figure();
        hold('on');

        for alg = algs
            ind = strcmp({c.alg}, alg{:}) & strcmp({c.flag}, 'mean') & [c.nu] == nu;
            plot([c(ind).nx], [c(ind).cpu_time], 'LineWidth', 2);
            
            ind_ci = strcmp({c.alg}, alg{:}) & strcmp({c.flag}, 'stddev') & [c.nu] == nu;
%             ci = 3 * [c(ind_ci).cpu_time];
%             ciplot(([c(ind).cpu_time] - ci).', ([c(ind).cpu_time] + ci).', [c(ind).nx], 'g');
        end

        legend(algs);
        title(sprintf('nu=%d', nu));
        xlabel('nx');
        ylabel('t [ns]');
        grid('on');
    end
end


function s = parse_benchmark(b)
    s.cpu_time = b.cpu_time;
    [s.alg, s.N, s.nx, s.nu, s.flag] = parse_benchmark_name(b.name);
end


% function s = parse_benchmark_name(str)
%     parts = strsplit(str, {'/', '_', 'Riccati'});
%     s.alg = parts{2};
%     s.N = str2double(parts{3});
%     s.nx = str2double(parts{4});
%     s.nu = str2double(parts{5});
%     s.value = parts{6};
% end

function [alg, N, nx, nu, flag] = parse_benchmark_name(str)
    parts = strsplit(str, {'/', '_', 'Riccati'});
    alg = parts{2};
    N = str2double(parts{3});
    nx = str2double(parts{4});
    nu = str2double(parts{5});
    flag = parts{6};
end