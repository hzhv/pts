function matlab_solve(filename)

if nargin < 1
    error('Usage: matlab_solve(''filename'')');
end

fid = fopen(filename, 'r');
n = fscanf(fid, '%d', 1);    
row_ptr_size = fscanf(fid, '%d', 1);
row_ptr = fscanf(fid, '%d', row_ptr_size);

ci_size = fscanf(fid, '%d', 1);
col_id  = fscanf(fid, '%d', ci_size);

v_size = fscanf(fid, '%d', 1);
val    = fscanf(fid, '%f', v_size);
fclose(fid);

rows = zeros(ci_size,1);
for i = 1:n
    for k = row_ptr(i)+1:row_ptr(i+1)
        rows(k) = i;
    end
end
cols = col_id + 1;
A = sparse(rows, cols, val, n, n);

density = nnz(A) / (n * n);

%figure('Position',[100 100 600 600]);
%subplot(2,1,1)
%spy(A);
%title(sprintf('Sparsity Pattern (n=%d, density=%.6f)', n, density));

%% Level-set assignment
L = A;
%pred = L - diag(diag(L)); 
%level = zeros(n,1);
%for k = 1:n
%    preds = find(pred(k,1:k-1));
%    if isempty(preds)
%        level(k) = 0;
%    else
%        level(k) = 1;
%    end
%end

%subplot(2,1,2)
%scatter(1:n, level, 10, level, 'filled');
%xlabel('Row index');
%ylabel('Level (0 or 1)');
%title('Level-set assignment for each row');
%ylim([-0.5 1.5]);
%yticks([0 1]);

%fname = sprintf('sp_n%d.png', n);
%print(fname, '-dpng');

% Start Solving
b = ones(n, 1);
x = zeros(n, 1);

total_time = 0;
for i = 1:1000
    tic
    x = L \ b;
    endtime = toc;
    total_time = total_time + endtime;
end
fprintf('Matlab solves this using %0.16f sec\n', total_time/1000.0);

end