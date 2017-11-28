function get_dvbs2_alist(R, alist_filename )
if ~exist(['scripts/sparse_matlab_matrix_to_alist.' mexext], 'file')
    if(ismac)
        mex  -v -Losx/lib  -Iinclude scripts/sparse_matlab_matrix_to_alist.cpp -outdir scripts -litpp_static_debug -lfftw3 -llapack -lblas
    elseif(isunix) %for linux, we can use a dynamic stock version of itpp installed in the systems default library paths as its in most distros repos
        mex  -v -L/usr/lib scripts/sparse_matlab_matrix_to_alist.cpp -outdir scripts -litpp
    else
        error('Rethink your life if you use windows for this kind of development ;) ');
end
H = dvbs2ldpc(R, 'sparse');

% degree_one_vn_idx = find(sum(H)==1);
% H(:,degree_one_vn_idx)=[];
[num_rows, num_cols] = size(H);
[row_idx, col_idx, ~] = find(H);
sparse_matlab_matrix_to_alist(num_rows, num_cols, row_idx-1, col_idx-1, alist_filename);

end