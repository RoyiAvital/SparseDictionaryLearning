function Display_D(D)

[a_len d_len] = size(D);
p_len         = sqrt(a_len);
p_pad_len     = p_len + 1;
lines_num     = ceil( sqrt(d_len) );
n             = 1:(lines_num*p_pad_len); n(p_pad_len:p_pad_len:end) = [];

D(:, end+1:lines_num^2) = min(D(:));
D2 = col2im( D, [p_len p_len], [ p_len * lines_num p_len * lines_num], 'distinct' );
D3 = min(D(:)) * ones( lines_num * p_pad_len , p_pad_len * lines_num );
D3(n,n) = D2;

imshow(D3,[]);


end
