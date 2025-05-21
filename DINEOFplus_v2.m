function  [Xr, rms, error_var]  = DINEOFplus_v2( Xo, d1, d2, idx, thrd, logFlag )
% DINEOF PLUS. This version includes calling for function getinterloc (i.e. connectivity mask)
% Note: this version include expected errror variance and only run DINEOF once 
% and doesn't include remove outlier 
% Inputs:
% - Xo:  Incomplete matrix, dimemsion = (space, time)
% - longitudeLength: length of longitude vector
% - latitudeLength: length of latitude vector
% - idx: index of observed pixels
% - thrd: minimum percentage of column/row based missing data
% - logFlag: 0 indicates no log transform, 1 indicates use of log transform
% Outputs:
% - Xr: recovered matrix
% - rms: 
% 
% Haipeng Zhao
% 2/8/2023
% 10/1/2023 1. Update DINEOF 2. run multiple times of DINEOF
% 1/30/2024  Add expected error variance (https://doi.org/10.5194/os-2-183-2006)
% 7/2/2024 Add input logFlat for indicating log-transform

interloc = getInterloc_v2( Xo, idx, [d1, d2], 8, 6 );

for i = 1 : numel( Xo )
    if isnan( Xo(i) )
        Xo(i) = 0;
    end
end

% remove low  PMD < 75%
pmd = nnz( ~Xo ) / length( Xo(:) ) * 100; % Percentage of missing data
disp( ['PMD = ', num2str( pmd ), '%'] );
[Xs, remainedColIdx, remainedRowIdx, pmd] = dvtime( Xo, thrd );
while pmd > 75
    thrd = thrd + 0.01;
    [Xs, remainedColIdx, remainedRowIdx, pmd] = dvtime( Xo, thrd );
    if isempty( Xs )
        break
    end
end

if ~isempty( Xs ) && nnz( Xs ) > 1000
    if logFlag == 1
        Xs(Xs > 0) = log10( Xs(Xs > 0) ); % log10 transform
    end

    idx_ob = find( Xs ~= 0 ); % Get positions of observed entries
    x = Xs(idx_ob);
    avg = mean( x(:) );
    Xs(idx_ob) = x - avg; % substract mean value from entire non-zero data
    clear x
    eof_max = min( [50, min( size( Xs ) )]);
    [Xr, neof, rms, ~] = DINEOF( Xs, eof_max, 1e-5 );
    for i = 5  % SET UP THE NUMBER OF LOOPs
        [Xr2, neof2, rms2, ~] = DINEOF( Xs, eof_max, 1e-5 );
        if rms2 < rms
            Xr = Xr2;
            rms = rms2;
            neof = neof2;
        end
    end

    %% Estimate error variance at each intrpolation by OI method
    error_var = nan( size( Xr ) );
    [U, D, V] = svd( Xr, 'econ');    
    Xn = single( U(:, 1 : neof) * D(1 : neof, 1 : neof) * V(:, 1 : neof)' ); 
    L = U(:, 1 : neof) * D(1 : neof, 1 : neof) / sqrt( size( Xs, 2 )); 
    I = eye( neof ); % 
    mu_sq = 1 / nnz( Xs ) * sum( (Xs(idx_ob) .^ 2 - Xn(idx_ob) .^ 2) );
    r = 50; % Calibrate error correlation
    clear U D V Xn
% This commented lines are used to determine calibration coefficient r
% using approach 2 (see reference)
%     [u, v] = ind2sub( size( Xr ), idx_valid );
%     uv = [u v];
%     [uv, index] = sortrows( uv, 2 );
%     v = unique( uv(:, 2) );
%     t = zeros( 10, 1 );
%     for r = 1 : 10              
%         oi = [];
%         for i = 1 : length( v )
%             idx = Xs(:, v(i)) > 0;
%             Lp = L(idx, :);
%             D = inv( (Lp' * Lp + r * mu_sq * I) );
%             u = uv(uv(:, 2) == v(i), 1);
%             Lb = L(u, :)';
%             oi = [oi; diag( Lb' * D * Lp' * Lp * Lb )];
%         end
%         t(r) = mae( oi -  Xr( idx_valid(index) ) );
%     end

    for i = 1 : size( Xr, 2 )
        idx = Xs(:, i) > 0; % Indicate locations of data covered points
        Lp = L(idx, :);
        C = mu_sq * inv( (Lp' * Lp + r * mu_sq * I) );
        for k = 1 : size( Xr, 1 )    
            ib = L(k, :)';
            error_var(k, i) = ib' * C * ib;
        end
    end
    error_var = single( error_var );
%     error_std = 10 .^ error_std;
  %% Insert non-interpolation row into reconstructed matrix and extract bin index with valide data 
    idx = Xr ~= 0;
    Xr(idx) = Xr(idx) + avg;
    if logFlag == 1
        Xr = 10 .^ Xr;          % The back-transformation from log10 space  
        Xr( Xr == 1 ) = 0;
    else
        Xr( Xr < 0 ) = 0;
    end
    
    Xd = Xo; 
    Xd(remainedRowIdx, remainedColIdx) = Xr;   

%     % Remove outlier to redue errors
%     entries_rct = log10( Xd(interloc) ); 
%     [~, idx] = rmoutliers( ( entries_rct ), 'percentiles', [0.01, 99.99] );
%     tmp = find( interloc > 0 );
%     idx = tmp(idx);
%     interloc(idx) = 0; 

    Xo(interloc) = Xd(interloc);
    Xr = Xo;

    error_var0 = single( zeros( size( Xo )) );
    error_var0(remainedRowIdx, remainedColIdx) = error_var;   
    error_var = nan( size( error_var0 ) );
    error_var(interloc) = error_var0(interloc);
    
else
    warning( 'Not enough data for reconstruction' )
    rms = 0;
    Xr = Xo;
    error_var = [];
end
Xr(Xr == 0) = nan;

end

