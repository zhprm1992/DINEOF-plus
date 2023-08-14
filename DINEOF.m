function [Xa_best, neof_best, rms_best] = DINEOF( Xo, neof_max, rms_d )
% Data interpolating Empirical Orthogonal Functions 
% Input:
% - X0, a gappy data matrix; 
% - drms, threshold of convergence
% - interloc, Positions of interpolated values.
% Output:
% - Xa_best, observed data with filled points; 
% - neof_best, number of eofs used in interpolation; 
% - rms_best, a vector of the RMS values from the iteration
%
% Haipeng Zhao
% Created 8/14/2023

    rms_best = inf;
    neof_best = 0;
    
    idx_gap = uint32( find( Xo == 0 ) );
    idx_known = uint32( find( Xo ~= 0 ) );

    num_valid = max( 30, floor( 0.01 * length( idx_known ) ) );
    idx_valid = idx_known( randperm( length( idx_known ), num_valid ) );

    Xa = Xo;
    Xa(idx_valid) = 0;
    Xo_idx_valid = Xo(idx_valid);
    clear Xo

    rms_pre = Inf;
    error_now = vpa( sqrt( sym( mean( (Xa(idx_valid) - Xo_idx_valid) .^ 2 ) ) ), 5 );
    error0 = error_now;
    neof = 1;
    
    % loop for determining optimal number of EOFs 
    Xr = zeros( size( Xa ), 'single' );
    Dmax = zeros( neof_max, 1, 'single' );

    while (neof < neof_max)
        % loop for interpolation until convergency
        while (rms_pre - error_now(end) > rms_d)  
            % replace missing points
            for i = 1 : length( idx_gap )
                Xa(idx_gap(i)) = Xr(idx_gap(i));
            end
            for i = 1 : length( idx_valid )
                Xa( idx_valid(i) ) = Xr( idx_valid(i) );
            end                        
            rms_pre = error_now(end);
            % calculation of EOFs to obatian restructed image
            [U, D, V] = svd( Xa, 'econ');                       
            x = U(:, 1 : neof) * D(1 : neof, 1 : neof) * V(:, 1 : neof)'; 
            Xr = x; 
            error = vpa( sqrt( sym( mean( (Xr(idx_valid) - Xo_idx_valid) .^ 2 ) ) ), 5 ); % using RMSE
            error_now(end + 1) = error;
        end
        if rms_best(end) >  error_now(end) || abs( rms_best(end) -  error_now(end) ) < rms_d
            rms_best(end + 1) = error_now(end);
            Xa_best = Xa;
            neof_best = neof;
            rms_pre = inf;
            error_now(end) = error0;          
            Dmax(neof) = D(neof, neof);
            neof = neof + 1;       
        else
            break;  
        end
    end
    Xa_best( idx_valid ) = Xo_idx_valid;
    rms_best = rms_best(end);
end

