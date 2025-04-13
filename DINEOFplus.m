function  [Xr, rms]  = DINEOFplus( Xo, longitudeLength, latitudeLength, idx, thrd )
% DINEOF PLUS. This version includes calling for function getinterloc (i.e. connectivity filter)
% Inputs:
% - Xo:  Incomplete matrix, dimemsion = (space, time)
% - longitudeLength: length of longitude vector
% - latitudeLength: length of latitude vector
% - idx: index of observed pixels
% - thrd: minimum percentage of column/row based missing data
% Outputs:
% - Xr: recovered matrix
% - rms
% 
% Haipeng Zhao
% 2/8/2023
% 10/1/2023 1. Update DINEOF 2. run multiple times of DINEOF

interloc = getInterloc_v2( Xo, idx, [longitudeLength, latitudeLength], 8, 6 );

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
    % log10 transform
    Xs(Xs > 0) = log10( Xs(Xs > 0) );
    
    idx_ob = find( Xs ~= 0 ); % Get positions of observed entries
    x = Xs(idx_ob);
    avg = mean( x(:) );
    Xs(idx_ob) = x - avg; % substract mean value from entire non-zero data
    clear x idx_ob

    eof_max = min( [50, min( size( Xs ) )]);
    [Xr, ~, rms] = DINEOF( Xs, eof_max, 1e-5 );
    for i = 3  % Note: SET UP THE LOOP NUMBER
        [Xr2, ~, rms2] = DINEOF( Xs, eof_max, 1e-5 );
        if rms2 < rms
            Xr = Xr2;
            rms = rms2;
        end
    end

    % The back-transformation from log10 space 
    idx = Xr ~= 0;
    Xr(idx) = Xr(idx) + avg;
    Xr = 10 .^ Xr;
    Xr( Xr == 1 ) = 0;
    
  %% insert non-interpolation row into reconstructed matrix and extract bin index with valide data 
    Xd = Xo; 
    Xd(remainedRowIdx, remainedColIdx) = Xr;   

    % Remove outlier to redue errors
    entries_rct = log10( Xd(interloc) ); 
    [~, idx] = rmoutliers( ( entries_rct ), 'percentiles', [0.01, 99.99] );
    tmp = find( interloc > 0 );
    idx = tmp(idx);
    interloc(idx) = 0; 

    Xo(interloc) = Xd(interloc);
    Xr = Xo;
else
    warning( 'Not enough data for reconstruction' )
    rms = 0;
    Xr = Xo;
end

Xr(Xr == 0) = nan;
%%
% figure
% imagesc( interloc > 0 )
% colormap( 'gray' )
% title( 'Locations for interpolation' )

end
