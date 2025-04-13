function [M, remainedColIdx, remainedRowIdx, pmd] = dvtime( M, thrd )
% data partition and discard pixels without 7 points in time series
% Input: region_part, obervation matrix; day_matrix, the date matrix xpt;
% thrd1: threshold for removing bins from maxtrix - line 30
% -M: the matrix; flag, tag on removed location; 
% -interp_loc: a matrix to mark the place where need to be interpolated
% -Mr: removed rows are stored 
% -Mc: removed cols are stored
%
% Haipeng
% 6/13/2020
% 6/27/2020 
% 7/5/2020 input is from 4.64km datasets
% 7/16/2020 revision for GSM datasets
% 4/29/2021 
% 2/10/2023 Remove outputs Mr and Mc
% 5/25/2023 Fix line 49: change size(M, 2) to size(M, 1)

M = single( M );
% Mr = zeros( size( M, 1 ), 1, 'logical' );  % variable to mark removed rows

remainedColIdx = zeros( size( M, 2 ), 1, 'logical' ); % zero represents the column removed from data matrix for no enough data through a year
remainedRowIdx = zeros( size( M, 1 ), 1, 'logical' ); % zero represents the row removed from data matrix for no enough data through a year

T = zeros( size( M ), 'single' );
t = 0;
% Discard pixels without any measurements in time series which could
% comprises of pixels located in land and perennial sea-ice covered region
for k = 1 : size( M, 1 )
    bin_ts = M(k, :);
    if nnz( bin_ts ) < 2 % Pixels was identified as  land or perenial ice cover NOTE: here may need adjustment based ond input
%         Mr(k) = 1; do nothing
    elseif nnz( bin_ts ) / size( M, 2 ) < thrd % No enough points in the row
%         Mr(k) = 1; do nothing
    else % Keep the row where sufficient points are available
        t = t + 1;        
        T(t, :) = bin_ts;
        remainedRowIdx(k) = 1;        
    end        
end
M = T(1 : t, :);

% Discard column where the percent of valid measurement < thrd
% Mc = zeros( size( M, 2 ), 1, 'logical' );
T = zeros( size( M ), 'single' );
t = 0;
for k = 1 : size( M, 2 )
    oneCol = M(:, k);
    if nnz( oneCol ) / size( M, 1 ) < thrd
%         Mc(k) = 1; do nothing
    else  % Keep the column where sufficient points are available
        t = t + 1;
        T(:, t) = oneCol;
        remainedColIdx(k) = 1;
    end
end
M = T(:, 1 : t);
pmd = nnz( ~M ) / length( M(:) ) * 100; % Percentage of missing data
disp( ['PMD = ', num2str( pmd ), '%'] );

if length( M ) < 5 
    warning( ['The matrix is too small '] );
    M = [];
end




