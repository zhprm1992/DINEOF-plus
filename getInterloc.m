function inter_loc = getInterloc( chl_gap, ind, mzsize )
% Identify location of interpolation. When two conditions are satisfied,
% the unobserved value is identified as a values that can be estimated.
%
% Input:
%   chl_gap - data matrix with gaps.
%   idx - index of pixels in orginal image
%   mzsize - size of original image

% Output;
%   inter_loc - Location of values that are able to be estimated (markded by 1)
%     
% Haipeng Zhao
% 2/8/2023
% 2/14/2023
% 4/7/2026

[u1, u2] = size( chl_gap );
pos_ob2D = ~isnan( chl_gap ); % Observed index

% Calculate temporal connectivity, the maximum number of unobsrved values
% is set 6.
neighbor1 = [pos_ob2D(:, 2 : end), zeros( u1, 1, 'logical' )]; % +1
neighbor2 = [pos_ob2D(:, 3 : end), zeros( u1, 2, 'logical' )]; % +2
neighbor3 = [pos_ob2D(:, 4 : end), zeros( u1, 3, 'logical' )]; % +3
neighbor4 = [zeros( u1, 1, 'logical' ), pos_ob2D(:, 1 : end - 1)]; % -1
neighbor5 = [zeros( u1, 2, 'logical' ), pos_ob2D(:, 1 : end - 2)]; % -2
neighbor6 = [zeros( u1, 3, 'logical' ), pos_ob2D(:, 1 : end - 3)]; % -3
tmpCnn = (neighbor1 | neighbor2 | neighbor3 | neighbor4 | neighbor5 | neighbor6) & ~pos_ob2D;
clear neighbor1 neighbor2 neighbor3 neighbor4 neighbor5 neighbor6

% Calculate spatial connectivity, the maximum number of unobsrved values
% is set 8.
[d1, d2] = ind2sub( mzsize, ind );
d1 = uint16( d1 );
d2 = uint16( d2 );
neighbor1 = [d1 + 1, d2 + 1]; % +1 +1
neighbor2 = [d1 + 1, d2]; % +1 +0
neighbor3 = [d1 + 1, d2 - 1]; % +1 -1
neighbor4 = [d1, d2 + 1]; % 0 +1
neighbor5 = [d1, d2 - 1]; % 0 -1
neighbor6 = [d1 - 1, d2 + 1]; % -1 +1
neighbor7 = [d1 - 1, d2]; % -1 0
neighbor8 = [d1 - 1, d2 - 1]; % -1 -1
for i = 1 : length( d1 )
    if neighbor1(i, 1) > mzsize(1)
        neighbor1(i, 1) = 1;
        neighbor2(i, 1) = 1;
        neighbor3(i, 1) = 1;
    end
    if neighbor6(i, 1) == 0
        neighbor6(i, 1) = mzsize(1);
        neighbor7(i, 1) = mzsize(1);
        neighbor8(i, 1) = mzsize(1);
    end
    if neighbor1(i, 2) > mzsize(2)
        neighbor1(i, 2) = mzsize(2);
        neighbor4(i, 2) = mzsize(2);
        neighbor6(i, 2) = mzsize(2);
    end
    if neighbor3(i, 2) == 0
        neighbor3(i, 2) = 1;
        neighbor5(i, 2) = 1;
        neighbor8(i, 2) = 1;
    end
end
neighbor1 = sub2ind( mzsize, neighbor1(:, 1), neighbor1(:, 2) );
[~, loc1] = ismember( neighbor1, ind );
clear neighbor1
neighbor2 = sub2ind( mzsize, neighbor2(:, 1), neighbor2(:, 2) );
[~, loc2] = ismember( neighbor2, ind );
clear neighbor2
neighbor3 = sub2ind( mzsize, neighbor3(:, 1), neighbor3(:, 2) );
[~, loc3] = ismember( neighbor3, ind );
clear neighbor3
neighbor4 = sub2ind( mzsize, neighbor4(:, 1), neighbor4(:, 2) );
[~, loc4] = ismember( neighbor4, ind );
clear neighbor4
neighbor5 = sub2ind( mzsize, neighbor5(:, 1), neighbor5(:, 2) );
[~, loc5] = ismember( neighbor5, ind );
clear neighbor5
neighbor6 = sub2ind( mzsize, neighbor6(:, 1), neighbor6(:, 2) );
[~, loc6] = ismember( neighbor6, ind );
clear neighbor6
neighbor7 = sub2ind( mzsize, neighbor7(:, 1), neighbor7(:, 2) );
[~, loc7] = ismember( neighbor7, ind );
clear neighbor7
neighbor8 = sub2ind( mzsize, neighbor8(:, 1), neighbor8(:, 2) );
[~, loc8] = ismember( neighbor8, ind );
clear neighbor8

sptCnn = zeros( size( tmpCnn ), 'logical' );
for i = 1 : u2
    for j = 1 : u1
        if ~pos_ob2D(j, i)
            % Get neighbors in spatial dimension
            if loc1(j) == 0 || loc2(j) == 0 || loc3(j) == 0 || loc4(j) == 0 || loc5(j) == 0 ||...
                    loc6(j) == 0 || loc7(j) == 0 || loc8(j) == 0
                    neighbor = [loc1(j); loc2(j); loc3(j); loc4(j); loc5(j); loc6(j); loc7(j); loc8(j)];
                    neighbor = neighbor( ( neighbor ) > 0 );
                    sptCnn(j, i) = any( pos_ob2D( neighbor, i ) );
            else
                sptCnn(j, i) = pos_ob2D(loc1(j), i) | pos_ob2D(loc2(j), i) | pos_ob2D(loc3(j), i) | pos_ob2D(loc4(j), i) ...
                    | pos_ob2D(loc5(j), i) | pos_ob2D(loc6(j), i) | pos_ob2D(loc7(j), i) | pos_ob2D(loc8(j), i);
            end
        end
    end
end
inter_loc = tmpCnn | sptCnn;

end

