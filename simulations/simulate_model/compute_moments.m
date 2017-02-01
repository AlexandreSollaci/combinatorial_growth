function [M, F] = compute_moments(summat)

% This function computes the model moments, M, and the value of the SMM value function.
%
% The input matrix, summat, has the following form: 
%			[year, fraction of new technolofy patents in year, same for new combinations, same for reuse]
%
% This function works for two different simluations, the first with only part of the whole period (typically 100 years)
% and the second with the full 200 years.
%
% Moments are:
% 100 year simulation: fraction of patents on 1880 and 1930, reuse fraction peak value and year
% 200 year simulation: fraction of patents on 1850, 1900, 1950 and 2000, reuse fraction peak value and year
%


	if summat(end,1) < 1950

	    index1 = find(summat(:,1) == 1880);
	    index2 = find(summat(:,1) == 1930);

	    M1 = summat(index1, 2);
	    M2 = summat(index1, 3);
	    M3 = summat(index1, 4);

	    M4 = summat(index2, 2);
	    M5 = summat(index2, 3);
	    M6 = summat(index2, 4);
	    [M7 , M8] = max(summat(:,4));

	    F = (M1/0.1 - 1)^2  + (M2/0.35 - 1)^2 + (M4/0.03 - 1)^2 + (M5/0.5 - 1)^2 + (M7/0.55 - 1)^2 + (M8/34 - 1)^2 ;
	    % M3 and M6 are not independent from the others.
	    M = [M1, M2, M3, M4, M5, M6, M7, M8];

	else

	    index1 = find(summat(:,1) == 1850);
	    index2 = find(summat(:,1) == 1900);
	    index3 = find(summat(:,1) == 1950);
	    index4 = find(summat(:,1) == 2000);

	    % patent type fractions, 1850
	    M1 = summat(index1, 2);
	    M2 = summat(index1, 3);
	    M3 = summat(index1, 4);
	    % patent type fractions, 1900
	    M4 = summat(index2, 2);
	    M5 = summat(index2, 3);
	    M6 = summat(index2, 4);
	    % patent type fractions, 1950
	    M7 = summat(index3, 2);
	    M8 = summat(index3, 3);
	    M9 = summat(index3, 4);
	    % patent type fractions, 200
	    M10 = summat(index4, 2);
	    M11 = summat(index4, 3);
	    M12 = summat(index4, 4);
	    % Peak of reuse fraction
	    [M13 , M14] = max(summat(:,4)); 

	    F = (M1/0.4 - 1)^2  + (M2/0.25 - 1)^2 + (M4/0.03 - 1)^2 + (M5/0.45 - 1)^2 + (M7/0.02 - 1)^2 + (M8/0.75 - 1)^2 + ...
	    	(M10/0.01 - 1)^2 + (M11/0.8 - 1)^2 + (M13/0.55 - 1)^2 + (M14/34 - 1)^2 ;
	    % reuse fractions are not independent from other fraction, so not included
	    M = [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M12, M13, M14];
	end

end