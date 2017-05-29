function [M, F] = compute_moments(summat)

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

	    F = (M1/0.1 - 1)^2  + (M2/0.35 - 1)^2 + (M4/0.03 - 1)^2 + (M5/0.5 - 1)^2 + (M7/0.55 - 1)^2 ; %+ (M8/34 - 1)^2 ;
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

	    F = (M1-0.4)^2  + (M2-0.25)^2 + (M4-0.03)^2 + (M5-0.45)^2 + (M7-0.02)^2 + (M8-0.75)^2 + ...
	    	(M10-0.01)^2 + (M11-0.8)^2 + (M13-0.55)^2; % + (M14/34 - 1)^2 ;
	    % reuse fractions are not independent from other fraction, so not included
	    M = [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M12, M13, M14];
	end

    
    T = table;
    T.moment = {'model'; 'data'};
    T.new_tech1850 = [M1 ; 0.4];
    T.new_comb1850 = [M2 ; 0.25];
    T.reuse1850 = [M3 ; 0.35];
    T.new_tech1900 = [M4 ; 0.03];
    T.new_combo1900 = [M5 ; 0.45];
    T.reuse1900 = [M6 ; 0.52];

    T2 = table;
    T2.moment = {'model'; 'data'};
    T2.new_tech1950 = [M7 ; 0.02];
    T2.new_comb1950 = [M8 ; 0.75];
    T2.reuse1950 = [M9 ; 0.23];
    T2.new_tech2000 = [M10 ; 0.01];
    T2.new_comb2000 = [M11 ; 0.80];
    T2.reuse2000 = [M12 ; 0.19];

    T3 = table;
    T3.moment = {'model'; 'data'};
    T3.reuse_peak = [M13; 0.55];
    %T3.peak_year = [M14 ; 34];

    display(T)
    display(T2)
    display(T3)

    data = [0.4; 0.25; 0.35; 0.03; 0.45; 0.52; 0.02; 0.75; 0.33; 0.01; 0.80; 0.19; 0.55];% 34];
    Summary_stats = [data, M(1:end-1)'];

    rowLabels = {'new tech 1850'; 'new comb 1850'; 'reuse 1850'; 'new tech 1900'; 'new comb 1900'; 'reuse 1900'; ...
        'new tech 1950'; 'new comb 1950'; 'reuse 1950'; 'new tech 2000'; 'new comb 2000'; 'reuse 2000'; 'reuse peak'};% 'peak year'};
    columnLabels = {'Data', 'Model'};
    matrix2latex(Summary_stats, 'tex_files/moments_table.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels)

end