function apchew


if 1
	c=get(gca,'ColorOrder');
	maxc = size(c, 1);
	clf
	
	str = {};
	h = [];

	fls = dir('Data/S2I10sp_*.txt');
	N = length(fls);
	
	for i=1:N
		nm = fls(i).name;		
		o = regexp(nm, 'S2I10sp_(\d+)_(\w+)_(\d+)_-(\d+).txt', 'tokens');
		freq          = str2double(o{1}{1});
		drug          = o{1}{2};
		concentration = str2double(o{1}{3});
		Eleak         = -str2double(o{1}{4});

		selected = 0;
		
		if strcmp('carbamazepine', drug)==1
			selected = selected | 1;
		end
		if strcmp('phenytoin', drug)==1
			selected = selected | 1;
		end
		if strcmp('none', drug)==1
			selected = selected | 1;
		end
		if concentration>20
			selected = 0;
		end
% 		if Eleak~=-70
% 			selected = 0;
% 		end
		if freq~=50
			selected = 0;
		end
		
		if ~selected, continue, end

		junk = load(['Data/' nm]);
		if isempty(junk), continue, end

		APs = junk(junk(:,2)<=499,1); % Granule cells only
		
		binsize = 100;
		tm = 0:binsize:max(APs)+binsize;
		dt = histc(APs, tm)/500*10;
		
        out(:, 1)   = tm(1:21)';
        out(:, end+1) = dt(1:21);
        fprintf(1, '%s %g\n', drug, Eleak);
        
		col = c(mod(i-1, maxc)+1, :);
		h(end+1) = line(tm, dt, ...
			'Color', col, ...
			'LineWidth', 2, ...
			'markersize', 2, ...
			'markerfacecolor', 'black');
		str{end+1} = sprintf('freq=%d [%s]=%gmM Eleak=%g', freq, drug, concentration, Eleak);
	end
	grid on
	set(gca, 'Xlim', [1 1900], 'Ylim', [0 50])
	legend(h, str)
	xlabel('Time (ms)')
	ylabel('Frequency (Hz)')
	toggleplot(gcf)
    
    save async.dat out -ASCII
end

if 0
	fls = dir('Data/M2I10sp_*.txt');

	for i=1:length(fls)
		nm = fls(i).name;
		A = sscanf(nm, 'M2I10sp_%d_%g.txt');
		freq = A(1);
		persistentNa = A(2);
	
		figure(i)
		figure('Name', sprintf('Freq=%d  PersistentNa=%g', freq, persistentNa*100), ...
			'NumberTitle', 'off')
		clf
		junk = dlmread(['Data/' nm], '', 2, 1);
		toggleplot(junk(:,1:6))
	end
end