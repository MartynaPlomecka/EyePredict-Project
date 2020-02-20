directory = 'V:\Neurometric\2018\Test-Retest';
datadir=dir(directory);
datadir=struct2table(datadir);

a=[];
for iSub=4:389
    
    sub = char(datadir.name(iSub));
    stages = [];
    for i=1:5
        try
        deg = [];
        D = importdata([directory, '\', sub, '\ARCHIVE\', sub, '_AS', num2str(i), '.asc'], '\t', 200); % Path to .asc file
       
        try
            D = D.textdata;
        catch
           
        end

        %%
        % looking for the OFFSET in the asc file
        for j=1:length(D)
            line = char(D(j));
            pattern = "OFFSET (\d\.\d\d) deg"; 
            t = regexp(line, pattern, 'tokens');
            if isempty(t)
                continue
            else
                deg = [deg, str2num(string(t{:}))];
            end
        end
        catch
        end
        stages = [stages, mean(deg)];
        
    end
    a=[a; {sub, stages}]
    
    
    
end
