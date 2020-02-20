directory = 'V:\Neurometric\2018\Test-Retest';
datadir=dir(directory);
datadir=struct2table(datadir);

%good 3D scan subjects
% [~, dreid] = %xlsread('C:\Users\wults\OneDrive\Dokumente\Masterarbeit\good3dscan.xlsx');
%dreidgood = %readtable('C:\Users\wults\OneDrive\Dokumente\Masterarbeit\good3dscan.xlsx');

a=[];
for iSub=4:300
    
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
            pattern = "OFFSET (\d\.\d\d) deg"; %the degrees after offset in a cell %we are only interested in <1.00 deg? --> this is still missing
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
