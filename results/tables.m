for index = 1:4
    
    switch index
        case 1
            fold  = '01 approximal';
            title = 'original parameters, approximal (i.e., original data)';
        case 2
            fold  = '01 proximal';
            title = 'original parameters, proximal';
        case 3
            fold  = '02 approximal';
            title = 'new parameters, approximal';
        case 4
            fold  = '02 proximal';
            title = 'new parameters, proximal';
    end
    
    load([fold, '/Experiments_DRvsFistavsSynvsAna.mat'])

    res = squeeze(max(res,[],5));
    transvec = [{'gab'},{'erb'},{'wav'}];
    synthesisvec = [{'synthsis'}, {'analysis'}];

    fprintf(['\n\n', title, '\n']);
    for m = 1:4
        fprintf('\n')
        fprintf('--------------------------------- %d ------------------------\n',m);
        fprintf('           DR / CP                | FISTA\n');
        fprintf('           GAB     WAV     ERB    | GAB     WAV     ERB\n');
        for k = 1:2
            fprintf([synthesisvec{k} ': ']);
            for j = [2 1]
                for l = [1 3 2]
                    fprintf(' %6.3f ',res(m,j,k,l));
                end
                fprintf('|');
            end
            fprintf('\n');
        end
    end

end