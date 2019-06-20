% Apply AJIVE with data from 2 datasets
clear
clc
addpath ../AJIVECode/

ioutput = [1, 1, 1, 1, 1, 1, 1, 1, 1];
paramstruct0 = struct('ioutput', ioutput, ...
                      'iplot', [0 0]);

% Do a loop over the generated data
% Case d = 2, no scale
folder = 'd2/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
error = zeros(100,1);

for i=1:100
    
    % Load the data
    filename = [folder 'rep' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = X';
    datablock{2} = Y'; 
    p1 = size(X,2);
    p2 = size(Y,2);
    
    vecr = vecr';

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);

    % Get back the estimation error
    signal = V * U';
    for d=1:2
        estimate = (outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d}) * sqrt(norms(d));
        if d==1
            index = 1:p1;
        else
            index = (p1+1):(p1+p2);
        end
        error(i) = error(i) + sum(sum((signal(index,:) - estimate).^2))/sum(sum(signal(index,:).^2));
    end    
    
end

savefile = [folder 'results.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');


% Case d = 2, scale
folder = 'd2_scale/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
error = zeros(100,1);

for i=1:100
    
    % Load the data
    filename = [folder 'rep' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = X';
    datablock{2} = Y'; 
    p1 = size(X,2);
    p2 = size(Y,2);
    
    vecr = vecr';

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);

    % Get back the estimation error
    signal = V * U';
    for d=1:2
        estimate = (outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d}) * sqrt(norms(d));
        if d==1
            index = 1:p1;
        else
            index = (p1+1):(p1+p2);
        end
        error(i) = error(i) + sum(sum((signal(index,:) - estimate).^2))/sum(sum(signal(index,:).^2));
    end    
    
end

savefile = [folder 'results.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');

% Case d = 2, varying p
folder = 'd2_p/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
error = zeros(100,1);

for i=1:100
    
    % Load the data
    filename = [folder 'rep' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = X';
    datablock{2} = Y'; 
    p1 = size(X,2);
    p2 = size(Y,2);
    
    vecr = vecr';

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);

    % Get back the estimation error
    signal = V * U';
    for d=1:2
        estimate = (outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d}) * sqrt(norms(d));
        if d==1
            index = 1:p1;
        else
            index = (p1+1):(p1+p2);
        end
        error(i) = error(i) + sum(sum((signal(index,:) - estimate).^2))/sum(sum(signal(index,:).^2));
    end    
    
end

savefile = [folder 'results.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');

% Case d = 2, AC model
folder = 'd2_AC/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
error = zeros(100,1);

for i=1:100
    
    % Load the data
    filename = [folder 'rep' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = X';
    datablock{2} = Y'; 
    p1 = size(X,2);
    p2 = size(Y,2);
    
    vecr = vecr';

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);

    % Get back the estimation error
    signal = V * U';
    for d=1:2
        estimate = (outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d}) * sqrt(norms(d));
        if d==1
            index = 1:p1;
        else
            index = (p1+1):(p1+p2);
        end
        error(i) = error(i) + sum(sum((signal(index,:) - estimate).^2))/sum(sum(signal(index,:).^2));
    end    
    
end

savefile = [folder 'results.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');

% Case d = 3
folder = 'd3/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);

for i=1:100
    
    % Load the data
    filename = [folder 'rep' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = X';
    datablock{2} = Y'; 
    datablock{3} = Z'; 
    p1 = size(X,2);
    p2 = size(Y,2);
    p3 = size(Z,2);
    
    vecr = vecr';

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    signal = V * U';
    for d=1:3
        estimate = (outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d}) * sqrt(norms(d));
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        error(i) = error(i) + sum(sum((signal(index,:) - estimate).^2))/sum(sum(signal(index,:).^2));
    end    
    
end

savefile = [folder 'results.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec','r3vec','error');

% Case d = 3, noise
folder = 'd3_noise/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);

for i=1:100
    
    % Load the data
    filename = [folder 'rep' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = X';
    datablock{2} = Y'; 
    datablock{3} = Z'; 
    p1 = size(X,2);
    p2 = size(Y,2);
    p3 = size(Z,2);
    
    vecr = vecr';

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    signal = V * U';
    for d=1:3
        estimate = (outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d}) * sqrt(norms(d));
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        if (sum(sum(signal(index,:).^2))>0)
            error(i) = error(i) + sum(sum((signal(index,:) - estimate).^2))/sum(sum(signal(index,:).^2));
        else
            error(i) = error(i) + sum(sum(estimate.^2));
        end    
    end    
    
end

savefile = [folder 'results.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec','r3vec','error');

% Case d = 3, perturb
folder = 'd3_perturb/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);

for i=1:100
    
    % Load the data
    filename = [folder 'rep' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = X';
    datablock{2} = Y'; 
    datablock{3} = Z'; 
    p1 = size(X,2);
    p2 = size(Y,2);
    p3 = size(Z,2);
    
    vecr = vecr';

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    signal = V * U';
    for d=1:3
        estimate = (outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d}) * sqrt(norms(d));
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        if (sum(sum(signal(index,:).^2))>0)
            error(i) = error(i) + sum(sum((signal(index,:) - estimate).^2))/sum(sum(signal(index,:).^2));
        else
            error(i) = error(i) + sum(sum(estimate.^2));
        end    
    end    
    
end

savefile = [folder 'results.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec','r3vec','error');




