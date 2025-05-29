clear
clc
fclose('all');
% This program repeatedly calls the KMC routine and stochastically
% generates a new isomer each time. The properties of the isomer are stored
% in a CSV file

% Re-initialize the random number generator in MATLAB
rng('shuffle')

% Define a list of pore sizes for which to generate isomers
pore_size_list  = [22];

% Number of isomers to generate for each size
Niso=300;

% Directory in which the isomers will be stored. Each pore size N will have
% a subdirectory called poreN, e.g., pore6 and pore8 in the current version
% of the code
basedir = './catalog_ES_0.8/';

% Cycle through all pore sizes
for j=pore_size_list
    
    % Create directory porej
    dirname = [basedir,'pore',num2str(j)];
    status = mkdir(dirname);
    mkdir([dirname,'/path']);
    
    % Obtain properties of nanopore isomer generated
    tf_list = zeros(Niso,1);
    tknock_list = tf_list;
    
    % If directory was not created, it already exists, so some isomers have
    % already been generated
    if (status == 0)
       data = csvread([dirname,'/Analysis.csv']);
       num_isomers_done  = data(end,1);
    else
        
    % If directory was created, no isomers have yet been generated
       num_isomers_done = 0;
    end
    
    % Generate the remaining number of isomers
    for i=(num_isomers_done+1):Niso
        pore_index = j;
        iso_index = i;

        fprintf("\nPore index: %d\n",j);
        fprintf("Isomer index: %d\n",i);
        
        % Run the KMC algorithm
        [tf,tknock_timeseries]=kmc_isomers_SiEtch(i,j,basedir);
        
        % Store the isomer properties in lists
        tf_list(i) = tf;
        fprintf("The formation time: %3.2f\n",tf);
        tknock_list(i) = tknock_timeseries(end);
       
        % Write isomer properties in CSV file
        dlmwrite([dirname,'/Analysis.csv'],[i,tf_list(i),tknock_list(i)],'-append');
    end
end