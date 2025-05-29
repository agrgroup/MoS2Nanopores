function [tf,tknock_timeseries] = kmc_isomers_SiEtch(isoindex,poresize,basedir)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global basedirectory nuperp N T R NA kb kb_eV area_site lattice_const num_dt isomer_index pore_size fid_reg;

% Initialize time
t0=0;

% Set isomer index, pore size, and base directory global variables
isomer_index = isoindex;
pore_size = poresize;
basedirectory = basedir;

% Set grid size
N = 18; % No of edge atoms (must be an even number)

% Number of KMC steps required to form pore of given size
num_dt = pore_size-8;

% Values of physical constants
R=8.314;        % Gas constant in J/mol-K
NA=6.022e23;    % Avogadro's number in mol^-1
kb=R/NA;        % Boltzmann constant in SI units
kb_eV=8.6173324e-5;% Boltzmann constant in eV/K

% Set attempt frequency for A and B sublattices
nuperp = 1e13*[1,1,1];

% Set temperature of system in Kelvin
T = 300.15;

% Set lattice size of MoS2
lattice_const = 2*2.39959*cosd(49.0278);
area_site = 0.50*sqrt(3)*(lattice_const^2)*1e-20;
      
% Initialize 2D material layer with a monoatom vacancy
[~,C,num_atoms]=init(N);
 
% Do the KMC algorithm
tic
[tf,tknock_timeseries]=kmc(C,num_atoms,t0,N);
toc

end

%% Function to perform the kmc
function [t,tknock_timeseries]=kmc(C0,num_atoms,t0,N)
% This function performs the KMC algorithm and returns final configuration
% of the system

% Inputs: initial surface configuration (C0), number of atoms on A and B sublattices (num_atoms), initial time (t0), size of
% grid (N)

% Outputs: tf - formation time of nanopore, tknock_timeseries - time
% required for atom to be knocked out at each timestep,

global num_dt;
global basedirectory isomer_index pore_size; % Only for printing mechanism

t=t0; % initial time
count=0; % counter for number of time steps
C=C0; % store initial state of system in variable C

% Determine the list of sites which can be etched
sites_list=sites(C,N);
r = sites_list(:,7); % vector of all rates
num=length(r);

% Initialize variables to store timeseries
tknock_timeseries = zeros(num_dt+1,1);
% movieviz(C,N,count);

% % Print
% fid_mech = fopen([basedirectory,'/pore',num2str(pore_size),...
%            '/pore',num2str(isomer_index),'_mech.txt'], 'wt');
while (num>0 && count<=num_dt)
        
    % Generate two uniformly-distributed random numbers u1 and u2
    random = rand(1,2);
    u1=random(1);
    u2=random(2);
      
    % Compute effective rate vector
    [num,~]=size(sites_list);
    rcum = cumsum(r);
    rcum=rcum/rcum(num); % normalize with total rate
    
    tknock_timeseries(count+1) = 1/sum(r);

    % Find which of the sites is going to have a reaction, based on
    % Gillespie's algorithm
    site_sel=0;
    for i=1:num       
        if (i==1 && u1<=rcum(i))
            site_sel=i;
        elseif (i>1 && u1>rcum(i-1) && u1<=rcum(i))
            site_sel=i;
        end     
    end
    
    % Increment the system time based on a Poisson process
    rnet=sum(r);
    if (rnet>0)
        dt = log(1/u2)/rnet;
        t=t+dt;   
        count=count+1;
    end

    if (site_sel>0)
        % Find (r,c) position of the lattice site selected to be etched
        pos_sel_r = sites_list(site_sel,4); 
        pos_sel_c = sites_list(site_sel,5);
        
        % Find sublattice at which etching is going to occurr
        sublattice_sel = sites_list(site_sel,1);
        
        % Find if selected lattice site is already etched
        site_status = sites_list(site_sel,6);
        
        if (site_status==0)
        %    % Printing mechanism of removal - comment out if not needed
        %    atom_sig = atomtype(pos_sel_r, pos_sel_c, sublattice_sel, C);
        %    fprintf(fid_mech,"%d \t %s \t %3.6f\n",count,atom_sig,dt);

           
           % If selected lattice site is 'not etched', then update site
           % status to 'etched'
           C(pos_sel_r,pos_sel_c,sublattice_sel) = sublattice_sel;
           num_atoms(sublattice_sel) = num_atoms(sublattice_sel)+1; 
        end           
        
        % Update list of available sites
        [sites_list, ~]=sites_change(C,sublattice_sel,[pos_sel_r,pos_sel_c],sites_list);
        
        % Compute the next step's rate vector
        r = sites_list(:,7);                    
    end
      
          
    if (mod(count,1000)==0)
        disp(['count = ', num2str(count)]);
    end
end
% % Close print file
% fclose(fid_mech);

% Extract and save the knocking time at
% each time step
tknock_timeseries(end+1) = 1/sum(r);

% Save XYZ file and antimolecule adjacency matrix at the last timestep
visualize(C,N,count+1);  

end

%% Function to classify the atomtype
function y=atomtype(i,j,sublattice,C)
% This function classifies a given lattice site to be an armchair, zigzag,
% singly-bonded, lone, or fully-bonded carbon atom

% Inputs: location (i,j) of lattice site, sublattice, state of the system (C)
% Output: atomtype which is a sequence of numbers 
%         (sublattice-N-NN(1)-NN(2)-NN(3)-NN(4)-NN(5)-bottom/top sulphur removed if applicable)
%         Where N=Nearest, NN=Next Nearest Array positions

% Initialize switch variable
check='';

check=strcat(check,int2str(sublattice));

% Store the row and the column position
rowpos = i;
colpos = j;

% Find the list of nearest neighbor lattice sites to the given lattice site
[nearest,~]=neighbor_sites(i,j,sublattice);

% Initialize number of vacancy neighbors, number of filled neighbors, and
% number of filled neighbors of filled neighbors, all to be zero
num_vacancy_neigh = 0;
num_atomic_neigh  = 0;
num_atomic_neigh_atomic_neigh = zeros(5,1);

% Loop through all nearest neighbor lattice sites
for i=1:size(nearest,1)
    if (C(nearest(i,1),nearest(i,2),nearest(i,3))==0) % nearest site is not a vacancy, so it is filled with an atom
        num_atomic_neigh = num_atomic_neigh+1;
        how_many_vacancy_neigh_atomic_neigh = count_vacancy_neighbors(nearest(i,1),nearest(i,2),nearest(i,3),C); % this function counts the number of vacancy neighbors
        % Adjustment of num_max_neigh due to sublattice
        if (nearest(i,3)==2) % If sublattice = 2 
            num_max_neigh = 6;
        else
            num_max_neigh = 3;
        end
        num_atomic_neigh_atomic_neigh(num_atomic_neigh) = num_max_neigh-(how_many_vacancy_neigh_atomic_neigh(1));
    else
        num_vacancy_neigh = num_vacancy_neigh+1;
    end
end
check = strcat(check,int2str(num_atomic_neigh));
for i = 1:size(num_atomic_neigh_atomic_neigh,1)
    check = strcat(check,int2str(num_atomic_neigh_atomic_neigh(i)));
end

% Accounting for top or bottom sulphur removed
if (sublattice == 1)
    if (C(rowpos,colpos,3)==0) % Bottom Sulphur is present
        check = strcat(check,int2str(0));
    elseif (C(rowpos,colpos,3)==3) % Bottom Sulphur is removed
        check = strcat(check,int2str(1));
    end
elseif (sublattice == 3)
    if (C(rowpos,colpos,1)==0) % Top Sulphur is present
        check = strcat(check,int2str(0));
    elseif (C(rowpos,colpos,1)==1) % Top Sulphur is removed
        check = strcat(check,int2str(1));
    end
end

y = check; % Storing the id of the atom type for processing.

end

%% Function to specify the etching barrier
function y=barrier(atomtype,C,pos_r,pos_c,sublattice)
% This function returns the etching barrier for a given atomtype
% Input: atomtype
% Output: etching barrier (in eV)

global N pore_size lattice_const

[y, success] = barrier_calc(atomtype);

if success == 0
    disp(atomtype)
    disp('Barrier not found.')
    disp(pos_r)
    disp(pos_c)
    disp(C)
    Visualize(C,N,pore_size)
    i=pos_r;
    j=pos_c;
    if (C(i,j,1)==0 && sublattice==1)
        % atom is present, i.e. vacancy not added
        coord_sublattice_1 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
        disp(['S ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   1.568']);           

    end

    if (C(i,j,2)==0 && sublattice==2)
        % atom is present, i.e. vacancy not added
        coord_sublattice_2 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];
        disp(['Mo ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000']);

    end

    if (C(i,j,3)==0 && sublattice==3)
        % atom is present, i.e. vacancy not added
        coord_sublattice_1 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
        fprintf(['S ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   -1.568']);           

    end
end

end

%% Function to specify the various sites
function sites_list=sites(C,N)
% Compute the positions of all possible relavent interactions
% The columns are 1. sublattice, 2. number of vacancy neighbors, 
% 3. number of vacancy next neighbors, 4. row-position (y),
% 5. column-position (x), 6. vacancy flag (1-occupied by vacancy,
% 0-occupied by atom), 7. etching rate

global nuperp T kb_eV;

% Initialize monoatom vacancy in lattice, if need other shape, change
% init() function
[initial_vacancy,~,~]=init(N);

[num,~]=size(initial_vacancy);
sites_list=zeros(num,7);
sites_list(:,[1,4,5])=initial_vacancy(:,[1,2,3]); % store sublattice, row location, column location in 1st, 4th, and 5th columns of sites_list

for i=1:num 
    % This has been made a loop so that if one wants to initialize with a
    % defect other than a monoatom vacancy, the code can be easily modified simply by changing `initial_vacancy`
    sublattice = initial_vacancy(i,1);
    pos = initial_vacancy(i,2:3);

    [nearest,next_nearest]=neighbor_sites(pos(1),pos(2),sublattice);
    [n_near, ~] = size(nearest);
    [n_nextnear, ~] = size(next_nearest);

    % Update vacancy flag (i.e., 6th column) to be 1
    [~,locb] =   ismember([sublattice,pos(1),pos(2)],sites_list(:,[1,4,5]),'rows');
    sites_list(locb,6) = 1;
    
    % Loop through nearest neighbor lattice sites
    for k=1:n_near
        posr = nearest(k,1);
        posc = nearest(k,2);
        neigh_sublattice = nearest(k,3);
        
        [lia,locb] =   ismember([neigh_sublattice,posr,posc],sites_list(:,[1,4,5]),'rows');      
        if (~any(lia)) % If lattice site is not in sites_list, add it there
            sites_list = [sites_list; neigh_sublattice,1,0,posr,posc,ismember([neigh_sublattice,posr,posc],initial_vacancy,'rows'),0];
        
        else % If lattice site is already there, update number of vacancy neighbors
            sites_list(locb,2) = sites_list(locb,2)+1;
        end
    end
    
    for k=1:n_nextnear
        posr = next_nearest(k,1);
        posc = next_nearest(k,2);
        neigh_sublattice = next_nearest(k,3);
            
        [lia,locb] =   ismember([neigh_sublattice,posr,posc],sites_list(:,[1,4,5]),'rows');
        if (~any(lia)) % If lattice site is not in sites_list, add it there
            sites_list = [sites_list; neigh_sublattice,0,1,posr,posc,ismember([neigh_sublattice,posr,posc],initial_vacancy,'rows'),0];
        else % If lattice site is already there, update number of vacancy next neighbors
            sites_list(locb,3) = sites_list(locb,3)+1;
        end
    end
end   

% Loop through all lattice sites in the sites list
for i=1:size(sites_list,1)
    sublattice = sites_list(i,1);
    num_vacancy_neigh = sites_list(i,2);
    vacancy_flag = sites_list(i,6);
    posr = sites_list(i,4);
    posc = sites_list(i,5);
    
    if (vacancy_flag==0) % site doesn't have vacancy
        if (num_vacancy_neigh>0) % site has at least one vacancy neighbor
            Ea = barrier(atomtype(posr,posc,sublattice,C));
            sites_list(i,7)=nuperp(sublattice)*exp(-Ea/(kb_eV*T));  % Rate constant for etching of atom
        else
            sites_list(i,7)=0; % site has no vacancy neighbor, so cannot be etched
        end
    else  % site already has vacancy, so etching rate is zero for lattice site
        sites_list(i,7) = 0;
    end
end
end

%% Function to implement sites update
function [updated_sites_list,modified_sites_index]=sites_change(C,sublattice,pos,sites_list)
% This function updates the sites_list for the simulation and returns an
% updated list of potential sites for etching

% Inputs: state of the system (C), properties of currently etched site
% (sublattice, pos), current sites list (sites_list)
% Outputs: updated sites list and list of modified lattice sites

global nuperp T kb_eV;

updated_sites_list = sites_list;

% Update properties of etched lattice site
[~,locb] =   ismember([sublattice,pos(1),pos(2)],updated_sites_list(:,[1,4,5]),'rows'); % find etched site in sites_list
original_state = updated_sites_list(locb,6);
updated_sites_list(locb,6) = ~original_state;
modified_sites_index = locb;

[nearest,next_nearest]=neighbor_sites(pos(1),pos(2),sublattice);
[n_near, ~] = size(nearest);
[n_nextnear, ~] = size(next_nearest);

found_neighbor = 0;

% Loop over nearest neighbors of the etched site
for k=1:n_near
    posr = nearest(k,1);
    posc = nearest(k,2);
    neigh_sublattice = nearest(k,3);

    % Try to find neighoring site in updated_sites_list
    [lia,locb] =   ismember([neigh_sublattice,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
    
    if (~any(lia)) % if not found, add it to the list
        updated_sites_list = [updated_sites_list; neigh_sublattice,1,0,posr,posc,C(posr,posc,neigh_sublattice)==neigh_sublattice,0];
        [num_sites,~] = size(updated_sites_list);
        modified_sites_index = [modified_sites_index; num_sites];
    
    else % if found, update the number of vacancy neighbors to the site
        updated_sites_list(locb,2) = updated_sites_list(locb,2)+1;
        modified_sites_index = [modified_sites_index; locb];
    end
end

% Loop over next nearest neighbors of the etched site
for k=1:n_nextnear
    posr = next_nearest(k,1);
    posc = next_nearest(k,2);
    neigh_sublattice = next_nearest(k,3);

    % Try to find next neighboring site in updated_sites_list
    [lia,locb] =   ismember([neigh_sublattice,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
    
    if (~any(lia)) % if not found, add it to the list
        updated_sites_list = [updated_sites_list; neigh_sublattice,0,1,posr,posc,C(posr,posc,neigh_sublattice)==neigh_sublattice,0];
        [num_sites,~] = size(updated_sites_list);
        modified_sites_index = [modified_sites_index; num_sites];
    
    else % if found, update the number of next nearest vacancy neighbors to the site
        updated_sites_list(locb,3) = updated_sites_list(locb,3)+1;
        modified_sites_index = [modified_sites_index; locb];
    end
end

% Loop over updated_sites_list to perform housekeeping check
for i=1:size(updated_sites_list,1)
    sublattice = updated_sites_list(i,1);
    num_neigh = updated_sites_list(i,2);
    num_next_neigh = updated_sites_list(i,3);
    vacancy_flag = updated_sites_list(i,6);
    posr = updated_sites_list(i,4);
    posc = updated_sites_list(i,5);
    
    if (~all(count_vacancy_neighbors(posr,posc,sublattice,C)==[num_neigh;num_next_neigh]))
        display(posr)
        display(posc)
        display(vacancy_flag)
        stored_neighs=[num_neigh,num_next_neigh];
        actual_neighs=count_vacancy_neighbors(posr,posc,sublattice,C)';
        display(sublattice)
        display(stored_neighs)
        display(actual_neighs)
        error('Neighbor coordination number does not match');
    end
end

% Loop over modified sites to assign them rates, and to remove those modified sites without nearest
% neighbor or next-nearest neighbor vacancy sites
for_removal=[];
for i=1:length(modified_sites_index)
    index = modified_sites_index(i);
    sublattice = updated_sites_list(index,1);
    num_neigh = updated_sites_list(index,2);
    num_next_neigh = updated_sites_list(index,3);
    vacancy_flag = updated_sites_list(index,6);
    pos_r = updated_sites_list(index,4);
    pos_c = updated_sites_list(index,5);
    
    if (vacancy_flag==0) % site doesn't have vacancy yet
        if (num_neigh==0 && num_next_neigh==0)
            for_removal = [for_removal;index];
        else
            if (num_neigh>0) % needs to have atleast one vacancy as neighbor
                Ea = barrier(atomtype(pos_r,pos_c,sublattice,C),C,pos_r,pos_c,sublattice);
                updated_sites_list(index,7)=nuperp(sublattice)*exp(-Ea/(kb_eV*T));  % Rate for etching
            else
                updated_sites_list(index,7)=0;
            end
        end
    else % site already has vacancy, so etching rate is set to zero
            updated_sites_list(index,7) = 0;
    end
end
updated_sites_list(for_removal,:) = [];


if (found_neighbor==0 && original_state==1)
    disp('Exception: Site with no neighbor / next-nearest neighbor removed');
end
end

%% Function to implement the Periodic Boundary Condition(PBC)
function y=pbc(x)
% This function implements PBC for lattice indices.
% Input: row/column before wrapping
% Output: row/column after periodic boundary wrapping
global N;
if (x==0 || x==N)
    y=N;
else
    y=mod(x,N);
end
end

%% Function to classify neighbour sites
function [nearest,next_nearest]=neighbor_sites(i,j,sublattice)
% This function returns the nearest and next-nearest lattice sites for a
% given site (i, j, sublattice)

% Inputs: row (i), column (j), and sublattice (1 or 2)
% Outputs: list of nearest neighbor and next-nearest neighbor lattice sites
% in the format (row, column, sublattice)

if (sublattice==1 || sublattice==3)
    
    if (sublattice==1)
        const=3;
    else
        const=1;
    end
    nearest = [pbc(i+1),    pbc(j-1+(~mod(i,2))),        2;
               pbc(i+1),    pbc(j+(~mod(i,2))),          2;
               i,           j,                           2];% if even row add 1, odd row add 0
    
    next_nearest = [pbc(i+1),    pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                pbc(j+(~mod(i,2))),          sublattice;
        i,                       pbc(j-1),                    sublattice;
        i,                       pbc(j+1),                    sublattice;
        pbc(i-1),                pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                pbc(j+(~mod(i,2))),          sublattice;
        i,                       j,                                const;
        pbc(i+1),                pbc(j-1+(~mod(i,2))),             const;
        pbc(i+1),                pbc(j+(~mod(i,2))),               const;
        i,                       pbc(j-1),                         const;
        i,                       pbc(j+1),                         const;
        pbc(i-1),                pbc(j-1+(~mod(i,2))),             const;
        pbc(i-1),                pbc(j+(~mod(i,2))),               const];
    
elseif (sublattice==2)

    nearest = [i,       j,                       1;
        i,              j,                       3;
        pbc(i-1),       pbc(j+(~mod(i,2))),      1;
        pbc(i-1),       pbc(j+(~mod(i,2))),      3;
        pbc(i-1),       pbc(j-1+(~mod(i,2))),    1;
        pbc(i-1),       pbc(j-1+(~mod(i,2))),    3];
    
    next_nearest = [pbc(i+1),       pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                   pbc(j+(~mod(i,2))),          sublattice;
        i,                          pbc(j-1),                    sublattice;
        i,                          pbc(j+1),                    sublattice;
        pbc(i-1),                   pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                   pbc(j+(~mod(i,2))),          sublattice];
end
end

%% Function to counting vacant neighbor sites
function y=count_vacancy_neighbors(i,j,sublattice,C)
% This function counts the number of vacancy neighbors to a given lattice
% site

% Input: row (i), column (j), sublattice (1,2 or 3), state of the system (C)
% Output: number of vacancy neighbors in the format (nearest, next-nearest)

[nearest,next_nearest]=neighbor_sites(i,j,sublattice);
[n_near, ~] = size(nearest);
[n_nextnear, ~] = size(next_nearest);

y=[0;0];
for k=1:n_near
    posr = nearest(k,1); % row number
    posc = nearest(k,2); % col position
    neigh_sublattice = nearest(k,3);
    % checking if the position has a vacancy. 1 if vacancy on sublattice 1
    % 2 if vacancy on sublattice 2.
    if (C(posr,posc,neigh_sublattice)==neigh_sublattice) % site contains vacancy
        y(1)=y(1)+1;
    end
end
for k=1:n_nextnear
    posr = next_nearest(k,1);
    posc = next_nearest(k,2);
    next_neigh_sublattice = next_nearest(k,3);
    if (C(posr,posc,next_neigh_sublattice)==next_neigh_sublattice) % site contains vacancy
        y(2)=y(2)+1;
    end
end
end

%% Function to initialise the system
function [initial_vacancy,C0,num_atoms]=init(N)
% This function initialises the state of the system.

% Input: The size of the lattice(N), Number of lattice site is N*N*3, the
% factor of 3 is due to there being three sublattices in MoS2
% Output: Initial state of the system (initial_vacancy and MS0), number of
% vacancies in each sublattice (num_atoms)

% Initialsising the system with a Mo monoatomic vacancy, can be tweaked as
% and when required. 0 signifies filled with an atom, 1 signifies vacancy.
C0 = zeros(N,N,3);
center = round(N/2);

% % S vac start
% C0(center,center,1) = 1;

% % Mo vac start
% C0(center,center,2) = 2;

% % MoS2 vac start
% C0(center,center,2) = 2;
% 
% C0(center,center,1) = 1;
% C0(center,center,3) = 3;

% MoS6 vac start
C0(center,center,2) = 2;

C0(center,center,1) = 1;
C0(center,center,3) = 3;

C0(center-1,center,1) = 1;
C0(center-1,center,3) = 3;

C0(center-1,center-1,1) = 1;
C0(center-1,center-1,3) = 3;
% C0;

initial_vacancy = [1 center center;3 center center;2 center center;1 center-1 center;3 center-1 center;1 center-1 center-1;3 center-1 center-1];

num_atoms = [3,1,3];
end

%% Function to visualize the results
function []=visualize(C,N,count)
% This function saves the current state of the system, given the surface
% state matrix, as an XYZ file and a directed adjacency matrix for the
% antimolecule of the formed nanopore

% Inputs: Surface state matrix (C), Size of grid (N), and current timestep (count)
% Outputs: number of armchair, zigzag, and unassigned atoms, and 5-membered
% rings

global lattice_const isomer_index pore_size basedirectory;

% Define box sizes in the x and y directions
box_x = N*2*(lattice_const/sqrt(3))*cos(pi/6);
box_y = N*(lattice_const/sqrt(3))*(2+2*cos(pi/3))/2;

% if (count<pore_size && count>=N*3) % What is this condition doing?
if (count<(pore_size-6) && count>=7) % What is this condition doing?
    fprintf("The count is %d and pore size is %d\n",count,pore_size);
else
    fid_xyz = fopen(      [basedirectory,'/pore',num2str(pore_size),'/pore',            num2str(isomer_index)  ,'.xyz'],'wt');
    fid1_xyz = fopen( [basedirectory,'/pore',num2str(pore_size),'/antimolecule',num2str(isomer_index),'.xyz'],'wt');
    fid_txt = fopen(      [basedirectory,'/pore',num2str(pore_size),'/pore',            num2str(isomer_index)  ,'.txt'],'wt');
    fid_adjmat_removed =  [basedirectory,'/pore',num2str(pore_size),'/adjmat_antimol_', num2str(isomer_index)  ,'.txt'];
end

% Create variables for number of total atoms, number of remaining atoms,
% and number of removed atoms
[m,n]=size(C(:,:,1));
num_total_atoms = m*n*3;
num_remaining_atoms = sum(sum(C(:,:,1)==0))+sum(sum(C(:,:,2)==0))+sum(sum(C(:,:,3)==0));
num_removed_atoms = m*n*3-num_remaining_atoms;

count_removed_atoms=0;
% count_rim_atoms=0;
no_of_atoms = N*N*3 - pore_size;
fprintf(fid_xyz,[num2str(no_of_atoms), '   \n']);
fprintf(fid_xyz,['Mo S2   \n']);

% Create file to write antimolecules
fprintf(fid1_xyz,[num2str(num_removed_atoms), '   \n']);
fprintf(fid1_xyz,['Mo S2   \n']);

% Cycle through all locations on the 2D lattice
for i=1:N
    for j=1:N       
        if (C(i,j,1)==0)
            % atom is present, i.e. vacancy not added
            coord_sublattice_1 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            fprintf(fid_xyz,['S ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   1.568\n']);           
            
        else
            % print removed atoms with different name
            coord_sublattice_1 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];           
            fprintf(fid1_xyz,['S ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   1.568\n']);
          
            % assign index to removed atoms
            fprintf(fid_txt,'%d %d %d\n',1,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,1) = count_removed_atoms;
        end
        
        if (C(i,j,2)==0)
            % atom is present, i.e. vacancy not added
            coord_sublattice_2 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];
            fprintf(fid_xyz,['Mo ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000\n']);
            
        else
            % print removed atoms with different name
            coord_sublattice_2 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];
            fprintf(fid1_xyz,['Mo ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000\n']);
            
            % assign index to removed atoms
            fprintf(fid_txt,'%d %d %d\n',2,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,2) = count_removed_atoms;
        end
        
        if (C(i,j,3)==0)
            % atom is present, i.e. vacancy not added
            coord_sublattice_1 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            fprintf(fid_xyz,['S ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   -1.568\n']);           
            
        else
            % print removed atoms with different name
            coord_sublattice_1 = [ (~mod(i,2)*(lattice_const/2)) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];           
            fprintf(fid1_xyz,['S ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   -1.568\n']);
          
            % assign index to removed atoms
            fprintf(fid_txt,'%d %d %d\n',3,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,3) = count_removed_atoms;
        end
    end
end



% make adjacency matrix out of removed atoms, including bond orientations
adjmat_removed = zeros(num_removed_atoms, num_removed_atoms);
adj_direction = 0;
for i=1:size(C,1)
    for j=1:size(C,2)
        curr_index_1 = C(i,j,1);
        if (curr_index_1>0)  % current atom is removed atom or rim atom
            [nearest,~]=neighbor_sites(i,j,1);
            [n_near, ~] = size(nearest);
            coord_A = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];

            for k=1:n_near
                posr = nearest(k,1);
                posc = nearest(k,2);
                neigh_spec = nearest(k,3);
                coord_B = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2   ];

                direction = coord_B-coord_A; % find direction of bond
                
                if (direction(1)>box_x/2.0) % right
                    direction(1) = -(box_x-direction(1));
                end
                if (direction(1)<-box_x/2.0) % left
                    direction(1) = (box_x+direction(1));                   
                end
                if (direction(2)>box_y/2.0) % top
                    direction(2) = -(box_y-direction(2));
                end
                if (direction(2)<-box_y/2.0) % down
                    direction(2) = (box_y+direction(2));
                end                  
                
                slope_direction = atan(abs(direction(2)/direction(1)))*180/pi;
                
                if (slope_direction > 85 && slope_direction < 95)
                    if (direction(2)>0)
                        adj_direction = 1; % top
                    else
                        adj_direction = 2; % bottom
                    end
                elseif (slope_direction>25 && slope_direction<35)
                    if (direction(2)>0) % top
                        if (direction(1)<0)
                            adj_direction = 3; % top left
                        else
                            adj_direction = 4; % top right
                        end
                    else                % bottom
                        if (direction(1)<0)
                            adj_direction = 5; % bottom left
                        else
                            adj_direction = 6; % bottom right
                        end
                    end
                end
                                   
                neigh_index = C(posr,posc,neigh_spec); 
                if (neigh_index>0) 
                   if (curr_index_1 <= num_removed_atoms && neigh_index <= num_removed_atoms)  % both current and neighboring atoms are removed atom only
                        adjmat_removed(curr_index_1,neigh_index)=adj_direction;
                   end
                end
            end
        end
        
        curr_index_2 = C(i,j,2);
        if (curr_index_2>0) 
            [nearest,~]=neighbor_sites(i,j,2);
            [n_near, ~] = size(nearest);
            coord_B = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];

            for k=1:n_near
                posr = nearest(k,1);
                posc = nearest(k,2);
                neigh_spec = nearest(k,3);
                coord_A = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];

                direction = coord_A-coord_B;
                
                if (direction(1)>box_x/2.0) % right
                    direction(1) = -(box_x-direction(1));
                end
                if (direction(1)<-box_x/2.0) % left
                    direction(1) = (box_x+direction(1));                   
                end
                if (direction(2)>box_y/2.0) % top
                    direction(2) = -(box_y-direction(2));
                end
                if (direction(2)<-box_y/2.0) % down
                    direction(2) = (box_y+direction(2));
                end                 
                                                 
                slope_direction = atan(abs(direction(2)/direction(1)))*180/pi;
                
                if (slope_direction > 85 && slope_direction < 95)
                    if (direction(2)>0)
                        adj_direction = 1; % top
                    else
                        adj_direction = 2; % bottom
                    end
                elseif (slope_direction>25 && slope_direction<35)
                    
                    if (direction(2)>0) % top
                        if (direction(1)<0)
                            adj_direction = 3; % top left
                        else
                            adj_direction = 4; % top right
                        end
                    else                % bottom
                        if (direction(1)<0)
                            adj_direction = 5; % bottom left
                        else
                            adj_direction = 6; % bottom right
                        end
                    end
                end
                                
                neigh_index = C(posr,posc,neigh_spec); 
                if (neigh_index>0)   % neighboring atom is removed or rim atom
                    if (curr_index_2 <= num_removed_atoms && neigh_index <= num_removed_atoms)  % both current and neighboring atoms are removed atom only
                        adjmat_removed(curr_index_2,neigh_index)=adj_direction;
                    end
                end
            end
        end
        
        curr_index_3 = C(i,j,3);
        if (curr_index_3>0)  % current atom is removed atom or rim atom
            [nearest,~]=neighbor_sites(i,j,3);
            [n_near, ~] = size(nearest);
            coord_A = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];

            for k=1:n_near
                posr = nearest(k,1);
                posc = nearest(k,2);
                neigh_spec = nearest(k,3);
                coord_B = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2   ];

                direction = coord_B-coord_A; % find direction of bond
                
                if (direction(1)>box_x/2.0) % right
                    direction(1) = -(box_x-direction(1));
                end
                if (direction(1)<-box_x/2.0) % left
                    direction(1) = (box_x+direction(1));                   
                end
                if (direction(2)>box_y/2.0) % top
                    direction(2) = -(box_y-direction(2));
                end
                if (direction(2)<-box_y/2.0) % down
                    direction(2) = (box_y+direction(2));
                end                  
                
                slope_direction = atan(abs(direction(2)/direction(1)))*180/pi;
                
                if (slope_direction > 85 && slope_direction < 95)
                    if (direction(2)>0)
                        adj_direction = 1; % top
                    else
                        adj_direction = 2; % bottom
                    end
                elseif (slope_direction>25 && slope_direction<35)
                    if (direction(2)>0) % top
                        if (direction(1)<0)
                            adj_direction = 3; % top left
                        else
                            adj_direction = 4; % top right
                        end
                    else                % bottom
                        if (direction(1)<0)
                            adj_direction = 5; % bottom left
                        else
                            adj_direction = 6; % bottom right
                        end
                    end
                end
                                   
                neigh_index = C(posr,posc,neigh_spec); 
                if (neigh_index>0) 
                   if (curr_index_3 <= num_removed_atoms && neigh_index <= num_removed_atoms)  % both current and neighboring atoms are removed atom only
                        adjmat_removed(curr_index_3,neigh_index)=adj_direction;
                   end
                end
            end
        end
    end
end
dlmwrite(fid_adjmat_removed, adjmat_removed);
fclose('all');
end