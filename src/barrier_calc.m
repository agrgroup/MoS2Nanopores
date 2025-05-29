function [Ea, success]=barrier_calc(atomtype)
% Initial values (dummy)
success = 0;
Ea = 1000;

% Energy of E_S or E_Mo (under eqm)
E_MoS2 = -23.0310123; 
% E_S = -5.678026435; % Mo rich

E_S_bulk = -4.32257173;
inpt = -0.8;
E_S = inpt + E_S_bulk;

E_Mo = E_MoS2 - 2*E_S;

% Intrinsic barrier
Ea0 = 1.077897;

% Case catalog
Cases_Mo = ['2433330'
'2433220'
'2533233'
'2333200'
'2233000'
'2533133'
'2522332'
'2522331'
'2422220'
'2433120'
'2422120'
'2433110'
'2522221'
'2422110'
'2333100'
'2322200'
'2321100'
'2222000'];

Cases_Mo_sorted = zeros(size(Cases_Mo));
for i = 1:length(Cases_Mo)
    Cases_Mo_sorted(i) = sort_fingerprint(Cases_Mo(i,:)); % After removing 2 from case name in begining
end

Cases_S = ['12660000'
'32550001'
'32350001'
'31500001'
'12460000'
'12650000'
'32540001'
'32540000'
'11400001'
'32550000'
'32440001'
'12440000'
'31300001'
'11600000'
'32430001'
'12360000'
'32420001'];

Cases_S_sorted = zeros(size(Cases_S));
for i = 1:length(Cases_S)
    Cases_S_sorted(i) = sort_fingerprint(Cases_S(i,:)); % After removing 1 or 3 from the case name in begining
end

% Initial Energies
E_init_Mo = [-1774.856644
-1762.706458
-1833.120478
-1755.468411
-1849.224751
-1843.377574
-1832.976746
-1821.080126
-1825.549370
-1739.109971
-1712.658595
-1714.757354
-1795.355431
-1692.170718
-1731.612020
-1820.689306
-1510.854845
-1815.221313];

E_init_S = [-1838.849668
-1832.818334
-1757.106765
-1746.811353
-1762.683011
-1832.976497
-1828.976408
-1832.541659
-1741.618088
-1814.722620
-1809.307592
-1749.334059
-1812.779254
-1750.736066
-1751.392015
-1821.126714
-1696.276958];

% Final Energies
E_final_Mo = [-1763.925307
-1751.214241
-1820.923430
-1746.824617
-1840.268139
-1834.309026
-1821.704054
-1809.407860
-1814.680537
-1727.797076
-1704.296495
-1703.268170
-1786.732310
-1679.605537
-1723.397850
-1810.986483
-1499.346692
-1805.118023];

E_final_S = [-1833.705007
-1827.599721
-1752.494970
-1740.484187
-1756.916882
-1829.205485
-1822.830996
-1828.356545
-1735.100784
-1810.160124
-1803.695788
-1744.331278
-1805.101711
-1746.811386
-1745.097012
-1817.163093
-1689.815726];

% Calculate barrier
first = str2num(atomtype(1));
sorted_atmtyp = sort_fingerprint(atomtype);

if first == 2
    idx = find(Cases_Mo_sorted == sorted_atmtyp);
    if (isempty(idx))
        fprintf('Mo atomtype not found.\n');
    else
        del_E = E_final_Mo(idx) + E_Mo - E_init_Mo(idx);
        Ea = Ea0 * (1 + ((del_E)/(4*Ea0)))^2;
        success = 1;
    end
else
    idx = find(Cases_S_sorted == sorted_atmtyp);
    if (isempty(idx))
        if sorted_atmtyp == 0000001
            Ea = 0;
            success = 1;
        else
            fprintf('S atomtype not found.\n');
        end
        
    else
        del_E = E_final_S(idx) + E_S - E_init_S(idx);
        Ea = Ea0 * (1 + ((del_E)/(4*Ea0)))^2;
        success = 1;
    end
end

end

function sorted_fp = sort_fingerprint(fp)
first = str2num(fp(1));
copy = fp(2:end);
a = [];
for i = 1:length(copy)
    a=[a; str2num(copy(i))];
end

if first == 2
    a = sort(a);
else
    last = a(end);
    a = sort(a(1:end-1));
    a = [a; last];
end

sorted_fp = '';
% sorted_fp = strcat(sorted_fp,first);
for i = 1 : length(a)
    sorted_fp = strcat(sorted_fp, num2str(a(i)));
end
sorted_fp = str2num(sorted_fp);

end