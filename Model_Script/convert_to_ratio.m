function [rela_ratio, crel_ratio] = convert_to_ratio(output)

untran_r = 0.130;
untran_c = 0.035;

%RelA Calculations
rela_cyto = sum(output(:,1:5),2);
rela_nuc = sum(output(:,6:8),2);
rela_total = (3.5*rela_cyto+rela_nuc)/4.5;
rela_ratio = rela_nuc/(rela_total(1)+untran_r);
rela_ratio = rela_ratio - rela_ratio(1);

%cRel Calculations
crel_cyto = sum(output(:,9:13),2);
crel_nuc = sum(output(:,14:16),2);
crel_total = (3.5*crel_cyto+crel_nuc)/4.5;
crel_ratio = crel_nuc/(crel_total(1)+untran_c);
crel_ratio = crel_ratio - crel_ratio(1);
end