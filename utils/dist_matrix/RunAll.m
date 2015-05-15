function [] = RunAll(hemi)
% RunAll('lh');

subjects = { '00796' '03820' '08747' '13081' '13753' '21923' '22529' '22872' '23197' '23269' '23278' '23302' '23305' '23353' '23426' '23428' '23429' '23464' '23514' '23574' '23576' '23602' '23629' '23649' '23650' '23651' '23668' '23688' '23697' '23700' '23705' '23708' '23713' '23734' '23741' '23780' '23850' '23860' '23861' '23884' '23886' '23946' '23947' '23948' '23960' '23963' '23965' '23983' '23985' '24002' '24010' '24024' '24050' '24057' '24058' '24061' '24102' '24134' '24246' '24247' '24275' '24444' '24445' '24587' '24614' '24691' '24700' '24703' '24707' '24715' '24719' '24720' '24730' '24731' '24732' '24756' '24757' '24759' '24766' '24772' '24773' '24774' '24790' '24843' '24844' '24877' '24899' '24915' '24916' '24918' '24944' '24945' '25019' '25036' '25081' '25165' '25171' '25177' '25188' '25189' '25192' '25194' '25195' '25197' '25199' '25200' '25201' '25260' '25264' '25274' '25283' '25326' '25524' '25552' '25862' '26091' '26120' '26121' '26139' '26141' '26142' '26265' '26266' '26271' '26284' '26297' '26311' '26320' '26331' '26336' '26350' '26361' '26362' '26372' '26407' '26409' '26410' '26413' '26419' '26435' '26436' '26449' '26450' '26471' '26489' '26493' '26496' '26498' '26499' '26500' '26524' '26526' '26571' '26589' '26591' '26617' '26618' '26619' '26626' '26635' '26642' '26659' '26674' '26687' '26704' '26719' '26723' '26724' '26727' '26728' '26729' '26739' '26753' '26768' '26776' '26782' '26789' '26793' '26795' '26801' '26802' '26803' '26804' '26805' '26806' '26820' '26824' '26839' '26841' '26842' '26843' '26844' '26854' '26855' '26856' '26857' '26858' '26902' '26917' '26926'};
% Problem: '24051' '24131'    

% for ReassembleDistMat: only select those that have complete distance files in ../temp/
if hemi == 'lh'
    meshlen = 9354;
elseif hemi == 'rh'
    meshlen = 9361;
end

% create a list with complete subjects and assemble distance matrix
complete = {};
for sub = 1:length(subjects)
    files = dir(['/scr/liberia1/data/lsd/surface/temp/*' subjects{sub} '_' hemi '*']);
    if length(files) == meshlen
        
        disp(['running ' subjects{sub}]);
        % DoDistIndividual2fsa5(subjects{sub}, hemi);
        ReassembleDistMat(subjects{sub}, hemi);
        complete{end+1} = ['mv /scr/liberia1/data/lsd/surface/temp/*' subjects{sub} '_' hemi '* /scr/liberia1/data/lsd/surface/distance/delete'];
    end
end

fid = fopen(['/scr/liberia1/data/lsd/surface/' hemi '_reassembled_subjects.txt'], 'w');
for x = 1 : length(complete)
    fprintf(fid,'%s \n', complete{x});
end
fclose(fid);