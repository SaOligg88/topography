function [surf] = prepSurf(surfType, hemi, sub, remove, outputPrefix)
% surfType = '32' or '164' or 'freesurfer'
% hemi = 'L' or 'R'
% remove = 'medial' or 'both'

switch surfType
    case '32'
        if hemi == 'L'
            surf = gifti('data/Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii');
            aparc = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/fsaverage_LR32k/' sub '.L.aparc.a2009s.32k_fs_LR.label.gii']);
        end
        if hemi == 'R'
            surf = gifti('data/Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii');
            aparc = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/fsaverage_LR32k/' sub '.R.aparc.a2009s.32k_fs_LR.label.gii']);
        end
        aparc = aparc.cdata;
        if strcmp(remove,'medial')
               noncortex = find(ismember(aparc, [0]));
        elseif strcmp(remove,'both')
               noncortex = find(ismember(aparc, [0 50 17 49 41 18])); % remove 35?
               %noncortex = sort([noncortex; 20207; 20256; 20255; 20206]);
        end
        %find(aparc == 0); [50 17 49 41]

        % To remove mesh problem
        if hemi == 'R'
            noncortex = sort([noncortex; 27806; 27810; 128; 129; 27779; 27792; 27791; 27798; ...
                27785; 27786; 27793; 27799; 27804; 27809; 27805; 27799; 27800; 27806; 27810; 27810 ]);
        end

    case '164'
        filename = ['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub '.L.midthickness.164k_fs_LR.surf.gii'];
        surf = gifti(filename);
        filename = ['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub '.L.very_inflated.164k_fs_LR.surf.gii'];
        surf_gii = gifti(filename);
        surfvi.coord = surf_gii.vertices'; surfvi.tri = surf_gii.faces;

        aparc = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub ...
            '.L.aparc.a2009s.164k_fs_LR.label.gii']);
        aparc = aparc.cdata;
        noncortex = find(aparc == 0);

    case 'freesurfer'
        surfp = SurfStatReadSurf([sub '/surf/' hemi '.pial']);
        surfw = SurfStatReadSurf([sub '/surf/' hemi '.smoothwm']);
        % find midpoint
        surf.coord = (surfp.coord + surfw.coord) ./ 2;
        surf.faces = surfp.tri;
        surf.vertices = surf.coord';
        c = read_label(sub,[hemi '.cortex']);
        cortex = c(:,1) + 1; clear c;
        noncortex = setdiff(1:length(surf.vertices),cortex);
end

index1 = ismember(surf.faces(:,1),noncortex);
index2 = ismember(surf.faces(:,2),noncortex);
index3 = ismember(surf.faces(:,3),noncortex);
index = index1 + index2 + index3;
surfN.faces = surf.faces(index == 0,:);
[incld r1 r2] = unique(surfN.faces);
surfN.vertices = surf.vertices(incld,:);
new = 1:length(incld);
surfN.faces = reshape(new(r2), size(surfN.faces));


%% Write out:
h = fopen([outputPrefix '.asc'],'w');
fprintf(h, '%5d %5d\n', [length(surfN.vertices) length(surfN.faces)]);
fprintf(h, '%4.3f\t%4.3f\t%4.3f\n',surfN.vertices');
fprintf(h, '%5i\t%5i\t%5i\n', (surfN.faces-1)');
fclose(h);

save([outputPrefix '_incld.mat'],'-v7.3','incld');


