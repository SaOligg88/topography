#!/usr/bin/python

from multiprocessing import Pool
import numpy as np, pandas as pd
import os


# Set defaults
dataDir = '/afs/cbs.mpg.de/projects/mar005_lsd-lemon-surf/probands'
fsDir = '/afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer'
out_file = '/scr/liberia1/data/individual_dist_label/res_individual_dist_label_lh_20150714.txt'


### !!! hemi needs to be changed within function !!! ###

# ----------------------------------- functions -----------------------------------

def run_individual_dist_label(subject):    
    
    hemi = 'lh'
    import os, glob, subprocess, h5py
    import numpy as np, pandas as pd, nibabel as nib
    import nipype.interfaces.freesurfer as fs
    from surfer import Brain
    from sklearn.utils.arpack import eigsh

    dataDir = '/afs/cbs.mpg.de/projects/mar005_lsd-lemon-surf/probands'
    fsDir = '/afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer'
    outDir = '/scr/liberia1/data/individual_dist_label'
    
    def img2disc(data, foci=False, labelfile=False, hemi='lh', filename='temp.png'):
        brain = Brain('fsaverage5', hemi, 'inflated', curv=False)
        brain.add_data(data, data.min(), data.max(), colormap="spectral", alpha=0.6)
        if labelfile:
            brain.add_label(labelfile, borders=True, color='grey')
        if foci:
            brain.add_foci(foci, coords_as_verts=True, scale_factor=.7, color='black')
        brain.save_montage(filename, order=['lat', 'med'], orientation='h', border_size=10)
    
    def runFiedler(conn):
        # https://github.com/margulies/topography
        # prep for embedding
        K = (conn + 1) / 2.
        v = np.sqrt(np.sum(K, axis=1))
        A = K/(v[:, None] * v[None, :])
        del K
        A = np.squeeze(A * [A > 0])
        # diffusion embedding
        n_components_embedding = 2
        lambdas, vectors = eigsh(A, k=n_components_embedding+1)
        del A
        lambdas = lambdas[::-1]
        vectors = vectors[:, ::-1]
        psi = vectors/vectors[:, 0][:, None]
        lambdas = lambdas[1:] / (1 - lambdas[1:])
        embedding = psi[:, 1:(n_components_embedding + 1 + 1)] * lambdas[:n_components_embedding+1][None, :]
        return embedding
    
    def runMasking(data, hemi):
        mask = np.zeros((10242))
        for label in [39, 40, 46, 47, 49, 50, 51, 68, 85, 86]:
            label = np.sort(nib.freesurfer.io.read_label(glob.glob('%s/fsaverage5/label/*%s*label*' % (fsDir, hemi))[label]))
            mask[label] = 1
        masked = data * mask
        return masked
    
    def runSmoothing(data, hemi, subject):
        temp1 = './temp1_%s.mgz' % subject
        temp2 = './temp2_%s.mgz' % subject
        img = np.expand_dims(data, axis=0)
        img = np.expand_dims(img, axis=0)
        img = nib.freesurfer.mghformat.MGHImage(img.astype(float32), affine=None)
        img.to_filename(temp1)
        smoothing = fs.SurfaceSmooth(subjects_dir=fsDir,
                                     subject_id='fsaverage5',
                                     in_file=temp1,
                                     out_file=temp2,
                                     hemi=hemi,
                                     fwhm=5,
                                     cortex=True,
                                     terminal_output='none')
        smoothing.run()
        out = nib.load(temp2).get_data().squeeze()
        os.remove(temp1)
        os.remove(temp2)
        return out
    
    def runExtrema(data, hemi, subject):
        temp1 = './temp_%s.mgz' % subject
        temp2 = './temp_%s.log' % subject        
        thmin = (abs(data).max() - 1.3*abs(data).std())
        cluster = np.array([x if x > thmin else 0 for x in abs(data)])
        cluster_img = np.expand_dims(cluster, axis=0)
        cluster_img = np.expand_dims(cluster_img, axis=0)
        cluster_img = nib.freesurfer.mghformat.MGHImage(cluster_img.astype(float32), affine=None)
        cluster_img.to_filename(temp1)
        cml = 'mri_surfcluster --in %s --subject fsaverage5 --hemi %s --thmin %s --annot aparc.a2009s --sum %s' % (temp1, hemi, thmin, temp2)
        subprocess.call(cml, shell=True)
        extrema_log = pd.read_csv(temp2, skiprows=34, skipinitialspace=21, header=None, dtype={0:np.str})
        extrema_vertices = [int(extrema_log[0].iloc[i][15:25]) for i in range(len(extrema_log))]
        os.remove(temp1)
        os.remove(temp2)
        return extrema_vertices
    
    
        
    # file names and location
    corr_file1 = '%s/%s/correlation_maps/%s_lsd_corr_1ab_fsa5_%s.npy' % (dataDir, subject, subject, hemi)
    corr_file2 = '%s/%s/correlation_maps/%s_lsd_corr_2ab_fsa5_%s.npy' % (dataDir, subject, subject, hemi)
    dist_file = '%s/%s/distance_maps/%s_%s_geoDist_fsa5.mat' % (dataDir, subject, subject, hemi)
    parietal_label_file = '%s/%s/labels/fsa5/%s.G_pariet_inf-Angular_fsa5.label' % (dataDir, subject, hemi)
    temporal_label_file = '%s/%s/labels/fsa5/%s.Pole_temporal_fsa5.label' % (dataDir, subject, hemi)
    V1_label_file = '%s/%s/labels/fsa5/%s.S_calcarine_fsa5.label' % (dataDir, subject, hemi)
    A1_label_file = '%s/%s/labels/fsa5/%s.G_temp_sup-G_T_transv_fsa5.label' % (dataDir, subject, hemi)
    fiedler_file = '%s/fiedler/%s_fiedler_%s' % (outDir, subject, hemi)
    peak_img_file = '%s/qc/%s_fiedler_dmnExtrema_%s.png' % (outDir, subject, hemi)
            
    
    
    try:
        #if not False in [os.path.isfile(i) for i in [corr_file, dist_file, parietal_label_file, temporal_label_file, V1_label_file, A1_label_file]]: 
        # read in data
        cort = np.sort(nib.freesurfer.io.read_label('%s/fsaverage5/label/%s.cortex.label' % (fsDir, hemi)))
        corr1 = np.load(corr_file1)
        corr2 = np.load(corr_file2)
        corr = (corr1+corr2) /2
        with h5py.File(dist_file, 'r') as f:
            dist = f['dataAll'][()]
        parietal_vertices = np.sort(nib.freesurfer.io.read_label(parietal_label_file))
        temppole_vertices = np.sort(nib.freesurfer.io.read_label(temporal_label_file))
        V1_vertices = np.sort(nib.freesurfer.io.read_label(V1_label_file))
        A1_vertices = np.sort(nib.freesurfer.io.read_label(A1_label_file))
                
        # local extrema in fiedler vector
        fiedler = np.zeros(len(corr))
        fiedler[cort] = runFiedler(corr[cort, :][:, cort])[:,0]
        del corr
        f_smoothed = runSmoothing(fiedler, hemi, subject)
        f_masked = runMasking(f_smoothed, hemi)
        f_extrema_vertices = runExtrema(f_masked, hemi, subject)
                
        # distances
        dist_extrema_2_parietal = [np.mean(dist[parietal_vertices, i]) for i in f_extrema_vertices]
        parietal_peak_vertex = f_extrema_vertices[dist_extrema_2_parietal.index(min(dist_extrema_2_parietal))]
        dist_extrema_2_temporal = [np.mean(dist[temppole_vertices, i]) for i in f_extrema_vertices]
        temporal_peak_vertex = f_extrema_vertices[dist_extrema_2_temporal.index(min(dist_extrema_2_temporal))]
        
        # save standardized fiedler
        if fiedler[parietal_peak_vertex] < 0:
            f_stand = -fiedler
        else:
            f_stand = fiedler
        
        np.save(fiedler_file, f_stand)
        img2disc(f_stand, foci=[parietal_peak_vertex, temporal_peak_vertex], hemi=hemi, filename=peak_img_file)

        # return results
        V1_vertices = nib.freesurfer.io.read_label('%s/%s/labels/fsa5/%s.S_calcarine_fsa5.label' % (dataDir, subject, hemi))
        A1_vertices = nib.freesurfer.io.read_label('%s/%s/labels/fsa5/%s.G_temp_sup-G_T_transv_fsa5.label' % (dataDir, subject, hemi))
        V1_parietal = dist[V1_vertices, parietal_peak_vertex].mean()
        A1_parietal = dist[A1_vertices, parietal_peak_vertex].mean()
        V1_temporal = dist[V1_vertices, temporal_peak_vertex].mean()
        A1_temporal = dist[A1_vertices, temporal_peak_vertex].mean()
        return subject, hemi, V1_parietal, A1_parietal, V1_temporal, A1_temporal

    except:
        return subject, hemi, None, None, None, None
        pass



# --------------------------------------------------------------------------------------------------------------------------------------------------

### run serially ###

subjects = [sub for sub in os.listdir(dataDir) if os.path.isdir(os.path.join(dataDir, sub))]
output_dict = {}
res = []

for subject in subjects:
    res.append(run_individual_dist_label(subject))

output_dict['subject'], output_dict['hemi'], output_dict['V1_parietal'], output_dict['A1_parietal'], output_dict['V1_temporal'], output_dict['A1_temporal'] = np.array(res).T    
pd.DataFrame(output_dict).to_csv(out_file, sep='\t', index=False, columns=['subject', 'hemi', 'V1_parietal', 'A1_parietal', 'V1_temporal', 'A1_temporal'])
    

### run in parallel ### (not good for qc screenshots though)

#p = Pool(20)
#res = p.map(run_individual_dist_label, subjects)    
#output_dict = {}
#output_dict['subject'], output_dict['hemi'], output_dict['V1_parietal'], output_dict['A1_parietal'], output_dict['V1_temporal'], output_dict['A1_temporal'] = np.array(res).T    
#pd.DataFrame(output_dict).to_csv(out_file, sep='\t', index=False, columns=['subject', 'hemi', 'V1_parietal', 'A1_parietal', 'V1_temporal', 'A1_temporal'])