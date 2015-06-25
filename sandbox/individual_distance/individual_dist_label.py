#!/usr/bin/python

import os, glob, subprocess, h5py
import numpy as np, scipy as sp, pandas as pd, nibabel as nib
import nipype.interfaces.freesurfer as fs
from surfer import Brain
from sklearn.utils.arpack import eigsh

# Set defaults
dataDir = '/afs/cbs.mpg.de/projects/mar005_lsd-lemon-surf/probands'
fsDir = '/afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer'
outDir = '/scr/liberia1/temp'
out_file = '%s/individual_dist_label_output_20150625.txt' % outDir
error_file = '%s/individual_dist_label_error_20150625.txt' % outDir
fwhm = 10 # for smoothing fiedler vector


# functions

def img2disc(data, foci=False, labelfile=False, hemi='lh', filename='temp.png'):
    brain = Brain('fsaverage5', hemi, 'inflated', curv=False)
    brain.add_data(data, data.min(), data.max(), colormap="spectral", alpha=0.6)
    if labelfile:
        brain.add_label(labelfile, borders=True, color='grey')
    if foci:
        brain.add_foci(foci, coords_as_verts=True, scale_factor=.7, color='black')
    brain.save_montage(filename, order=['lat', 'med'], orientation='h', border_size=10)

def runFiedler(conn):
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

def runSmoothing(data, hemi, fwhm):
    img = np.expand_dims(data, axis=0)
    img = np.expand_dims(img, axis=0)
    img = nib.freesurfer.mghformat.MGHImage(img.astype(float32), affine=None)
    img.to_filename('./temp1.mgz')
    smoothing = fs.SurfaceSmooth(subjects_dir=fsDir,
                                 subject_id='fsaverage5',
                                 in_file='./temp1.mgz',
                                 out_file='./temp2.mgz',
                                 hemi=hemi,
                                 fwhm=fwhm,
                                 cortex=True,
                                 terminal_output='none')
    smoothing.run()
    out = nib.load('./temp2.mgz').get_data().squeeze()
    os.remove('./temp2.mgz')
    os.remove('./temp1.mgz')
    return out

def runExtrema(data, hemi):
    thmin = (abs(data).max() - 1.5*abs(data).std())
    cluster = np.array([x if x > thmin else 0 for x in abs(data)])
    cluster_img = np.expand_dims(cluster, axis=0)
    cluster_img = np.expand_dims(cluster_img, axis=0)
    cluster_img = nib.freesurfer.mghformat.MGHImage(cluster_img.astype(float32), affine=None)
    cluster_img.to_filename('./temp.mgz')
    cml = 'mri_surfcluster --in ./temp.mgz --subject fsaverage5 --hemi %s --thmin %s --annot aparc.a2009s --sum ./temp.log' % (hemi, thmin)
    subprocess.call(cml, shell=True)
    extrema_log = pd.read_csv('./temp.log', skiprows=34, skipinitialspace=21, header=None, dtype={0:np.str})
    extrema_vertices = [int(extrema_log[0].iloc[i][15:24]) for i in range(len(extrema_log))]
    os.remove('./temp.log')
    os.remove('./temp.mgz')
    return extrema_vertices



output_dict = {'subject': [], 'dist_A1_parietal': [], 'dist_V1_parietal': [], 'hemisphere': []}
error_dict = {'subject': [], 'hemisphere': []}

for subject in [sub for sub in os.listdir(dataDir) if os.path.isdir(os.path.join(dataDir, sub))]:
    for hemi in ['lh', 'rh']:

        # file names and location
        corr_file = '%s/%s/correlation_maps/%s_lsd_corr_1ab_fsa5_%s.npy' % (dataDir, subject, subject, hemi)
        dist_file = '%s/%s/distance_maps/%s_%s_geoDist_fsa5.mat' % (dataDir, subject, subject, hemi)
        parietal_label_file = '%s/%s/labels/fsa5/%s.G_pariet_inf-Angular_fsa5.label' % (dataDir, subject, hemi)
        V1_label_file = '%s/%s/labels/fsa5/%s.S_calcarine_fsa5.label' % (dataDir, subject, hemi)
        A1_label_file = '%s/%s/labels/fsa5/%s.G_temp_sup-G_T_transv_fsa5.label' % (dataDir, subject, hemi)
        img1_name = '%s/%s_extrema_parietal_label_%s.png' % (outDir, subject, hemi)
        img2_name = '%s/%s_parietal_peak_%s.png' % (outDir, subject, hemi)
        

        if not False in [os.path.isfile(i) for i in [corr_file, dist_file, parietal_label_file, V1_label_file, A1_label_file]]: 
            
            print '\n\n' + subject          
            
            try:
                # read in data
                cort = np.sort(nib.freesurfer.io.read_label('%s/fsaverage5/label/%s.cortex.label' % (fsDir, hemi)))
                corr = np.load(corr_file)
                with h5py.File(dist_file, 'r') as f:
                        dist = f['dataAll'][()]
                parietal_vertices = np.sort(nib.freesurfer.io.read_label(parietal_label_file))
                V1_vertices = np.sort(nib.freesurfer.io.read_label(V1_label_file))
                A1_vertices = np.sort(nib.freesurfer.io.read_label(A1_label_file))
                # fiedler vector
                fiedler = np.zeros(len(corr))
                fiedler[cort] = runFiedler(corr[cort, :][:, cort])[:,0]
                del corr
                # smoothing, finding local extrema
                f_smoothed = runSmoothing(fiedler, hemi, fwhm)
                f_extrema_vertices = runExtrema(f_smoothed, hemi)
                img2disc(f_smoothed, foci=f_extrema_vertices, labelfile=parietal_label_file, filename=img1_name, hemi=hemi)
                # get local peak that is closest to parietal label
                dist_extrema_2_parietal = [np.mean(dist[parietal_vertices, i]) for i in f_extrema_vertices]
                parietal_peak_vertex = f_extrema_vertices[dist_extrema_2_parietal.index(min(dist_extrema_2_parietal))]
                img2disc(f_smoothed, foci=parietal_peak_vertex, filename=img2_name, hemi=hemi)
                # save results
                output_dict['subject'].append(subject)
                output_dict['hemisphere'].append(hemi)
                output_dict['dist_V1_parietal'].append(dist[V1_vertices, parietal_peak_vertex].mean())
                output_dict['dist_A1_parietal'].append(dist[A1_vertices, parietal_peak_vertex].mean())
            
            except:
                error_dict['subject'].append(subject)
                error_dict['hemisphere'].append(hemi)
                pass
            
pd.DataFrame(output_dict).to_csv(out_file, sep='\t', index=False, columns=['subject', 'hemisphere', 'dist_V1_parietal', 'dist_A1_parietal'])
pd.DataFrame(error_dict).to_csv(error_file, sep='\t', index=False, columns=['subject', 'hemisphere'])
