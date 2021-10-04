from datetime import datetime
import sdds
import numpy as np
import scipy.io as sio
import scipy.optimize as op
import os
import shutil
import gzip
import time
import glob
import pickle
from zip_mat import save_zip
import timestamp_helpers as th


class MountainRange(object):
    def __init__(self, complete_path):
        #print(f'start of mountain range class {complete_path}')
        if complete_path.endswith('.mat.gz'):
            temp_filename = complete_path.split('.gz')[0]
            with open(temp_filename, "wb") as tmp:
                shutil.copyfileobj(gzip.open(complete_path), tmp)
            dict_mr = sio.loadmat(temp_filename)
            os.remove(temp_filename)
        elif complete_path.endswith('.mat'):
            dict_mr = sio.loadmat(complete_path)
        else:
            print('Unknown file extension for MountainRange file. Should be ' +
                  '.mat or .mat.gz')
        self.value = dict_mr['value']
        self.trigger_stamp = dict_mr['triggerStamp']
        self.SC_numb = np.int(np.squeeze(dict_mr['SCNumber']))
        self.first_trigger_t_stamp_unix = dict_mr['first_trigger_t_stamp_unix']
        self.sample_interval = float(np.squeeze(dict_mr['sampleInterval']))
        self.first_sample_time = dict_mr['firstSampleTime']
        self.sensitivity = dict_mr['sensitivity']
        self.offset = dict_mr['offset']
        self.SPSuser = dict_mr['SPSuser']
        self.t_stamp_unix = dict_mr['t_stamp_unix']

        self.time_axis = np.float_(range(self.value.shape[1]))*self.sample_interval-self.value.shape[1]*self.sample_interval/2.

    def get_average_profile(self):
        average_profile = np.mean(self.value[-20:,:],axis=0)
        baseline = np.mean(average_profile[0:20])
        average_profile -= baseline
        return average_profile

    def fit_sigma_t_first_bunch(self):
        average_profile = self.get_average_profile()
        mask_average_profile = self.time_axis < self.time_axis[np.argmax(average_profile)] + 1.2
        mask_average_profile2 = self.time_axis > self.time_axis[np.argmax(average_profile)] - 1.2
        mask_average_profile = mask_average_profile & mask_average_profile2

        self._fitfunc = lambda p, x: p[0]*np.exp(-((x-p[1])**2)/(2.*p[2]**2))
        self._errfunc = lambda p, x, y: self._fitfunc(p, x) - y # Distance to the target function
        sigma0 = 1.
        ampl = np.max(average_profile[mask_average_profile])
        mu0 = self.time_axis[np.argmax(average_profile)]
        p0 = [ampl, mu0, sigma0]
        p1, success = op.leastsq(self._errfunc, p0[:], args=(self.time_axis[mask_average_profile], 
                average_profile[mask_average_profile]))
        if success:
            sigma_t = p1[2]
        else:
            sigma_t = 0.
            ampl = 0.
            mu0 = 0.
        return sigma_t, ampl, mu0 # return 1*sigma_t in ns

def sdds_to_dict_new(in_complete_path): # 4/10/2021 (Natalia)
    
    try:
        sdds_data = sdds.read(in_complete_path)
    except IndexError:
        print('Failed to open data file. (save_mr_mat)')
        return
    
    timeStamps =s dds_data.values['triggerStamp']
    time_string_filename = str(datetime.fromtimestamp(timeStamps[-1]*1e-9))[:19]
    t_stamp_unix = time.mktime(time.strptime(time_string_filename, '%Y-%m-%d %H:%M:%S'))

    dict_meas = {
        'sampleInterval': sdds_data.values['sampleInterval'],
        'firstSampleTime': sdds_data.values['firstSampleTime'],
        'cycleTag': sdds_data.values['cycleTag'],
        'offset': sdds_data.values['offset'],
        'CTime': sdds_data.values['CTime'],
        'diffStamps': sdds_data.values['diffStamps'],
        'SCNumber': sdds_data.values['SCNumber'],
        'sensitivity': sdds_data.values['sensitivity'],
        'triggerDistLatency': sdds_data.values['triggerDistLatency'],
        'value_units': sdds_data.values['value_units'],
        'triggerError':  sdds_data.values['triggerError'],
        'triggerStamp':  sdds_data.values['triggerStamp'],
        'value':  sdds_data.values['value'],
        'first_trigger_t_stamp_unix': float(sdds_data.values['triggerStamp'][0])/1e9,
        'SPSuser': in_complete_path.split('SPS.USER.')[-1].split('.')[0],
        't_stamp_unix': t_stamp_unix}

    return dict_meas



def sdds_to_file(in_complete_path, mat_filename_prefix='SPSmeas_', outp_folder='mr/'):

    dict_meas = sdds_to_dict_new(in_complete_path)
    us_string = dict_meas['SPSuser']
    t_stamp_unix = dict_meas['t_stamp_unix']

    out_filename = mat_filename_prefix + us_string + ('_%d'%t_stamp_unix)
    out_complete_path = outp_folder + us_string +'/'+ out_filename

    if not os.path.isdir(outp_folder + us_string):
        print(f'I create folder {outp_folder}{us_string}')
        #print 'I create folder: '+ outp_folder + us_string
        os.makedirs(outp_folder + us_string)

    sio.savemat(out_complete_path+'.mat', dict_meas, oned_as='row')
    save_zip(out_complete_path)



def make_pickle(pickle_name='mr_overview.pkl', mat_folder='./mr'):
    if os.path.isfile(pickle_name):
        with open(pickle_name) as fid:
            beams = pickle.load(fid)
            print(f'\n Updating file {pickle_name}')
           
    else:
        beams = {}
        #print(f'\n Creating file {pickle_name}')

    SPSuser_list = os.listdir(mat_folder)[:-1]
    
    for SPSuser in SPSuser_list:

        if not(SPSuser in beams.keys()):
            beams[SPSuser] = {}
            beams[SPSuser]['sigma_t'] = np.array([])
            beams[SPSuser]['mu0'] = np.array([])
            beams[SPSuser]['ampl'] = np.array([])
            beams[SPSuser]['timestamp_float'] = np.array([])

        list_mr_files = os.listdir(mat_folder +'/'+ SPSuser) # exclude .ipynb.checkpoints file
        files2ignore = ['.ipynb_checkpoints']
        N_cycles = len(list_mr_files)
        
        # for ii in xrange(N_cycles):
        for ii in range(N_cycles): # range() behaves like xrange() in python2
            
            filename_mr = list_mr_files[ii]
            if filename_mr not in files2ignore:
                tstamp_mat_filename = ((filename_mr.split('_')[-1]).split('.mat')[0])

                if tstamp_mat_filename in beams[SPSuser]['timestamp_float']:
                    continue
                try:
                    #print(f'{SPSuser}, {ii}, {N_cycles-1}') #'%s %d/%d'%(SPSuser, ii, N_cycles - 1)
                    curr_mr = MountainRange(mat_folder +'/'+ SPSuser +'/'+ filename_mr)
                    SPSuser = curr_mr.SPSuser[0]
                    sigmat, ampl, mu0 = curr_mr.fit_sigma_t_first_bunch()
                    beams[SPSuser]['ampl'] = np.append(beams[SPSuser]['ampl'], ampl)
                    beams[SPSuser]['mu0'] = np.append(beams[SPSuser]['mu0'], mu0)
                    beams[SPSuser]['sigma_t'] = np.append(beams[SPSuser]['sigma_t'], sigmat)
                    beams[SPSuser]['timestamp_float'] = np.append(beams[SPSuser]['timestamp_float'], curr_mr.t_stamp_unix)

                except IOError as err:
                    print(err)

        #print(beams)
        ind_sorted = np.argsort(beams[SPSuser]['timestamp_float'])
        for kk in beams[SPSuser].keys():
            beams[SPSuser][kk] = np.take(beams[SPSuser][kk], ind_sorted)
       
        with open(pickle_name, 'wb') as fid:
            pickle.dump(beams, fid)
