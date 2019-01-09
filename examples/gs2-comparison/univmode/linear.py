#!/usr/bin/python

from sys import argv
import netCDF4
import numpy
from scipy import stats
import matplotlib.pyplot as plt

def get_phase( cmplx_array ):
    n = len(cmplx_array)
    phases = numpy.empty( n )
    phases[0] = 0.0
    for j in range(1,n):
        delta_phase = numpy.angle( cmplx_array[j]/cmplx_array[j-1] )
        phases[j] = phases[j-1] + delta_phase
    return phases

filename = argv[1]

GS2Output = netCDF4.Dataset(filename,'r', format="NETCDF4")

n_kx = len(GS2Output.dimensions["kx"])
n_ky = len(GS2Output.dimensions["ky"])

t_cut = 35
for i_t in range(len(GS2Output["t"])):
    if GS2Output["t"][i_t] > t_cut:
        it_cut = i_t
        break


t_vals = GS2Output["t"][it_cut:]
theta_zero_idx = ( len(GS2Output["theta"]) - 1 ) / 2
print "Averaging over the range t= {} to {}".format(t_vals[0],t_vals[-1])

for i_kx in range(n_kx):
    akx = GS2Output["kx"][i_kx]
    for i_ky in range(n_ky):
        aky = GS2Output["ky"][i_ky]
        log_phi_vals = numpy.log(GS2Output["phi2_by_mode"][it_cut:,i_ky,i_kx]) / 2.0;
        gamma, blah, r_value, blog, err = stats.linregress( t_vals, log_phi_vals );
        if r_value > 0.95:
            print "Mode at k_y = {}, k_x = {} has growth rate {} -- exponential with r={}".format(aky,akx,gamma,r_value)
        else:
            print "Mode at k_y = {}, k_x = {} is not displaying exponential time behaviour on the range t = {} to {}".format(aky,akx,t_vals[0],t_vals[-1])
            plt.plot( t_vals,log_phi_vals )
            plt.show()
            continue
        # Now do frequency
        phi_array = numpy.empty( len(t_vals), dtype=complex )
        phi_array.real = GS2Output["phi_t"][it_cut:,i_ky,i_kx,theta_zero_idx,0]
        phi_array.imag = GS2Output["phi_t"][it_cut:,i_ky,i_kx,theta_zero_idx,1]
        phase_array = get_phase( phi_array )
        
        omega, blah, r_value, blog, err = stats.linregress( t_vals, phase_array );
        print "Mode at k_y = {}, k_x = {} has frequency {} -- growth of the phase is linear with r={}".format(aky,akx,omega,r_value)


