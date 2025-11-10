#!/usr/bin/env python
# coding: utf-8

# In[1]:


#START OF MAIN PROGRAM
#collect all file names in array
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime
import matplotlib.dates as date
import glob
import os
os.scandir(r"C:\Users\Joe\Documents\Uni\Year 3\PX3350 Physics Project\Test Data\TOI-2046b (Solved)\exoplanet transit-20240330T124750Z-001\exoplanet transit\light frames")
Files = glob.glob('*.fits')


# In[19]:


#world coords of reference stars
ref1 = [16.1389341,74.3619563]
ref2 = [16.1780348,74.3538810]
ref3 = [16.1795937,74.3109115]
ref4 = [16.1598312,74.2938497]
ref5 = [16.0957021,74.2858666]

TOI_2046b = [16.1851793,74.3313322]

#2-d reference array containing all star positions
ref = [[16.1389341,74.3619563],
       [16.1780348,74.3538810],
       [16.1795937,74.3109115],
       [16.1598312,74.2938497],
       [16.0957021,74.2858666],
       [16.1851793,74.3313322]]


# In[20]:


#define aperture function

def aperture(ra,dec,image,ap_rad,ap_bkg):
    #open image
    file = fits.open(image)
    #pixel values from image
    img_data = file[0].data
    #convert ra and dec to pixel positions for given image
    w = WCS(file[0].header)
    x_pos,y_pos = w.world_to_pixel_values(ra,dec)
    #to nearest pixel(whole number)
    x=int(x_pos)
    y=int(y_pos)
    
    #define sum and count for aperture and background
    ap_sum = 0
    N_ap = 0
    bkg_sum = 0
    N_bkg = 0  
    
    #loop over background radius (includes aperture and background)
    for i in range(-ap_bkg,ap_bkg+1):
        for j in range(-ap_bkg,ap_bkg+1):
            #if pixel position is inside aperture but not inside background, sum the pixel values and pixel number
            if (i**2+j**2<=ap_rad**2):
                ap_sum+=img_data[y+i,x+j]
                N_ap+=1
            #if pixel position is inside backgoudn but not inside aperture, sum and count pixels
            if ((i**2+j**2<=ap_bkg**2)&(i**2+j**2>ap_rad**2)):
                bkg_sum+=img_data[y+i,x+j]
                N_bkg+=1
             
    #average value for background pixel
    avg_bkg = bkg_sum/N_bkg
    #flux = aperture sum minus background for each pixel in aperture
    flux = ap_sum-(avg_bkg*N_ap)
    file.close()
    return flux
        
#apsum: sum of pixel magnitudes in central aperture 
#N_ap: number of pixels in central aperture
#N_back: number of pixels in background 
#back_total: sum of pixels magnitudes in background 


# In[21]:


import matplotlib.pyplot as plt
import numpy as np

#time array point for each image
time = np.arange(0,len(Files),1)

#create array(6,63) of zeros to fill
#6 columns for each star
#63 rows for each image
flux = np.zeros((6,len(Files)))

#for i over each star
#for j over all image files
for i in range(0,len(ref)):
    for j in range(0,len(Files)):
        #fill flux array using aperture function
        flux[i][j] = aperture(ref[i][0],ref[i][1],Files[j],15,30)    


# In[24]:


print(int(len(Files)/3+0.5))


# In[140]:


from scipy.optimize import curve_fit
from datetime import datetime
import matplotlib.dates as date
from astropy.time import Time

count = 0 
new_flux = np.zeros((len(ref),int(len(Files)/3+0.5)))
for i in range(0,len(ref)):
    count = 0
    for j in range(0,len(Files)):
        if(j%3==0):
        #fill flux array using aperture function
            new_flux[i,count] = np.mean(flux[i,j:j+3])
            count +=1


ref_frame = 59 #select the reference frame

multi_factor = np.zeros((len(ref)-1,int(len(Files)/3+0.5)))

#loop over all reference stars 
for i in range(0,len(ref)-1):
    multi_factor[i] = (new_flux[i]/new_flux[i,ref_frame])#fill multi factor array

#zeros array to fill (length 63)
multi_factor_avg = np.zeros(int(len(Files)/3+0.5))


for j in range(0,int(len(Files)/3+0.5)):
    for i in range(0,len(ref)-1):
        multi_factor_avg[j] += (multi_factor[i,j]/(len(ref)-1))


counter = 0

time = np.zeros(int(len(Files)/3+0.5)).astype(datetime)
for i in range(0,len(Files)):
    if (i%3==0):
        file = fits.open(Files[i])
        hdr = file[0].header
        time[counter] = datetime.fromisoformat(hdr['DATE-AVG'])  
        counter+=1
        
    
#date = date.date2num(time)


#zeros array to fill with normalised flux (same size as flux array)       
normalised = np.zeros((len(ref),int(len(Files)/3+0.5)))

#apply this average to the target star
normalised[len(ref)-1] = (new_flux[len(ref)-1]/new_flux[len(ref)-1,ref_frame])/multi_factor_avg
normalised_flux = normalised[len(ref)-1]

#plot target star flux 
plt.figure(figsize=(9, 5))
plt.title("Transit Light Curve")
plt.xlabel("Time [dd hh:mm]")
plt.ylabel("Relative Flux")
plt.ylim(0.97,1.01)
plt.plot(time,normalised_flux,'.', markersize=10)
plt.grid()


# In[142]:


std = np.std(multi_factor,axis=0)
std_err = std/np.sqrt(len(ref)-1)

#plot target star flux 
plt.figure(figsize=(9, 5))
plt.title("Transit Light Curve")
plt.xlabel("Time [dd hh:mm]")
plt.ylabel("Relative Flux")
plt.ylim(0.97,1.01)

plt.errorbar(time,normalised_flux,yerr=std_err,fmt ='.',markersize=10)
#plt.plot(date,normalised_flux,'.', markersize=10)
plt.grid()


# In[143]:


from astropy.timeseries import BoxLeastSquares
import batman
from scipy.optimize import minimize
#import corner
import numpy as np
import matplotlib.pyplot as plt
#import emcee
import astropy.units as u


# In[144]:


#max_power = np.argmax(periodogram.power)
period_bls = 1.5
epoch_bls = 0
depth_bls = 0.0150
#duration_bls = periodogram.duration[max_power]
t = np.linspace(-0.03, 0.03, 66)  #times at which to calculate light curve

depth_value = np.zeros(350)
scaled_seperation_value = np.zeros(350)
inclination_value = np.zeros(350)


# In[145]:


def log_likelihood(parameters, times, observed_flux):
    global count

    depth_value[count] = parameters[2]
    scaled_seperation_value[count] = parameters[3]
    inclination_value[count] = parameters[4]
    # Generate the model flux using the batman package and the provided parameters.
    model_flux = f_batman(times, *parameters)
    
    # Calculate the variance of the observed flux (squared error).
   
    # Compute the log-likelihood using the Gaussian log-likelihood formula.
    log_likelihood_value = -0.5 * np.sum(
        (observed_flux - model_flux) ** 2 
    )
    
    count +=1
    return log_likelihood_value


def f_batman(
    times,
    time_of_conjunction,
    orbital_period,
    planet_radius,
    semi_major_axis,
    orbital_inclination,
    baseline_flux=0.0,
    eccentricity=0,
    longitude_of_periastron=90,
    limb_darkening_coefficients=[0.3, 0.28],
    limb_darkening_model="quadratic",
):
   
    # Initialize the parameters for the batman transit model
    params = batman.TransitParams()
    params.t0 = time_of_conjunction
    params.per = orbital_period
    params.rp = planet_radius
    params.a = semi_major_axis
    params.inc = orbital_inclination
    params.ecc = eccentricity
    params.w = longitude_of_periastron
    params.u = limb_darkening_coefficients
    params.limb_dark = limb_darkening_model

    # Initialize the batman transit model and calculate the light curve
    transit_model = batman.TransitModel(params, times)
    flux_model = transit_model.light_curve(params)

    # Add the baseline flux to the model output and return
    return np.array(flux_model) + baseline_flux


# In[146]:


# Constants for conversion
SOLAR_RADIUS_IN_AU = 215.0  # Solar radius in astronomical units

# Calculation of the semi-major axis in solar radii assuming a 1 Msol, 1 Rsol host

semi_major_axis = 5
# Initial guess for the transit model parameters
initial_guess = np.array(
    [
        epoch_bls,  # Time of central transit
        period_bls,  # Orbital period
        depth_bls**0.5,  # Approximation of planet radius
        semi_major_axis,  # Semi-major axis
        86,  # Orbital inclination (degrees)
        0.0,  # Baseline flux offset
    ]
)

count = 0
# Define the negative log likelihood function for optimization
negative_log_likelihood = lambda *args: -log_likelihood(*args)

# Perform optimization to find the best-fit parameters
solution = minimize(
    negative_log_likelihood,
    initial_guess,
    args=(t, normalised_flux), 
)
print(count)
# Parameter labels for output
parameter_labels = [
    "t0 (Epoch)",
    "per (Period)",
    "rp (Planet Radius)",
    "a (Semi-major Axis)",
    "inc (Inclination)",
    "baseline (Flux Offset)",
]


# Print the header of the table
header = f"{'Parameter':<25} {'Initial Value':<15} {'Fitted Value':<15}"
print(header)
print("-" * len(header))
#solution.x[1]=1.75
# Iterate over each parameter and print the values in a table-like format
for idx in range(len(solution.x)):
    parameter_label = parameter_labels[idx]
    initial_value = initial_guess[idx]
    fitted_value = solution.x[idx]

    row = f"{parameter_label:<25} {initial_value:<15.4f} {fitted_value:<15.4f}"
    print(row)


# In[147]:


std_dev_depth = np.std(depth_value)
std_dev_scaled_seperation = np.std(scaled_seperation_value)
std_dev_inclination = np.std(inclination_value)


print(std_dev_depth)
print(std_dev_scaled_seperation)
print(std_dev_inclination)


# In[148]:


# Constants
NUMBER_OF_POINTS = 10000
PHASE_FOLD_THRESHOLD = 0.5

# Create a time array for the transit model
model_times = np.linspace(np.min(t), np.max(t), NUMBER_OF_POINTS)

# Compute the initial model flux and the fitted model flux
initial_model_flux = f_batman(model_times, *initial_guess)
fitted_model_flux = f_batman(model_times, *solution.x)

# Phase-fold the observed data time-array
phase_observed = (t - solution.x[0]) % solution.x[1] / solution.x[1]
phase_observed[phase_observed > PHASE_FOLD_THRESHOLD] -= 1

# Phase-fold the BLS model time-array and sort
phase_bls_model = (model_times - epoch_bls) % period_bls / period_bls
phase_bls_model[phase_bls_model > PHASE_FOLD_THRESHOLD] -= 1
sorted_bls_model_idx = np.argsort(phase_bls_model)

# Phase-fold the fitted model time-array and sort
phase_fitted_model = (model_times - solution.x[0]) % solution.x[1] / solution.x[1]
phase_fitted_model[phase_fitted_model > PHASE_FOLD_THRESHOLD] -= 1
sorted_fitted_model_idx = np.argsort(phase_fitted_model)

# Plotting
# Time series plot of data and models
plt.figure(figsize=(10, 5))
plt.plot(t, normalised_flux, "o", color="xkcd:dusty red", label="Data (Time Series)", alpha=0.4)
#plt.plot(model_times, initial_model_flux, "k", label="Initial BLS Model", alpha=0.6)
plt.plot(model_times, fitted_model_flux, "--k", label="Fitted Model")
plt.xlabel("Time [days]")
plt.ylabel("Normalized Flux")
plt.legend()
plt.grid

# Phase-folded plot of data and fitted model
plt.figure(figsize=(12, 5))
#plt.errorbar(date,normalised_flux,yerr=std_err,fmt ='.',markersize=10)
plt.errorbar(
    phase_observed * solution.x[1],
    normalised_flux,
    yerr=std_err,
    fmt = "o",
    color="xkcd:dusty red",
    label="Data",
    markersize=5,
    alpha=0.9,
)
plt.plot(
    phase_fitted_model[sorted_fitted_model_idx] * solution.x[1],
    fitted_model_flux[sorted_fitted_model_idx],
    "--k",
    label="Fitted Model",
)
plt.xlabel("Time from Central Transit [days]")
plt.ylabel("Normalized Flux")
plt.xlim(-0.03, 0.035)
plt.legend()
plt.grid()


# In[136]:


start_time = time[normalised_flux==normalised_flux[7]]
end_time = time[normalised_flux==normalised_flux[59]]

t = np.arange(66)
y = np.linspace(1.025,0.995,66)


start_point = np.full(66,start_time)
end_point = np.full(66,end_time)


difference = time[normalised_flux==normalised_flux[7]]-time[normalised_flux==normalised_flux[0]]
print(difference)


# In[137]:


plt.figure(figsize=(12, 5))
plt.plot(start_point, y)
plt.plot(end_point,y)
plt.plot(time,normalised_flux,'.')


# In[ ]:




