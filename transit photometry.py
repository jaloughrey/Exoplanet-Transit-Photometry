
#START OF MAIN PROGRAM

from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime
import matplotlib.dates as date
from astropy.time import Time
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares
import batman
from scipy.optimize import minimize
import astropy.units as u
from scipy.optimize import curve_fit




folder = r"C:\Users\Joe\Documents\Uni\Year 3\PX3350 Physics Project\Test Data\TOI-2046b (Solved)\exoplanet transit-20240330T124750Z-001\exoplanet transit\light frames\LIGHTS"

full_paths = glob.glob(folder + r"\*.fits")
Files = [os.path.basename(f) for f in full_paths]

#read coordinates file
coords = pd.read_csv('coordinates.csv')
refs = coords[coords['type']=='ref'][['ra','dec']].values
target = coords[coords['type']=='target'][['ra','dec']].values[0]
ref = np.vstack([refs, target])

#read paramter file
param = pd.read_csv('parameters.csv')
target_star = param['target_star'].values[0]
ap_rad = param['ap_rad'].values[0]
ap_bkg = param['ap_bkg'].values[0]
reference_frame = param['reference_frame'].values[0]
group_size = param['group_size'].values[0]



#define aperture function
def aperture(ra,dec,image,ap_rad,ap_bkg):
    #open image
    file = fits.open(os.path.join(folder, image))
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



#time array point for each image
time = np.arange(0,len(Files),1)

#create array(6,63) of zeros to fill
#6 columns for each star
#63 rows for each image
flux = np.zeros((len(ref),len(Files)))

#for i over each star
#for j over all image files
for i in range(0,len(ref)):
    for j in range(0,len(Files)):
        #fill flux array using aperture function
        flux[i][j] = aperture(ref[i][0],ref[i][1],Files[j],ap_rad,ap_bkg)    





# --- Dynamic binning based on target_points ---


num_frames = len(Files)
num_bins = int(np.ceil(num_frames / group_size))


print("Output bins =", num_bins)

# --- Create binned flux array (replaces new_flux) ---

new_flux = np.zeros((len(ref), num_bins))

for i in range(len(ref)):  # each star
    count = 0
    for j in range(0, num_frames, group_size):
        new_flux[i, count] = np.mean(flux[i, j:j + group_size])
        count += 1


# --- Normalisation section (mostly unchanged) ---

ref_frame = reference_frame   # from your parameters.csv

multi_factor = np.zeros((len(ref)-1, num_bins))

# loop over all reference stars (not including target)
for i in range(len(ref)-1):
    multi_factor[i] = new_flux[i] / new_flux[i, ref_frame]

# average reference star correction
multi_factor_avg = np.zeros(num_bins)

for j in range(num_bins):
    multi_factor_avg[j] = np.mean(multi_factor[:, j])


# --- Time binning ---

time = np.zeros(num_bins, dtype=datetime)

count = 0
for j in range(0, num_frames, group_size):
    hdr = fits.getheader(os.path.join(folder, Files[j]))
    time[count] = datetime.fromisoformat(hdr["DATE-AVG"])
    count += 1


# --- Final normalisation of target star flux ---

normalised = np.zeros((len(ref), num_bins))

target_index = len(ref) - 1   # target is last row
normalised[target_index] = (
    (new_flux[target_index] / new_flux[target_index, ref_frame]) /
    multi_factor_avg
)

normalised_flux = normalised[target_index]



#calculate error bars for each data point
std = np.std(multi_factor,axis=0)
std_err = std/np.sqrt(len(ref)-1)

# --- Plot light curve ---

#plot target star flux 
plt.figure(figsize=(9, 5))
plt.title("Transit Light Curve")
plt.xlabel("Time [dd hh:mm]")
plt.ylabel("Relative Flux")
plt.ylim(0.97,1.01)

plt.errorbar(time,normalised_flux,yerr=std_err,fmt ='.',markersize=10)
#plt.plot(date,normalised_flux,'.', markersize=10)
plt.grid()


# In[35]:


import batman
import numpy as np
import emcee

def f_batman(times, t0, period, rp, a_rs, inc, baseline=0.0, ecc=0.0, omega=90.0, u=[0.3, 0.28], ld_model="quadratic"):
  
    params = batman.TransitParams()
    params.t0 = t0
    params.per = period
    params.rp = rp
    params.a = a_rs
    params.inc = inc
    params.ecc = ecc
    params.w = omega
    params.u = u
    params.limb_dark = ld_model
    
    # Create the model and compute the light curve
    model = batman.TransitModel(params, times)
    flux = model.light_curve(params)
    
    return np.array(flux) + baseline

def generate_model(params, period, times):
    
    t0, rp, a_rs, inc, baseline = params[:5]
    
    flux = f_batman(
        times=times,
        t0=t0,
        period=period,
        rp=rp,
        a_rs=a_rs,
        inc=inc,
        baseline=baseline
    )
    
    return flux

def log_prior(params, t_min=None, t_max=None):
    
    t0, rp, a_rs, inc, baseline, ln_sigma = params
    
    # ----- PRIOR 1: t0 (transit midpoint) -----
    if t_min is not None and t_max is not None:
        if not (t_min - 0.2 <= t0 <= t_max + 0.2):
            return -np.inf
    else:
        # fallback broad prior
        if not (-1.0 <= t0 <= 1.0):
            return -np.inf
    
    # ----- PRIOR 2: rp (radius ratio) -----
    if not (0.0 < rp < 0.5):
        return -np.inf
    
    # ----- PRIOR 3: a/R* -----
    if not (1.0 < a_rs < 100.0):
        return -np.inf
    
    # ----- PRIOR 4: inclination (degrees) -----
    if not (75.0 < inc < 90.0):
        return -np.inf
    
    # ----- PRIOR 5: baseline -----
    if not (-0.05 < baseline < 0.05):
        return -np.inf
    
    # ----- PRIOR 6: ln(sigma) -----
    if not (-20.0 < ln_sigma < -1.0):
        return -np.inf
    
    # If all priors satisfied → returning 0 means uniform priors
    return 0.0

def log_likelihood(params, period, times, flux):
    
    # unpack parameters
    t0, rp, a_rs, inc, baseline, ln_sigma = params
    sigma = np.exp(ln_sigma)  # noise must be positive
    
    # compute model flux
    model_flux = generate_model([t0, rp, a_rs, inc, baseline], 
                                period,
                                times,)
    
    # Gaussian log likelihood
    residuals = flux - model_flux
    inv_var = 1.0 / (sigma**2)
    
    likelihood = -0.5 * np.sum(residuals**2 * inv_var + np.log(2 * np.pi * sigma**2))
    
    return likelihood

def log_posterior(params, period, times, flux):
    lp = log_prior(params, t_min=np.min(times), t_max=np.max(times))
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, period, times, flux)

# wrapper to freeze t and flux
def log_prob_wrapper(params):
    return log_posterior(params, period, t, normalised_flux)

def summarize_parameters(samples, labels):
    results = {}

    for i, label in enumerate(labels):
        param = samples[:, i]
        median = np.percentile(param, 50)
        lower  = np.percentile(param, 16)
        upper  = np.percentile(param, 84)
        
        results[label] = {
            "median": median,
            "minus_1sigma": median - lower,
            "plus_1sigma": upper - median,
        }

    return results


def run_emcee(initial_guess, period, log_posterior, nwalkers=32, burn_steps=2000, nsteps = 10000, rng=None):
    #ndim: number of parameters for the MCMC
    ndim = len(initial_guess)
    print("Number of parameters:", ndim)
    print("Number of walkers:", nwalkers)
    
    if rng is None:
        rng = np.random.default_rng()
    
    #initialise walkers
    p0_center = np.array(initial_guess, dtype=float)
    scales = np.maximum(np.abs(p0_center) * 1e-4, 1e-4) #small spread
    p0 = p0_center + rng.normal(scale=scales, size=(nwalkers, ndim))
    
    print("Initial walker positions (p0) shape:", p0.shape)
    
    #create sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior)
    
    print("Sampler created with", nwalkers, "walkers and", ndim, "parameters.")

    #run burn-in
    print("Running burn-in for", burn_steps, "steps...")
    
    state = sampler.run_mcmc(p0, burn_steps, progress=True)
    
    burnin_pos = state.coords #burn in end point
    
    print("Burn-in complete. Resetting sampler to discard burn-in samples.")
    sampler.reset()  #discard burn-in

    #run main sampling 

    print(f"Running main MCMC run for {nsteps} steps...")

    final_state = sampler.run_mcmc(burnin_pos, nsteps, progress=True)

    #extract the sample chains 
    samples = sampler.get_chain(flat=True)
    
    print("Main MCMC run complete.")
    print("Total samples collected:", len(samples))

    return samples
    


# In[36]:


t = np.linspace(-0.03, 0.03, 66)

#period = database lookup
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

result = NasaExoplanetArchive.query_criteria(
    table="pscomppars",
    select="pl_name, pl_orbper",
    where=f"pl_name = '{target_star}'"
)

period = result["pl_orbper"][0].value
print(period)


"""
params = [
    t0,         # time of transit center
    rp,         # planet radius / stellar radius
    a_rs,       # semi-major axis / stellar radius
    inc,        # inclination in degrees
    baseline,   # flux offset
    ln_sigma    # log of observational noise
]
"""


#produce initial parameters 
# params order: [t0, period, rp, a/Rs, inc, baseline, ln_sigma]

initial_params = [0.0, 0.12, 5.0, 86, 0.0, np.log(0.0005)]
print("initial params:", initial_params)


samples = run_emcee(initial_params, period, log_prob_wrapper)


labels = ["t0", "rp", "a/Rs", "inc", "baseline", "ln_sigma"]
results = summarize_parameters(samples, labels)

#print
for key, val in results.items():
    m = val["median"]
    l = val["minus_1sigma"]
    u = val["plus_1sigma"]
    print(f"{key:10s} = {m:.6f} -{l:.6f} +{u:.6f}")


# In[28]:


import corner


fig = corner.corner(
    samples,
    labels=labels,
    show_titles=True,
    title_fmt=".4f",
    quantiles=[0.16, 0.5, 0.84],
    title_kwargs={"fontsize": 12},
)

plt.show()


# In[37]:


# Compute median parameter vector
median_params = np.median(samples, axis=0)

print("Median parameter vector:")
for name, val in zip(labels, median_params):
    print(f"{name:10s} = {val:.6f}")


t_model = np.linspace(np.min(t), np.max(t), 2000)

best_model = generate_model(median_params, period, t_model)

plt.figure(figsize=(10,5))

plt.scatter(t, normalised_flux, s=12, color="black", alpha=0.6, label="Data")
plt.plot(t_model, best_model, color="red", lw=2, label="Best-fit model")

plt.xlabel("Time [days]")
plt.ylabel("Normalized Flux")
plt.legend()
plt.grid()
plt.show()


# In[38]:


plt.figure(figsize=(10,5))

plt.scatter(t, normalised_flux, s=12, color="black", alpha=0.6, label="Data")

# Draw 50 random posterior samples
idx = np.random.choice(len(samples), 50, replace=False)

for i in idx:
    m = generate_model(samples[i],period, t_model)
    plt.plot(t_model, m, color="red", alpha=0.05)

plt.plot(t_model, best_model, color="red", lw=2, label="Median model")

plt.xlabel("Time [days]")
plt.ylabel("Normalized Flux")
plt.legend()
plt.grid()
plt.show()


# In[39]:


# Extract median parameters
t0_med, per_med, *_ = median_params

# Compute phase
phase_obs = ((t - t0_med + 0.5*per_med) % per_med) / per_med - 0.5

# Fine time grid spanning the entire transit window
t_model = np.linspace(np.min(t), np.max(t), 3000)

# Generate model flux for these times
model_flux = generate_model(median_params, period, t_model)

# Compute model phases
phase_model = ((t_model - t0_med + 0.5*per_med) % per_med) / per_med - 0.5

# Sort by phase to make a nice continuous curve
sort_idx = np.argsort(phase_model)
phase_model_sorted = phase_model[sort_idx]
model_flux_sorted = model_flux[sort_idx]

plt.figure(figsize=(10, 5))

# Plot observed data
plt.errorbar(phase_obs, normalised_flux,yerr=std_err,fmt = "o", color="xkcd:dusty red", alpha=0.9, label="Data",markersize=5,)

# Plot model
plt.plot(phase_model_sorted, model_flux_sorted,"--k", lw=2, label="Best-fit Model")

plt.xlabel("Orbital Phase")
plt.ylabel("Normalized Flux")
plt.title("Phase-Folded Transit Light Curve")
#plt.xlim(-0.04, 0.04)   # zoom around transit
plt.legend()
plt.grid()
plt.show()


# In[ ]:




