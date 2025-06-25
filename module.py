import numpy as np
def compute_A(x,f,n_IRE=2.4,theta=45):
    '''
    Compute the absorbance spectrum of water using refractive index data from Segelstein et al. (1981).

    Inputs:
        x (np.ndarray): new wavenumber axis in cm⁻¹ for interpolation.
        f (float): s-polarization fraction (0 ≤ f ≤ 1).
    n_IRE (float): refractive index of the internal reflection element (IRE) (default is 2.4).
    theta (float, optional): angle of incidence in degrees (default is 45°).

    Returns:
        absorbance_interp (np.ndarray): Interpolated absorbance values on the new wavenumber axis `x`.
    '''
    # Load water refractive index data
    h2o_cm, n_real, n_imag = obtain_reference_data()
    n_water = n_real + 1j * n_imag
    
    # Convert angle to radians
    theta_rad = np.radians(theta)
    
    # Define constants
    n_air = 1.0
    sin_theta_sq = np.sin(theta_rad) ** 2
    cos_theta = np.cos(theta_rad)

    # Common terms
    cos_t2_air = 1j * np.sqrt((n_IRE / n_air)**2 * sin_theta_sq - 1)
    cos_t2_h2o = 1j * np.sqrt((n_IRE / n_water)**2 * sin_theta_sq - 1)

    # Reflectance for air
    r_s_air = (n_IRE * cos_theta - n_air * cos_t2_air) / (n_IRE * cos_theta + n_air * cos_t2_air)
    r_p_air = (n_IRE * cos_t2_air - n_air * cos_theta) / (n_IRE * cos_t2_air + n_air * cos_theta)
    R_air = f * np.abs(r_s_air)**2 + (1 - f) * np.abs(r_p_air)**2

    # Reflectance for water
    r_s_h2o = (n_IRE * cos_theta - n_water * cos_t2_h2o) / (n_IRE * cos_theta + n_water * cos_t2_h2o)
    r_p_h2o = (n_IRE * cos_t2_h2o - n_water * cos_theta) / (n_IRE * cos_t2_h2o + n_water * cos_theta)
    R_h2o = f * np.abs(r_s_h2o)**2 + (1 - f) * np.abs(r_p_h2o)**2

    # Compute absorbance and interpolate
    absorbance = -np.log10(R_h2o / R_air)
    absorbance_interp = np.interp(x, h2o_cm[::-1], absorbance[::-1])

    return(absorbance_interp)
#-------------------------------

def A_to_R(A,n_IRE=2.4,f=0.5,theta=45):
    '''
    Computes reflectance from absorbance for ATR

    Inputs:
        A (np.ndarray): absorbance spectrum
        n_IRE (float, optional): refraction index of the internal reflection element (IRE) (default is 2.4)
        f (float, optional): s-polarization fraction (default is 0.5)
        theta (float, optional): angle of incidence in degrees (default is 45°)

    Returns:
        R (np.ndarray): Reflectance spectrum.
    '''
    
    theta_rad = np.radians(theta)
    n_air = 1.0
    
    sin_theta = np.sin(theta_rad)
    cos_t2_air = 1j * np.sqrt((n_IRE / n_air)**2 * sin_theta**2 - 1)

    
    # Fresnel reflection coefficients (air-IRE interface)
    r_s_air = (n_IRE * np.cos(theta_rad) - n_air * cos_t2_air) / (n_IRE * np.cos(theta_rad) + n_air * cos_t2_air)
    r_p_air = (n_IRE * cos_t2_air - n_air * np.cos(theta_rad)) / (n_IRE * cos_t2_air + n_air * np.cos(theta_rad))
    
     # Reflectance for s and p polarizations
    R_s_air = np.abs(r_s_air)**2
    R_p_air = np.abs(r_p_air)**2

    # Total reflectance for air interface
    R_air = f * R_s_air + (1 - f) * R_p_air

    # Convert absorbance to reflectance
    R = R_air * 10**(-A)

    return R


def get_pol_angle(f):
    """
    Calculate the polarization angle (in radians) from the polarization fraction.

    The polarization angle phi is given by:
        phi = arccos(sqrt(f))

    Inputs:
        f (float): The polarization fraction, a scalar between 0 and 1 (inclusive).

    Returns:
    phi (float): The polarization angle phi in radians, within the range [0, pi/2].


    """

    phi = np.arccos(np.sqrt(f))
    return phi



def R_S_45 (phi_rad, R):

    """
    Compute the s-polarized reflectance at 45° incidence (R_S_45) given the polarization angle
    and experimental absorbance.


    Inputs:
    phi_rad (float): Polarization angle in radians.

    R (np.ndarray): Experimental reflectance

    Returns:
    R_S (np.ndarray): Computed reflectance for s-polarized light at 45° incidence.

    """
    cot_phi = np.cos(phi_rad) / np.sin(phi_rad)
    csc_phi = 1 / np.sin(phi_rad)
    SR = np.sqrt(np.cos(phi_rad)**4 + 4*R*np.sin(phi_rad)**2)
    R_S = 0.5 * (-(cot_phi**2) + csc_phi**2 * SR)
    return( R_S)



def sumatory(x, x_list, R_S):
    """
    Compute the summation used in the phi_x calculation.

    Parameters:
    ----------
    x : float
        Wavenumber at which to evaluate the summation (cm⁻¹).
    x_list : array-like of float
        List of wavenumber values (cm⁻¹).
    R_S : numpy array
        Corresponding reflectance values R_s at each x in x_list.

    Returns:
    -------
    float
        The computed summation value used in phi_x.
    """
    x_arr = np.array(x_list)
    R_S_arr = np.array(R_S)

    # Avoid division by zero where x == x_list[i]
    mask = x_arr != x
    terms = x_arr[mask] * np.log(R_S_arr[mask]) / (x_arr[mask]**2 - x**2)

    return np.sum(terms)


def phi_x(x, x_list, R_S, n_inf, S, n_IRE):
    """
    Compute the phase shift phi(x) for wavenumber  x.

    Parameters:
    ----------
    x : float
        Wavenumber at which to compute phi (cm⁻¹).
    x_list : array-like of float
        List of wavenumber values (cm⁻¹).
    R_S : array-like of float
        Reflectance values R_s at each wavenumber.
    n_inf : float
        High-frequency limit of the refractive index of the sample.
    S : float
        Oscillator strength term.
    n_IRE : float
        Refractive index of the internal reflection element (IRE).

    Returns:
    -------
    float
        Phase shift phi(x) in radians.
    """
    d = (x_list[1] - x_list[0]) / 2  

    term_1 = -(2 * d / np.pi) * sumatory(x, x_list, R_S)

    theta = np.pi / 4  # 45 degrees incidence
    term_2 = np.pi - 2 * np.arctan(
        np.sqrt(n_IRE**2 * np.sin(theta)**2 - n_inf**2) / (n_IRE * np.cos(theta))
    )

    term_3 = -(S**2 / x**2) * (d / np.pi) * sumatory(4000, x_list, R_S)

    return term_1 + term_2 + term_3


def compute_phi_values(x, R_S, n_inf, S, n_IRE): 
    """
    Compute phi values for each element in the input array `x`.

    This function evaluates `phi_x` at each point in `x`, using the full array 
    and additional parameters. The results are returned as a NumPy array.

    Parameters
    ----------
    x : array_like
        1D array of input values where `phi_x` is evaluated.
    R_S : float
        A parameter required by the `phi_x` function.
    n_inf : float
        A parameter required by the `phi_x` function.
    S : float
        A parameter required by the `phi_x` function.
    n_IRE : float
        A parameter required by the `phi_x` function.

    Returns
    -------
    phi_values : np.ndarray
        Array of computed phi values corresponding to each entry in `x`.
    """
    return np.array([phi_x(xi, x, R_S, n_inf, S, n_IRE) for xi in x])


def compute_n_sample(phi, R_S, n_IRE):
    """
    Compute the complex refractive index of the sample from phase and R_S.

    Parameters:
    ----------
    n_IRE : float
        Refractive index of the internal reflection element (IRE).
    phi : array-like of float
        Phase values φ(x) in radians.
    R_S : array-like
        Reflectance values R_S at each wavenumber.

    Returns:
    -------
    np.ndarray of complex
        Complex refractive index of the sample at the given points.
    """
    phi = np.asarray(phi)

    theta = np.pi / 4  # 45 degrees incidence
    numerator = 1 + np.sqrt(R_S) * np.exp(1j * phi)
    denominator = 1 - np.sqrt(R_S) * np.exp(1j * phi)
    ratio_squared = (numerator / denominator) ** 2

    n_sample = n_IRE * np.sqrt(np.sin(theta)**2 + np.cos(theta)**2 * ratio_squared)

    return n_sample

def extract_eps_NP(exp_R, exp_x, f_s, f, eps_inf_np, S=0, n_IRE=2.4 ):
    """
    Extract the nanoparticle dielectric function (ε_np) from experimental reflectance data applying the LLL formula

    Parameters
    ----------
    exp_R : array-like
        Experimental reflectance values.
    exp_x : array-like
        Spectral axis or spatial positions corresponding to the reflectance values.
    f_s : float
        polarization degree
    f : float
        filling factor of nanoparticles in the effective medium.
    S : float, optional
        oscillator term
    n_IRE : float, optional
        Refractive index of the internal reflection element, by default 2.4.

    Returns
    -------
    np.ndarray
        Extracted complex dielectric function of the nanoparticles.
    """
    # Step 1: Estimate effective medium properties
    eps_inf_eff = (eps_inf_np**(1/3) * f + (1 - f))**3 #here we apply the LLL formula, if another EMT has to be applyied, this line shoud be modified.
    n_inf_eff = np.sqrt(eps_inf_eff)

    # Step 2: Compute phase angle and reflectance correction
    phi_rad = get_pol_angle(f_s)
    R_S = R_S_45(phi_rad, exp_R)

    # Step 3: Compute phase values
    phi_values = compute_phi_values(exp_x, R_S, n_inf_eff, S, n_IRE)

    # Step 4: Extract effective refractive index and dielectric function
    n_eff = compute_n_sample(phi_values, R_S, n_IRE)
    eps_eff = n_eff**2

    # Step 5: Retrieve nanoparticle permittivity
    eps_np = ((eps_eff**(1/3) - (1 - f)) / f)**3

    return eps_np


def obtain_reference_data():
    """
    Load and return the refractive index data of water from Segelstein et al. (1981).

    This function reads a text file containing the wavelength-dependent real and
    imaginary parts of the refractive index of water. The wavelengths are converted 
    from microns to wavenumbers in cm-1.

    Returns:
        h2o_cm (np.ndarray): Wavenumber values (in cm⁻¹).
        h2o_n_re (np.ndarray): Real part of the refractive index.
        h2o_n_im (np.ndarray): Imaginary part of the refractive index.
    """
    nname = "data/segelstein81_index.txt"
    from importlib_resources import files
    ref = files('miepython').joinpath(nname)
    h2o = np.genfromtxt(ref, delimiter='\t', skip_header=4)

    h2o_lam, h2o_n_re, h2o_n_im = h2o[:, 0], h2o[:, 1], h2o[:, 2]
    h2o_cm = 10000 / h2o_lam

    return h2o_cm, h2o_n_re, h2o_n_im