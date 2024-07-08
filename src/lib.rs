use std::iter::successors;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

pub const PI: f64 = 3.14159265358979323846;
pub const TWOPI: f64 = 6.28318530717958647692;

/// Jonswap spectrum
///
/// * `hs`       : [m]       significant wave height
/// * `tp`       : [s]       spectral peak period
/// * `gamma`    : [-]       JONSWAP peak enhancement factor
/// * `nharms`   : [-]       number of harmonics
/// * `seed`     : [-]       seed for random phase generator
/// * `w0`, `w1` : [rad/s]   boundaries of angular frequency domain
/// * `w`        : [rad/s]   vector with angular frequency domain for the spectrum
/// * `t`        : [s]       vector with period domain for the spectrum
/// * `pm`       : [m²s/rad] vector with Pierson-Morkowitz wave spectrum
/// * `js`       : [m²s/rad] vector with JONSWAP wave spectrum
/// * `amp`      : [m]       vector with amplitude for each harmonic
/// * `phi`      : [rad]     vector with the phases for each harmonic
/// * `td`       : [s]       vector with the time domain for wave elevation time trace
/// * `eta`      : [m]       vector with wave elevation time trace
pub struct JonswapSpectrum {
    hs: f64,      
    tp: f64,
    wp: f64,
    gamma: f64,
    nharms: usize,
    seed: u64,
    w0: f64,
    w1: f64,
    to: f64,
    tf: f64,
    ts: f64,
    w: Option<Vec<f64>>,
    t: Option<Vec<f64>>,
    pm: Option<Vec<f64>>,
    js: Option<Vec<f64>>,
    amp: Option<Vec<f64>>,
    phi: Option<Vec<f64>>,
    td: Option<Vec<f64>>,
    eta: Option<Vec<f64>>,
}
impl Default for JonswapSpectrum {
    fn default() -> Self {
	Self {
	    hs: 4.5,
	    tp: 10.0,
	    wp: TWOPI / 10.0,
	    gamma: gamma(4.5, 10.0),
	    nharms: 200,
	    seed: 12343,
	    w0: PI / 10.0,
	    w1: 2.0 * TWOPI,
	    to: 0.0,
	    tf: 120.0,
	    ts: 0.1,
	    w: None,
	    t: None,
	    pm: None,
	    js: None,
	    amp: None,
	    phi: None,
	    td: None,
	    eta: None,
	}
    }
}
impl JonswapSpectrum {
    pub fn new(hs: f64, tp: f64) -> Self {
	assert!(hs > 0.0);
	assert!(tp > 0.0);
	if !jonswap_is_valid(hs, tp) {
	    eprintln!("WARNING: Hs-Tp pair is outside of the validity of JONSWAP model.")
	}
	let wp = TWOPI / tp;
	let gamma = gamma(hs, tp);
	Self { hs, tp, wp, gamma, ..Self::default() }
    }
    pub fn update_gamma(&mut self, gamma: f64) {
	assert!(gamma > 0.0);
	self.gamma = gamma;
    }
    pub fn set_nharms(&mut self, nharms: usize) {
	self.nharms = nharms;
    }
    pub fn set_seed(&mut self, seed: u64) {
	self.seed = seed;
    }
    pub fn set_duration(&mut self, duration: f64) {
	self.tf = duration;
    }
    pub fn set_timestep(&mut self, timestep: f64) {
	self.ts = timestep;
    }
    pub fn set_w0_w1(&mut self, w0: f64, w1: f64) {
	self.w0 = w0;
	self.w1 = w1;
    }
    pub fn calculate_spectrum(&mut self) {
	self.w = Some(freq_domain(self.w0, self.w1, self.nharms));
	if let Some(w) = &self.w {
	    self.t = Some(period_domain(w));
	    self.pm = Some(pierson_moskowitz(self.hs, self.tp, w));
	    if let Some(pm) = &self.pm {
		self.js = Some(jonswap(self.wp, self.gamma, w, pm));
	    }
	    if let Some(js) = &self.js {
		self.amp = Some(jonswap_component_amplitude(js, w));
	    }
	}
	self.phi = Some(phases(self.nharms, self.seed));	       
    }
    pub fn calculate_time_realisation(&mut self) {
	
	self.td = Some(time_domain(self.to, self.tf, self.ts));
	if let (
	    Some(amp), Some(w), Some(phi), Some(td)
	) = (
	    &self.amp, &self.w, &self.phi, &self.td
	) {
	    self.eta = Some(wave_elevation(amp, w, phi, td));
	} else {
	    eprintln!("Should call `calculate_spectrum` first.");
	}
    }
    
    pub fn show_spectrum(&self, gnuplot_friendly: bool) {
	
	if let (
	    Some(t), Some(w), Some(pm), Some(js), Some(amp), Some(phi)
	) = (
	    &self.t, &self.w, &self.pm, &self.js, &self.amp, &self.phi
	) {
	    if !gnuplot_friendly {
		println!("Wave Hs            : {:8.2} m", self.hs);
		println!("Wave Tp            : {:8.2} s", self.tp);
		println!("Wave Gamma         : {:8.3} -", self.gamma);
		if let (Some(js), Some(w)) = (&self.js, &self.w) {
		    let m0 = spectral_moment(0, js, w);
		    let m2 = spectral_moment(2, js, w);
		    let m4 = spectral_moment(4, js, w);
		    let hm0 = 4.0 * m0.sqrt();
		    println!("Spectral moment m0 : {:8.4} m²", m0);
		    println!("Spectral moment m2 : {:8.4} m²(rd/s)²", m2);
		    println!("Spectral moment m4 : {:8.4} m²(rd/s)^4", m4);
		    println!("Hm0                : {:8.4} m", hm0);
		    
		}
		println!("\nSpectrum\n========");
	    }
	    println!("{:>10} {:>10} {:>12} {:>12} {:>12} {:>10}",
		     "T", "w", "PM", "JS", "amp", "phi");
	    if !gnuplot_friendly {
		println!("{:>10} {:>10} {:>12} {:>12} {:>12} {:>10}",
			 "[s]", "[rd/s]", "[m²s/rd]", "[m²s/rd]", "[m]", "[rd]");
	    }
	    for i in 0..self.nharms as usize {
		println!("{:10.3} {:10.6} {:12.8} {:12.8} {:12.5} {:10.5}",
			 t[i], w[i], pm[i], js[i], amp[i], phi[i]);
	    }
	} else {
	    eprintln!("Should call `calculate_spectrum` first.");
	}
    }
    
    pub fn show_time_realisation(&self, gnuplot_friendly: bool) {
	
	if let (Some(td), Some(eta)) = (&self.td, &self.eta) {
	    if !gnuplot_friendly {
		println!("\nTime History\n============");
	    }
            println!("{:10} {:10}", "Time", "Elevation");
            for j in 0..td.len() {
		println!("{:10.3} {:10.6}", td[j], eta[j]);
	    }
	} else {
	    eprintln!("Should call `calculate_time_realisation` first.");
	}
    }
}



/// Return gamma factor for a JONSWAP spectrum according to DNV
pub fn gamma(hs: f64, tp: f64) -> f64 {
    let chk = tp / hs.sqrt();
    if chk < 3.6 {
	5.0
    } else if chk > 5.0 {
	1.0
    } else {
	(5.75 - 1.15 * chk).exp()
    }
}

/// Return `true` is `hs` and `tp` are within a valid range for a JONSWAP spectrum
pub fn jonswap_is_valid(hs: f64, tp: f64) -> bool {
    let ratio = tp / hs.sqrt();
    3.6 < ratio && ratio < 5.0
}

/// Return an array of angular frequencies [rd/s] from w1 to w2, of size length.
pub fn freq_domain(w1: f64, w2: f64, length: usize) -> Vec<f64> {

    // Geometric progression f[i+1] = r * f[i]
    // Better discretisation in the region of interest.
    let r = (w2 / w1).powf(1.0 / f64::from(length as i32- 1));
    successors(Some(w1), |&w| {
	if w < w2 {
	    Some(w * r)
	} else {
	    None
	}
    }).collect()
}

/// Returns an array of periods [s] corresponding to the angular frequencies from the array w.
pub fn period_domain(w: &Vec<f64>) -> Vec<f64> {
    w.iter().map(|w| TWOPI / w).collect()
}

/// Return Pierson Moskowitz spectrum.
///
/// * `hs`     : significant wave height
/// * `tp`     : peak period
/// * `w`      : angular frequencies, or NULL to calculate it.
///
/// Note that the spectra are being calculated with radians, so
/// the dimension is [m**2 * s / rad]
///
pub fn pierson_moskowitz(hs: f64, tp: f64, w: &Vec<f64>) -> Vec<f64> {
    let wp = TWOPI / tp;
    let pm_cte = 0.3125 * hs.powi(2) * wp.powi(4);
    w
	.iter()
	.map(|w| pm_cte * w.powi(-5) * (-1.25 * (w/wp).powi(-4)).exp())
	.collect()
}

/// Return JONSWAP spectrum.
///
/// * `wp`     : spectral peak angular frequency
/// * `gamma`  : gamma parameter
/// * `w`      : angular frequencies
/// * `pm`     : Pierson Moskowitz spectrum
///
/// Note that the spectra are being calculated with radians, so
/// the dimension is [m**2 * s / rad]
///
pub fn jonswap(wp: f64, gamma: f64, w: &Vec<f64>, pm: &Vec<f64>) -> Vec<f64> {
    let norm = 1.0 - 0.287 * gamma.ln();
    pm
	.iter()
	.zip(w)
	.map( |(pi, &wi)| {
	    let sigma = if wi < wp { 0.07 } else { 0.09 };
	    norm * pi * gamma.powf(
		(-0.5 * ((wi - wp) / wp / sigma).powi(2)).exp()
	    )
	})
	.collect()
}

/// Return the n-th spectral moment.
/// 
/// * `n`      : order of the spectral moment to be calculated
/// * `s`      : spectrum array
/// * `w`      : angular frequency array
///
/// Note that the spectra are being calculated with radians, so
/// the dimension of the n-th moment is [(rad / s)**n * m**2]
///
pub fn spectral_moment(n: i32, s: &Vec<f64>, w: &Vec<f64>) -> f64 {
    let mut m: f64 = 0.0;
    // TODO: is there an iterator to replace this loop?
    for i in 1..s.len() {
	m += ((w[i] + w[i-1]) / 2.0).powi(n) * (w[i] - w[i-1]) * (s[i] + s[i-1]) / 2.0;
    }
    m
}

/// Return the amplitude of a sinusoidal signal for each frequency component
/// of a discretizes spectrum.
///     
/// * `s`      : spectrum array
/// * `w*       : angular frequency array
/// 
/// Ref: https://ceprofs.civil.tamu.edu/jzhang/ocen300/statistics-spectrum.pdf
///
pub fn jonswap_component_amplitude(s: &Vec<f64>, w: &Vec<f64>) -> Vec<f64> {
    let mut amp = Vec::with_capacity(s.len());
    for i in 1..s.len() {
	amp.push(
	    (2.0 * s[i] * (w[i] - w[i - 1])).sqrt()
	);
    }
    amp.push(0.0);
    amp
}

/// Return randomised phases between -PI and PI
///
pub fn phases(length: usize, seed: u64) -> Vec<f64> {
    let mut phi: Vec<f64> = Vec::with_capacity(length);
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..length {
	phi.push( PI * (2.0 * rng.gen::<f64>() - 1.0) );
    }
    phi
}


/// Return an array of times domain [s] to be used to calculate the wave elevation.
/// 
/// * `to`     : start time
/// * `tf`     : final time
/// * `ts`     : time step
///
pub fn time_domain(to: f64, tf: f64, ts: f64) -> Vec<f64> {
    successors(Some(to), |&t| {
	if (t - tf).abs() < ts / 2.0 {
	    None
	} else {
	    Some(t + ts)
	}
    }).collect()
}

/// Return an array of times domain [s] to be used to calculate the wave elevation.
///
/// Note that since input is time-step, the array length must be passed
/// as a pointer and it is where the length will be stored.
///
pub fn wave_elevation(amp: &Vec<f64>, w: &Vec<f64>, phi: &Vec<f64>, time: &Vec<f64>) -> Vec<f64> {
    time.iter().map(
	|&ti|
	amp.iter().zip(w).zip(phi).map(
	    |((&ampj, &wj), &phij)|
	    ampj * (wj * ti + phij).cos()
	).sum()
    ).collect()
}

// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////    U N I T    T E S T S  //////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    //use float_cmp::{approx_eq, assert_approx_eq};

    const ATOL: f64 = 1.0e-12;

    // frequency domain usedfor various tets
    // generated with `freq_domain(0.1, 3.1, 10);`
    const ANG_FREQ_DOMAIN: [f64; 10] = [
	    0.1, 0.1464558942252767, 0.2144932895332543, 0.31413806523913934,
	    0.46007371254796453, 0.6738050698075504, 0.9868272403218976,
	    1.445266659272055, 2.116678209776671, 3.100
    ];

    // Pierson-Moskowitz spectrum used in various tests
    // Generated with `pierson_moskowitz(2.5, 8.0, &ANG_FREQ_DOMAIN);`
    const PM_SPECTRUM: [f64; 10] = [
	0.0, 0.0, 4.2212229443328884e-95, 1.492713582894374e-19,
	8.840837064858753e-04, 5.324821611024473e-01, 4.809147991546064e-01,
	1.056832048498532e-01, 1.7081457182245e-02, 2.582524193397351e-03
    ];

    /// Return `true` is `abs(x - y) < atol`
    fn approx_eq(x: f64, y: f64, atol: f64) -> bool {
	println!("{x:e} vs {y:e}");
	(x - y).abs() < atol
    }

    /// Return `true` if both vectors have the same length and
    /// `abs(ui - vi) < atol` for all elements.
    fn vec_compare(u: &Vec<f64>, v: &Vec<f64>, atol: f64) -> bool {
	(u.len() == v.len() ) &&
	 u
	 .iter()
	 .zip(v)
	 .all(|(x, y)| approx_eq(*x, *y, atol))
	}

    #[test]
    fn test_js_valid() {
	assert_eq!(false, jonswap_is_valid(2.5, 5.5));
	assert_eq!(false, jonswap_is_valid(2.5, 8.0));
	assert_eq!(true, jonswap_is_valid(2.5, 6.0));
	assert_eq!(true, jonswap_is_valid(2.5, 7.5));

    }
    
    #[test]
    fn test_gamma() {
        assert!(approx_eq(5.0, gamma(1.0, 1.0), ATOL));
	assert!(approx_eq(1.0, gamma(1.0, 6.0), ATOL));
	let x1 = 3.1581929096897676272507006280068;
	let x2 = gamma(4.0, 8.0);
	assert!(approx_eq(x1, x2, ATOL));
    }

    #[test]
    fn test_freq_domain() {
	assert!(vec_compare(
	    &ANG_FREQ_DOMAIN.to_vec(),
	    &freq_domain(0.1, 3.1, 10),
	    ATOL
	));	 
    }

    #[test]
    fn test_period_domain() {
	let t1 = period_domain(&vec![TWOPI, PI]);
	let t2 = vec![1.0, 2.0];
	assert!(vec_compare(&t1, &t2, ATOL));
    }

    #[test]
    fn test_pm() {
	let pm2 = pierson_moskowitz(2.5, 8.0, &ANG_FREQ_DOMAIN.to_vec());
	assert!(vec_compare(&PM_SPECTRUM.to_vec(), &pm2, ATOL));
    }

    #[test]
    fn test_js_eq_pm() {
	let hs = 2.5;
	let tp = 8.0;
	let wp = TWOPI / tp;
	let g = gamma(hs, tp);
	let pm = pierson_moskowitz(hs, tp, &ANG_FREQ_DOMAIN.to_vec());
	let js = jonswap(wp, g, &ANG_FREQ_DOMAIN.to_vec(), &pm);
	assert_eq!(1.0, g);
	assert!(vec_compare(&pm, &js, ATOL));
    }

    #[test]
    fn test_js() {
	let hs = 3.5;
	let tp = 7.5;
	let wp = TWOPI / tp;
	let g = gamma(hs, tp);
	let pm = pierson_moskowitz(hs, tp, &ANG_FREQ_DOMAIN.to_vec());
	let js1 = jonswap(wp, g, &ANG_FREQ_DOMAIN.to_vec(), &pm);
	let js2 = vec![
	    0.0, 0.0, 1.300718043606472e-123, 1.4410938828945305e-25,
	    6.619940674687616e-5, 0.4714107603150374, 0.8324180119808562,
	    0.17473456030720524, 0.028960852623827464, 0.004402536499113709
	];
	assert!(jonswap_is_valid(hs, tp));	
	assert!(vec_compare(&js1, &js2, ATOL));
    }

    #[test]
    fn test_moments() {
	let hs = 3.5;
	let tp = 7.5;
	let wp = TWOPI / tp;
	let g = gamma(hs, tp);
    	let w = freq_domain(0.1, 6.3, 1_000);
	let pm = pierson_moskowitz(hs, tp, &w);
	let js = jonswap(wp, g, &w, &pm);
	let m0 = spectral_moment(0, &js, &w);
	let m2 = spectral_moment(2, &js, &w);
	let m4 = spectral_moment(4, &js, &w);
	assert!(jonswap_is_valid(hs, tp));	
	assert!(approx_eq(m0, 0.7669760588544421, ATOL));
	assert!(approx_eq(m2, 0.8830789568568432, ATOL));
	assert!(approx_eq(m4, 2.4409619089802894, ATOL));
    }

    #[test]
    fn test_phases() {
	let phi: Vec<f64> = phases(1_000, 1123);
	let min = phi.iter().fold( PI + 1.0, |a, &b| a.min(b));
	let max = phi.iter().fold(-PI - 1.0, |a, &b| a.max(b));
	assert_eq!(1_000, phi.len());
	assert!(max < PI);
	assert!(min > -PI);
    }
    
    #[test]
    fn test_time_domain() {
	let td1 = time_domain(0.0, 1.0, 0.1);
	let td2 =vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
	println!("{:?}", td1);
	assert_eq!(11, td1.len());
	assert!(vec_compare(&td1, &td2, ATOL));
    }
    
    #[test]
    fn test_amplitudes() {
	assert!(false);
    }

    #[test]
    fn test_wave_elevation() {
	assert!(false);
    }

}
