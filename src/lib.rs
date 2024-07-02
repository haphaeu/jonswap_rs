use std::iter::successors;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;


pub const PI: f64 = 3.14159265358979323846;
pub const TWOPI: f64 = 6.28318530717958647692;


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
pub fn freq_domain(w1: f64, w2: f64, length: i32) -> Vec<f64> {

    // Geometric progression f[i+1] = r * f[i]
    // Better discretisation in the region of interest.
    let r = (w2 / w1).powf(1.0 / f64::from(length - 1));
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
/// hs     : significant wave height
/// tp     : peak period
/// length : length of the array
/// w      : angular frequencies, or NULL to calculate it.
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
/// hs     : significant wave height
/// tp     : peak period
/// fgamma : function(hs, tp) retuning gamma parameter
/// length : length of the array
/// w      : angular frequencies
/// pm     : Pierson Moskowitz spectrum
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
/// n      : order of the spectral moment to be calculated
/// s      : spectrum array
/// w      : angular frequency array
/// length : length of the s and w arrays
pub fn spectral_moment(n: i32, s: &Vec<f64>, w: &Vec<f64>) -> f64 {
    let mut m: f64 = 0.0;
    // TODO: is there an iterator to replace this loop?
    for i in 1..s.len() {
	m += ((w[i] - w[i-1]) / 2.0).powi(n) * (w[i] - w[i-1]) * (s[i] + s[i-1]);
    }
    m / 2.0
}

/// Return the amplitude of a sinusoidal signal for each frequency component
/// of a discretizes spectrum.
///     
/// s      : spectrum array
/// w      : angular frequency array
/// length : length of the arrays
/// 
/// Ref: https://ceprofs.civil.tamu.edu/jzhang/ocen300/statistics-spectrum.pdf
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
/// to     : start time
/// tf     : final time
/// ts     : time step
pub fn time_domain(to: f64, tf: f64, ts: f64) -> Vec<f64> {
    successors(Some(to), |&t| {
	if t <= tf {
	    Some(t + ts)
	} else {
	    None
	}
    }).collect()
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
    fn test_moments() {
	let m2 = spectral_moment(2, &PM_SPECTRUM.to_vec(), &ANG_FREQ_DOMAIN.to_vec());
	assert!(false);
    }

    #[test]
    fn test_phases() {
	let phis: Vec<f64> = phases(10, 1123);
	assert!(false);
    }
    
    #[test]
    fn test_time_domain() {
	assert!(false);
    }
    
    #[test]
    fn test_amplitudes() {
	assert!(false);
    }

}
