use std::env;
use std::process;

use jonswap::{
    TWOPI,
    gamma,
    jonswap,
    pierson_moskowitz,
    period_domain,
    freq_domain,
    spectral_moment,
    jonswap_component_amplitude,
    phases,
    time_domain,
    wave_elevation,
};

fn main() {

    let mut args = env::args();
    args.next();
    let hs: f64 = args.next().unwrap_or_else(|| {
	eprintln!("Problem reading Hs value");
	process::exit(1);
    }).parse().unwrap();
    let tp: f64 = args.next().unwrap_or_else(|| {
	eprintln!("Problem reading Tp value");
	process::exit(1);
    }).parse().unwrap();
	
    let wp = TWOPI / tp;    
    let length = 50;
    let w1 = 0.1;
    let w2 = 6.3;
    let gamma = gamma(hs, tp);
    
    let w = freq_domain(w1, w2, length);
    let t = period_domain(&w);
    let pm = pierson_moskowitz(hs, tp, &w);
    let js = jonswap(wp, gamma, &w, &pm);

    let m0 = spectral_moment(0, &js, &w);
    let m2 = spectral_moment(2, &js, &w);
    let m4 = spectral_moment(4, &js, &w);

    let amp = jonswap_component_amplitude(&js, &w);
    let phi = phases(js.len(), 112345);
    let td = time_domain(0.0, 60.0, 0.5);
    let eta = wave_elevation(&amp, &w, &phi, &td);
    
    println!("gamma: {}", gamma);
    println!("m0: {m0}");
    println!("m2: {m2}");
    println!("m4: {m4}");

    println!("{}\t{}\t{}\t{}\t{}\t{}", "w", "t", "pm", "js", "amp", "phi");
    for i in 0..length as usize {
	println!("{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}", w[i], t[i], pm[i], js[i], amp[i], phi[i]);
    }

    println!("{:}\t{:}", "td", "eta");
    for i in 0..td.len() as usize {
	println!("{:.3}\t{:.3}", td[i], eta[i]);
    }    
}
