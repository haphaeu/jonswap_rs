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
};

fn main() {
    let length = 10;
    let w1 = 0.1;
    let w2 = 3.1;

    let hs = 2.5;
    let tp = 8.0;
    let wp = TWOPI / tp;
    
    let gamma = gamma(hs, tp);
    
    let w = freq_domain(w1, w2, length);
    let t = period_domain(&w);
    let pm = pierson_moskowitz(hs, tp, &w);
    let js = jonswap(wp, gamma, &w, &pm);

    let m2 = spectral_moment(2, &js, &w);
    let m4 = spectral_moment(4, &js, &w);

    let amps = jonswap_component_amplitude(&js, &w);
    let phis = phases(js.len(), 112345);
    let td = time_domain(0.0, 2.0, 0.1);
    
    println!("gamma: {}", gamma);
    println!("w: {:?}", w);
    println!("t: {:?}", t);
    println!("pm: {:?}", pm);
    println!("js: {:?}", js);
    println!("m2: {m2}");
    println!("m4: {m4}");
    println!("amps: {amps:?}");
    println!("phis: {phis:?}");
    println!("length: {}", td.len());
    println!("td: {td:?}");
    
}
