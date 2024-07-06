use std::env;
use std::str::FromStr;
use std::process;

use jonswap::{
    TWOPI,
    gamma,
    JonswapSpectrum
};

fn help() {
    eprintln!("Use:");
    eprintln!(
	"  {} hs tp [-n nharms] [-g gamma] [-s seed] [-d duration] [-t timestep] [-h]",
	env::args().collect::<Vec<_>>()[0]
    );
}
fn help_full() {
    help();
    println!("where:");
    println!("  hs: significant wave height in meters.");
    println!("  tp: spectral peak period in seconds");
    println!("  -n: number of harmonics used to discretise the spectrum.");
    println!("  -g: value for gamma. If ommited, DNV is used.");
    println!("  -s: seed number for phase randomisation.");
    println!("  -d: duration - timetrace will be shown.");
    println!("  -t: time step. default is 0.1 seconds.");	
    println!("  -h: show this help.");
}

fn arg_parse<T: FromStr>(arg: &str, name: &str) -> T {
	match arg.parse() {
		Ok(n) => n,
		Err(_) => {
			eprintln!("Error parsing {name}");
			help();
			process::exit(1);
		}
	}
}

fn main() {

    let mut args: Vec<_> = env::args().collect();
    
    if args.len() < 3 {
	eprintln!("Error: Invalid number of arguments.");
	help();
	process::exit(1);
    }

    // Hs and Tp - mandatory arguments
    let hs: f64 = match args[1].parse() {
	Ok(_hs) => _hs,
	Err(_) => {
	    eprintln!("Error parsing Hs value.");
	    help();
	    process::exit(1);
	}
    };
    let tp: f64 = match args[2].parse() {
	Ok(_tp) => _tp,
	Err(_) => {
	    eprintln!("Error parsing Tp value.");
	    help();
	    process::exit(1);
	}
    };

    let mut js = JonswapSpectrum::new(hs, tp);

    let mut iter = args[3..].iter();
    while let Some(arg) = iter.next() {
//    for arg in iter {
	match arg.as_str() {
	    "-h" => {
		help_full();
		return;
	    },
	    "-n" => js.set_nharms(arg_parse::<usize>(iter.next().expect("error"), "nharms")),
	    "-g" => js.update_gamma(arg_parse::<f64>(iter.next().expect("error"), "gamma")),
	    "-s" => js.set_seed(arg_parse::<u64>(iter.next().expect("error"), "seed")),
	    "-d" => js.set_duration(arg_parse::<f64>(iter.next().expect("error"), "duration")),
	    "-t" => js.set_timestep(arg_parse::<f64>(iter.next().expect("error"), "timetstep")),
	    inv => {
		eprintln!("Invalid argument {:?}", inv);
		help();
		process::exit(1);
	    }
	};
    }

    js.calculate_spectrum();
    js.show_spectrum();

    js.calculate_time_realisation();
    js.show_time_realisation();

}
