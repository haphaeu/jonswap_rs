use std::env;
use std::str::FromStr;
use std::fmt::Debug;
use std::process;

use jonswap::JonswapSpectrum;

fn help() {
    eprintln!("Use:");
    eprintln!(
	"  {} hs tp [-n nharms] [-g gamma] [-s seed] [-d duration] [-t timestep] [-p{{s|t}}] [-h]",
	env::args().collect::<Vec<_>>()[0]
    );
}

fn help_full() {
    help();
    println!("where:");
    println!("  hs      : significant wave height in meters.");
    println!("  tp      : spectral peak period in seconds");
    println!("  -n      : number of harmonics used to discretise the spectrum.");
    println!("  -g      : value for gamma. If ommited, DNV is used.");
    println!("  -s      : seed number for phase randomisation.");
    println!("  -d      : duration - timetrace will be shown.");
    println!("  -t      : time step. default is 0.1 seconds.");
    println!("  -p{{s|t}} : gnuplot friendly output of spectrum (-ps) or timetrace (-pt),");
    println!("            to be piped to `gnuplot -p -e \"plot '-' using i:j w l\"` where");
    println!("            `i:j` are the index of the columns to be used as x and y axis.");
    println!("            Note this only works for 1 plot. For multiple plots, need to");
    println!("            save the output to a temporary file.");
    println!("  -h      : show this help.");
}

// Parse a string into a value, eg, `"2.5" => 2.5`.
// Exit process in case of parse error.
fn arg_parse<T: FromStr>(arg: &str, name: &str) -> T
where <T as FromStr>::Err: Debug  {
    match arg.parse() {
	Ok(val) => val,
	Err(e) => {
	    eprintln!("Error parsing {name} - {e:?}, got \"{arg}\"");
	    help();
	    process::exit(1);
	}
    }
}

// Call `arg_parse()` for the `.next()` item of the iteator `iter`.
fn parse_next<'a, T: FromStr>(
    iter: &mut impl Iterator<Item = &'a String>,
    name: &str
) -> T where <T as FromStr>::Err: Debug  {
    match iter.next() {
	Some(arg) => arg_parse(arg, name),
	None => {
	eprintln!("Error: incomplete argument for {name}");
	process::exit(1);
	}
    }
}

fn main() {

    let args: Vec<_> = env::args().collect();

    // Show help message and exit, if it is evoqued
    // anywhere amongs the parameters' list
    if args.contains(&String::from("-h")) {
	help_full();
	process::exit(0);
    }

    // Minimal use: `jonswap hs tp` - 3 arguments
    // Hs and Tp - mandatory arguments
    if args.len() < 3 {
	eprintln!("Error: Invalid number of arguments.");
	help();
	process::exit(1);
    }

    let hs = arg_parse::<f64>(&args[1], "Hs");
    let tp = arg_parse::<f64>(&args[2], "Tp");    
    let mut js = JonswapSpectrum::new(hs, tp);

    // Iterate and parse the other arguments
    let mut gnuplot_spectrum = false;
    let mut gnuplot_timetrace = false;
    let mut iter = args[3..].iter();
    
    while let Some(arg) = iter.next() {
	match arg.as_str() {
	    "-n" => js.set_nharms(parse_next::<usize>(&mut iter, "nharms")),
	    "-g" => js.update_gamma(parse_next::<f64>(&mut iter, "gamma")),
	    "-s" => js.set_seed(parse_next::<u64>(&mut iter, "seed")),
	    "-d" => js.set_duration(parse_next::<f64>(&mut iter, "duration")),
	    "-t" => js.set_timestep(parse_next::<f64>(&mut iter, "timetstep")),
	    "-ps" => gnuplot_spectrum = true,
	    "-pt" => gnuplot_timetrace = true,
	    inv => {
		eprintln!("Invalid argument {:?}", inv);
		help();
		process::exit(1);
	    }
	};
    }

    js.calculate_spectrum();
    if !gnuplot_timetrace {
	js.show_spectrum(gnuplot_spectrum);
    }

    if !gnuplot_spectrum {
	js.calculate_time_realisation();
	js.show_time_realisation(gnuplot_timetrace);
    }

}
