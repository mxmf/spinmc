use anyhow::Result;
use clap::Parser;
use colored::*;
use mimalloc::MiMalloc;
use tracing_subscriber::FmtSubscriber;

use spinmc::runner::run;

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Path to input file (TOML format)
    #[arg(short, long, default_value = "config.toml")]
    input: String,
}

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
fn main() -> Result<()> {
    let subscriber = FmtSubscriber::builder()
        .with_max_level(tracing::Level::INFO)
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");
    let args = Args::parse();

    let content = std::fs::read_to_string(&args.input)?;
    if let Err(e) = run(&content) {
        eprintln!("{}", format!("Error: {e}").red().bold());
        std::process::exit(1);
    }
    Ok(())
}
