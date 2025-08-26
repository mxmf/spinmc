use anyhow::Result;
use clap::{Parser, Subcommand};
use mimalloc::MiMalloc;
use tracing_subscriber::FmtSubscriber;

#[derive(Parser, Debug)]
#[command(author, version, about, arg_required_else_help = true)]
struct Args {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Run the simulation
    Run {
        /// Path to input file (TOML format)
        #[arg(short, long, default_value = "config.toml")]
        input: String,
    },
}
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
fn main() -> Result<()> {
    use colored::*;
    use spinmc::runner::run;

    let subscriber = FmtSubscriber::builder()
        .with_max_level(tracing::Level::INFO)
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");

    let args = Args::parse();
    if let Some(Commands::Run { input }) = &args.command {
        let content = std::fs::read_to_string(input)?;
        if let Err(e) = run(&content) {
            eprintln!("{}", format!("Error: {e}").red().bold());
            std::process::exit(1);
        }
    }

    Ok(())
}
