use clap::{Args, Parser, Subcommand};

#[derive(Debug, Args, Clone)]
pub struct GlobalOpts {
    #[arg(short, long, default_value = "500")]
    pub width: i64,
    #[arg(short, long, default_value = "500")]
    pub height: i64,
    #[arg(long, default_value = "0.01")]
    pub step_size: f64,
    #[arg(long, default_value = "15000")]
    pub max_steps: usize,
    #[arg(long, default_value = "15000")]
    pub max_radius: f64,
    #[arg(long, default_value = "0.00001")]
    pub epsilon: f64,
    #[arg(short, long, value_delimiter = ',', default_value = "0.0,0.8,-18.0")]
    pub camera_position: Vec<f64>,
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct App {
    #[clap(flatten)]
    pub global_opts: GlobalOpts,
    #[arg(short, long)]
    pub config_file: String,
    #[command(subcommand)]
    pub action: Action,
}

#[derive(Subcommand, Debug, Clone)]
pub enum Action {
    Render {
        #[arg(long, default_value = "render.png")]
        filename: String,
    },
    RenderRay {
        #[arg(short, long)]
        row: i64,
        #[arg(short, long)]
        col: i64,
        #[arg(long, default_value = "rendered-ray.csv")]
        filename: String,
    },
    RenderRayAt {
        #[arg(short, long, value_delimiter = ',')]
        position: Vec<f64>,
        #[arg(short, long, value_delimiter = ',')]
        direction: Vec<f64>,
        #[arg(long, default_value = "rendered-ray-at.csv")]
        filename: String,
    },
}
