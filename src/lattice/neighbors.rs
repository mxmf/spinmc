use serde::Deserialize;
use serde::Serialize;

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum InteractionType {
    Exchange,
    DzyaloshinskiiMoriya,
    DoubleExchange,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct NeighborDescriptor {
    pub from_sub: usize,
    pub to_sub: usize,
    pub offset: [isize; 3],
    pub interaction: InteractionType,
    pub coupling: f64,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Config {
    pub neighbors: Vec<NeighborDescriptor>,
}
