use super::Structure;
use serde::Serialize;
use std::fmt;

use crate::config::Deserialize;

#[derive(Debug, Deserialize, Serialize)]
pub struct Exchange {
    #[serde(default)]
    pub from_sublattice: Option<usize>,
    #[serde(default)]
    pub to_sublattice: Option<usize>,
    #[serde(default)]
    pub offsets: Option<Vec<[isize; 3]>>,
    #[serde(default)]
    pub neighbor_order: Option<usize>,
    pub strength: f64,
}

#[derive(Debug)]
pub struct ParsedExchange {
    pub from_sub: usize,
    pub to_sub: usize,
    pub offset: [isize; 3],
    pub strength: f64,
}

impl Exchange {
    pub fn parse(
        &self,
        structure: &Option<Structure>,
        pbc: [bool; 3],
    ) -> anyhow::Result<Vec<ParsedExchange>> {
        let mut exchange_params = vec![];

        let (from_sub, to_sub, offsets, neighbor_order, strength) = (
            &self.from_sublattice,
            &self.to_sublattice,
            &self.offsets,
            &self.neighbor_order,
            &self.strength,
        );

        match (from_sub, to_sub, offsets, neighbor_order, structure) {
            (_, _, Some(_), Some(_), _) => anyhow::bail!(
                "Invalid configuration: do not set both `offsets` and `neighbor_order`; only one should be specified.",
            ),

            (_, _, None, None, _) => anyhow::bail!(
                "Missing configuration: you must specify either `offsets` or `neighbor_order`.",
            ),

            (_, _, _, Some(_), None) => anyhow::bail!(
                "Incomplete configuration: when using `neighbor_order`, `structure` must be set.",
            ),

            (from_sub, to_sub, Some(offsets), None, _) => {
                if let (Some(from_sub), Some(to_sub)) = (from_sub, to_sub) {
                    for offset in offsets {
                        exchange_params.push(ParsedExchange {
                            from_sub: *from_sub,
                            to_sub: *to_sub,
                            offset: *offset,
                            strength: *strength,
                        });
                    }
                } else {
                    anyhow::bail!(
                        "Incomplete configuration: when using `offsets`, both `from_sublattice` and `to_sublattice` must be specified."
                    )
                }
            }

            (from_sub, to_sub, None, Some(neighbor_order), Some(structure)) => {
                let atoms = crate::lattice::Atoms {
                    cell: structure.cell,
                    positions: structure.positions.clone(),
                    pbc,
                    tolerance: structure.tolerance.unwrap_or(0.0001),
                };

                let neighbors = match (from_sub, to_sub) {
                    (Some(from), Some(to)) => {
                        atoms.find_neighbors_from_to(*from, *to, *neighbor_order)
                    }
                    (Some(from), None) => atoms.find_neighbors_from(*from, *neighbor_order),
                    (None, None) => atoms.find_neighbors_all(*neighbor_order),
                    (None, Some(_)) => anyhow::bail!(
                        "Invalid configuration: `from_sublattice` must be specified when using `neighbor_order`."
                    ),
                };

                for neighbor in neighbors {
                    exchange_params.push(ParsedExchange {
                        from_sub: neighbor.from,
                        to_sub: neighbor.to,
                        offset: neighbor.offset,
                        strength: *strength,
                    });
                }
            }
        }

        Ok(exchange_params)
    }

    pub fn validate(&self, sublattices: usize) -> anyhow::Result<()> {
        if let Some(from) = self.from_sublattice {
            if from >= sublattices {
                anyhow::bail!(
                    "from_sublattice index ({from}) is out of range, must be less than sublattices count ({sublattices})"
                );
            }
        }
        if let Some(to) = self.to_sublattice {
            if to >= sublattices {
                anyhow::bail!(
                    "to_sublattice index ({to}) is out of range, must be less than sublattices count ({sublattices})"
                );
            }
        }

        Ok(())
    }
}

impl fmt::Display for ParsedExchange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (from_sub, to_sub, strength, offset) =
            (self.from_sub, self.to_sub, self.strength, self.offset);

        write!(
            f,
            "  {from_sub:<4} | {to_sub:<3} | {:>3} {:>3} {:>3}  | {strength:>8.12}",
            offset[0], offset[1], offset[2]
        )?;
        Ok(())
    }
}
