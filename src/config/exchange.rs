use crate::lattice::Structure;
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
    #[serde(default)]
    pub distance_range: Option<[f64; 2]>,
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

        match (&self.offsets, self.neighbor_order, self.distance_range) {
            (Some(offsets), None, None) => {
                let (Some(from_sub), Some(to_sub)) = (self.from_sublattice, self.to_sublattice)
                else {
                    anyhow::bail!(
                        "Incomplete configuration: when using `offsets`, both `from_sublattice` and `to_sublattice` must be specified."
                    );
                };
                for offset in offsets {
                    exchange_params.push(ParsedExchange {
                        from_sub,
                        to_sub,
                        offset: *offset,
                        strength: self.strength,
                    });
                }
            }

            (None, Some(neighbor_order), None) => {
                let atoms = Self::atoms_from_structure(structure, pbc, "neighbor_order")?;

                let neighbors = match (self.from_sublattice, self.to_sublattice) {
                    (Some(from), Some(to)) => {
                        atoms.find_neighbors_from_to(from, to, neighbor_order)
                    }
                    (Some(from), None) => atoms.find_neighbors_from(from, neighbor_order),
                    (None, None) => atoms.find_neighbors_all(neighbor_order),
                    (None, Some(_)) => anyhow::bail!(
                        "Invalid configuration: `from_sublattice` must be specified when using `neighbor_order`."
                    ),
                };

                for neighbor in neighbors {
                    exchange_params.push(ParsedExchange {
                        from_sub: neighbor.from,
                        to_sub: neighbor.to,
                        offset: neighbor.offset,
                        strength: self.strength,
                    });
                }
            }

            (None, None, Some([min_distance, max_distance])) => {
                let atoms = Self::atoms_from_structure(structure, pbc, "distance_range")?;

                let distances = match (self.from_sublattice, self.to_sublattice) {
                    (Some(from), Some(to)) => {
                        atoms.calc_distance_range_from_to(from, to, min_distance, max_distance)
                    }
                    (Some(from), None) => {
                        atoms.calc_distance_range_from(from, min_distance, max_distance)
                    }
                    (None, None) => atoms.calc_distance_range_all(min_distance, max_distance),
                    (None, Some(_)) => anyhow::bail!(
                        "Invalid configuration: `from_sublattice` must be specified when using `distance_range`."
                    ),
                };

                for distance in distances {
                    exchange_params.push(ParsedExchange {
                        from_sub: distance.neighbor.from,
                        to_sub: distance.neighbor.to,
                        offset: distance.neighbor.offset,
                        strength: self.strength,
                    });
                }
            }

            (None, None, None) => anyhow::bail!(
                "Missing configuration: you must specify one of `offsets`, `neighbor_order`, or `distance_range`.",
            ),

            _ => anyhow::bail!(
                "Invalid configuration: only one of `offsets`, `neighbor_order`, or `distance_range` may be specified.",
            ),
        }

        Ok(exchange_params)
    }

    fn atoms_from_structure(
        structure: &Option<Structure>,
        pbc: [bool; 3],
        field: &str,
    ) -> anyhow::Result<crate::lattice::Atoms> {
        let Some(structure) = structure else {
            anyhow::bail!(
                "Incomplete configuration: when using `{field}`, `structure` must be set."
            );
        };
        Ok(crate::lattice::Atoms {
            cell: structure.cell,
            positions: structure.positions.clone(),
            pbc,
            tolerance: structure.tolerance.unwrap_or(0.0001),
        })
    }

    pub fn validate(&self, sublattices: usize) -> anyhow::Result<()> {
        let specified = self.offsets.is_some() as usize
            + self.neighbor_order.is_some() as usize
            + self.distance_range.is_some() as usize;
        if specified != 1 {
            anyhow::bail!(
                "exactly one of `offsets`, `neighbor_order`, or `distance_range` must be specified"
            );
        }

        if let Some([min_distance, max_distance]) = self.distance_range {
            if min_distance < 0.0 {
                anyhow::bail!("distance_range minimum ({min_distance}) must be non-negative");
            }
            if max_distance < min_distance {
                anyhow::bail!(
                    "distance_range maximum ({max_distance}) must be greater than or equal to minimum ({min_distance})"
                );
            }
        }

        if let Some(from) = self.from_sublattice
            && from >= sublattices
        {
            anyhow::bail!(
                "from_sublattice index ({from}) is out of range, must be less than sublattices count ({sublattices})"
            );
        }
        if let Some(to) = self.to_sublattice
            && to >= sublattices
        {
            anyhow::bail!(
                "to_sublattice index ({to}) is out of range, must be less than sublattices count ({sublattices})"
            );
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
