use anyhow::Context;
use chemfiles;
use std::io::Read;

#[derive(Debug)]
pub struct Structure {
    pub cell: [[f64; 3]; 3],
    pub positions: Vec<[f64; 3]>,
    pub tolerance: Option<f64>,
    pub magnetic_indices: Option<Vec<u64>>,
    pub frame: Option<chemfiles::Frame>,
}

fn is_vasp_format(format: &Option<String>, filename: &str) -> bool {
    if format.as_ref().is_some_and(|f| {
        matches!(
            f.to_ascii_lowercase().as_str(),
            "poscar" | "vasp" | "contcar"
        )
    }) {
        return true;
    }

    if format.is_none() {
        use std::path::Path;
        let path = Path::new(filename);
        if path
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("vasp"))
        {
            return true;
        }
        if let Some(stem) = path.file_stem().or_else(|| path.file_name()) {
            let stem = stem.to_string_lossy();
            if matches!(stem.to_ascii_lowercase().as_str(), "poscar" | "contcar") {
                return true;
            }
        }
    }

    false
}

pub fn load_from_file(stru_file: &str, format: Option<String>) -> anyhow::Result<Structure> {
    if is_vasp_format(&format, stru_file) {
        let mut file = std::fs::File::open(stru_file)
            .with_context(|| format!("failed to open POSCAR `{stru_file}`"))?;
        let mut content = String::new();
        file.read_to_string(&mut content)
            .with_context(|| format!("failed to read POSCAR `{stru_file}`"))?;
        return parse_poscar_from_str(&content, stru_file);
    }

    let mut trajectory = if let Some(format) = &format {
        chemfiles::Trajectory::open_with_format(stru_file, 'r', format.as_str()).with_context(
            || format!("failed to read structure from structure file `{stru_file}`"),
        )?
    } else {
        chemfiles::Trajectory::open(stru_file, 'r').with_context(|| {
            format!("failed to read structure from structure file `{stru_file}`")
        })?
    };

    let mut frame = chemfiles::Frame::new();
    trajectory
        .read(&mut frame)
        .with_context(|| format!("failed to read structure from structure file `{stru_file}`"))?;

    Ok(Structure {
        positions: frame.positions().to_vec(),
        cell: frame.cell().matrix(),
        tolerance: None,
        magnetic_indices: None,
        frame: Some(frame),
    })
}

fn parse_poscar_from_str(content: &str, label: &str) -> anyhow::Result<Structure> {
    let lines: Vec<&str> = content.lines().collect();

    if lines.len() < 7 {
        anyhow::bail!("POSCAR `{label}` is too short, expected at least 7 lines");
    }

    let scale_parts: Vec<f64> = lines[1]
        .split_whitespace()
        .map(|s| s.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()
        .with_context(|| format!("invalid scale factor in POSCAR `{label}`"))?;

    // NOTE: A negative scale factor interpreted as target cell volume (VASP convention)
    // is NOT supported. This is a rarely-used feature and adds complexity (requires
    // computing volume and solving for scale). Users should provide positive scale
    // or set scale=1 and give lattice vectors directly.
    let scale: [f64; 3] = match scale_parts.len() {
        1 => [scale_parts[0], scale_parts[0], scale_parts[0]],
        3 => [scale_parts[0], scale_parts[1], scale_parts[2]],
        n => anyhow::bail!("scale factor must be 1 or 3 numbers, got {n} in POSCAR `{label}`"),
    };

    let mut lattice: [[f64; 3]; 3] = [[0.0; 3]; 3];
    for i in 0..3 {
        let parts: Vec<&str> = lines[2 + i].split_whitespace().collect();
        if parts.len() != 3 {
            anyhow::bail!("invalid lattice vector line {} in POSCAR `{label}`", i + 1);
        }
        for j in 0..3 {
            lattice[i][j] = parts[j]
                .parse::<f64>()
                .with_context(|| format!("invalid lattice vector in POSCAR `{label}`"))?
                * scale[j];
        }
    }

    // Line 5: element names or atom counts
    let (_group_symbols, group_counts, coord_start) = if lines[5]
        .trim()
        .chars()
        .next()
        .is_some_and(|c| c.is_ascii_digit())
    {
        (vec![], parse_group_counts(lines[5], label)?, 6)
    } else {
        let symbols: Vec<String> = lines[5].split_whitespace().map(|s| s.to_string()).collect();
        (symbols, parse_group_counts(lines[6], label)?, 7)
    };

    if coord_start >= lines.len() {
        anyhow::bail!("POSCAR `{label}` missing coordinate data");
    }

    let coord_mode_line = lines[coord_start].trim();
    let first_char = coord_mode_line
        .chars()
        .next()
        .map(|c| c.to_ascii_uppercase());

    let (is_direct, positions_start) = match first_char {
        Some('S') => {
            let mode_line = lines
                .get(coord_start + 1)
                .map(|l| l.trim())
                .unwrap_or("")
                .chars()
                .next()
                .map(|c| c.to_ascii_uppercase());
            let is_direct = !matches!(mode_line, Some('C') | Some('K'));
            (is_direct, coord_start + 2)
        }
        Some('C') | Some('K') => (false, coord_start + 1),
        _ => (true, coord_start + 1),
    };

    let total_atoms: usize = group_counts.iter().sum();
    if positions_start + total_atoms > lines.len() {
        anyhow::bail!("POSCAR `{label}` has insufficient coordinate lines");
    }

    let mut positions = Vec::with_capacity(total_atoms);
    let range = positions_start..positions_start + total_atoms;
    for i in range {
        let parts: Vec<&str> = lines[i].split_whitespace().collect();
        if parts.len() < 3 {
            anyhow::bail!("invalid coordinate line {} in POSCAR `{label}`", i + 1);
        }
        let x: f64 = parts[0]
            .parse()
            .with_context(|| format!("invalid coordinate in POSCAR `{label}`"))?;
        let y: f64 = parts[1]
            .parse()
            .with_context(|| format!("invalid coordinate in POSCAR `{label}`"))?;
        let z: f64 = parts[2]
            .parse()
            .with_context(|| format!("invalid coordinate in POSCAR `{label}`"))?;

        if is_direct {
            let cart = [
                x * lattice[0][0] + y * lattice[1][0] + z * lattice[2][0],
                x * lattice[0][1] + y * lattice[1][1] + z * lattice[2][1],
                x * lattice[0][2] + y * lattice[1][2] + z * lattice[2][2],
            ];
            positions.push(cart);
        } else {
            positions.push([x * scale[0], y * scale[1], z * scale[2]]);
        }
    }

    Ok(Structure {
        cell: lattice,
        positions,
        tolerance: None,
        magnetic_indices: None,
        frame: None,
    })
}

fn parse_group_counts(line: &str, label: &str) -> anyhow::Result<Vec<usize>> {
    line.split_whitespace()
        .map(|s| {
            s.parse::<usize>()
                .with_context(|| format!("invalid atom count in POSCAR `{label}`: '{s}'"))
        })
        .collect()
}

#[cfg(test)]
#[path = "structure_tests.rs"]
mod tests;
