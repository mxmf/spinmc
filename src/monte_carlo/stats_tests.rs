use super::*;
use crate::spin::IsingSpin;

// --- maybe() ---

#[test]
fn maybe_true() {
    assert_eq!(maybe(true, || 42), Some(42));
}

#[test]
fn maybe_false() {
    assert_eq!(maybe(false, || 42), None);
}

// --- StatsConfig Display ---

#[test]
fn stats_config_display_no_flags() {
    let cfg = StatsConfig {
        energy: false,
        heat_capacity: false,
        magnetization: false,
        susceptibility: false,
        magnetization_abs: false,
        susceptibility_abs: false,
        group_magnetization: false,
        group_susceptibility: false,
        group_magnetization_abs: false,
        group_susceptibility_abs: false,
        group_num: 0,
    };
    let s = format!("{cfg}");
    assert!(s.contains("#T(K)"));
}

#[test]
fn stats_config_display_with_energy() {
    let cfg = StatsConfig {
        energy: true,
        heat_capacity: true,
        magnetization: false,
        susceptibility: false,
        magnetization_abs: false,
        susceptibility_abs: false,
        group_magnetization: false,
        group_susceptibility: false,
        group_magnetization_abs: false,
        group_susceptibility_abs: false,
        group_num: 0,
    };
    let s = format!("{cfg}");
    assert!(s.contains("Energy"));
    assert!(s.contains("C"));
}

// --- Stats::result() formulas ---

fn make_stats(overrides: impl FnOnce(&mut Stats<IsingSpin>)) -> Stats<IsingSpin> {
    let mut s = Stats {
        energy_sum: 0.0,
        energy2_sum: 0.0,
        m_sum: IsingSpin::zero(),
        m_2_sum: 0.0,
        m_abs_sum: 0.0,
        steps: 1,
        size: 4.0,
        kb: 1.0,
        t: 2.0,
        stats_config: StatsConfig {
            energy: true,
            heat_capacity: true,
            magnetization: true,
            susceptibility: true,
            magnetization_abs: true,
            susceptibility_abs: true,
            group_magnetization: false,
            group_susceptibility: false,
            group_magnetization_abs: false,
            group_susceptibility_abs: false,
            group_num: 0,
        },
        partial_m_sum: vec![],
        partial_m_2_sum: vec![],
        partial_m_abs_sum: vec![],
        partial_size: vec![],
    };
    overrides(&mut s);
    s
}

#[test]
fn result_energy() {
    let stats = make_stats(|s| {
        s.energy_sum = 10.0;
        s.steps = 2;
    });
    let r = stats.result();
    // E = sum / n / size = 10 / 2 / 4 = 1.25
    assert!((r.energy.unwrap() - 1.25).abs() < 1e-10);
}

#[test]
fn result_specific_heat() {
    let stats = make_stats(|s| {
        s.energy_sum = 10.0;
        s.energy2_sum = 100.0;
        s.steps = 2;
        s.kb = 1.0;
        s.t = 1.0;
    });
    let r = stats.result();
    // e_avg = 10/2 = 5, e2_avg = 100/2 = 50
    // Cv = (50 - 25) / (1*1) / 4 = 25/4 = 6.25
    assert!((r.specific_heat.unwrap() - 6.25).abs() < 1e-10);
}

#[test]
fn result_magnetization() {
    let stats = make_stats(|s| {
        s.m_sum = IsingSpin::along_z(8.0).unwrap();
        s.steps = 2;
    });
    let r = stats.result();
    // |m_sum/n| / size = |8/2| / 4 = 1.0
    assert!((r.magnetization.unwrap() - 1.0).abs() < 1e-10);
}

#[test]
fn result_susceptibility() {
    let stats = make_stats(|s| {
        s.m_sum = IsingSpin::along_z(8.0).unwrap();
        s.m_2_sum = 64.0;
        s.steps = 2;
        s.kb = 1.0;
        s.t = 1.0;
    });
    let r = stats.result();
    // m_avg = 8/2 = 4, m_avg.norm_sqr() = 16
    // χ = (64/2 - 16) / (1*1) / 4 = (32-16)/4 = 4.0
    assert!((r.susceptibility.unwrap() - 4.0).abs() < 1e-10);
}

#[test]
fn result_magnetization_abs() {
    let stats = make_stats(|s| {
        s.m_abs_sum = 6.0;
        s.steps = 2;
    });
    let r = stats.result();
    assert!((r.magnetization_abs.unwrap() - 6.0 / 2.0 / 4.0).abs() < 1e-10);
}

#[test]
fn result_susceptibility_abs() {
    let stats = make_stats(|s| {
        s.m_abs_sum = 6.0;
        s.m_2_sum = 100.0;
        s.steps = 2;
        s.kb = 1.0;
        s.t = 1.0;
    });
    let r = stats.result();
    // m_abs_avg = 3, χ_abs = (50 - 9) / 1 / 4 = 10.25
    assert!((r.susceptibility_abs.unwrap() - 10.25).abs() < 1e-10);
}

#[test]
fn result_disabled_observables_return_none() {
    let stats = make_stats(|s| {
        s.stats_config.energy = false;
        s.stats_config.heat_capacity = false;
        s.stats_config.magnetization = false;
        s.stats_config.susceptibility = false;
        s.stats_config.magnetization_abs = false;
        s.stats_config.susceptibility_abs = false;
    });
    let r = stats.result();
    assert!(r.energy.is_none());
    assert!(r.specific_heat.is_none());
    assert!(r.magnetization.is_none());
    assert!(r.susceptibility.is_none());
    assert!(r.magnetization_abs.is_none());
    assert!(r.susceptibility_abs.is_none());
}

#[test]
fn result_single_step_variance_zero() {
    let stats = make_stats(|s| {
        s.energy_sum = 5.0;
        s.energy2_sum = 25.0;
        s.m_sum = IsingSpin::along_z(4.0).unwrap();
        s.m_2_sum = 16.0;
        s.steps = 1;
    });
    let r = stats.result();
    // Cv should be 0 (variance = 0)
    assert!((r.specific_heat.unwrap() - 0.0).abs() < 1e-10);
    // χ should be 0 (variance = 0)
    assert!((r.susceptibility.unwrap() - 0.0).abs() < 1e-10);
}

// --- StatResult Display ---

#[test]
fn stat_result_display_contains_temperature() {
    let r = StatResult {
        t: 300.0,
        energy: Some(-42.0),
        ..Default::default()
    };
    let s = format!("{r}");
    assert!(s.contains("300"));
    assert!(s.contains("-42"));
}
