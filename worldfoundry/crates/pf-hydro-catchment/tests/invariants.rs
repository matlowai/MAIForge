use pf_hydro_catchment::*;
fn line(bed: &[f64], source: &[f64]) -> Mesh {
    let n = bed.len();
    let mut edges = vec![vec![NONE; 8]; n];
    for i in 0..n {
        if i > 0 {
            edges[i][0] = i - 1;
        }
        if i + 1 < n {
            edges[i][1] = i + 1;
        }
    }
    Mesh {
        bed: bed.to_vec(),
        source: source.to_vec(),
        point: (0..n).map(|i| [i as f64 * 0.001, 0.]).collect(),
        cell: (0..n).map(|i| i / 2).collect(),
        edges,
        fixed: vec![None; n],
        physical_len: n,
        radius: 6_371_000.,
    }
}
#[test]
fn joint_boundary_fills_pit_without_mutating_bed_or_hiding_uphill_bed() {
    let mut m = line(&[0., 5., 1., 8.], &[0., 2., 3., 7.]);
    m.fixed[0] = Some((0., Terminal::WaterBody(17)));
    let r = route(&m).unwrap();
    assert_eq!(r.head, vec![0., 5., 5., 8.]);
    assert_eq!(m.bed, vec![0., 5., 1., 8.]);
    assert_eq!(r.receiver, vec![NONE, 0, 1, 2]);
    assert_eq!(r.audit.uphill_boundary_heads, 0);
    assert_eq!(r.audit.uphill_boundary_bed_links, 1);
    assert_eq!(r.audit.filled_nodes, 1);
    assert_eq!(r.discharge[0], 12.);
    assert!(r.owner.iter().all(|t| *t == Terminal::WaterBody(17)));
}
#[test]
fn two_terminals_have_independent_source_budgets_and_divide_can_cross_cell() {
    let mut m = line(&[0., 2., 7., 3., 1.], &[0.2, 1., 2., 4., 0.3]);
    m.fixed[0] = Some((0., Terminal::WaterBody(10)));
    m.fixed[4] = Some((1., Terminal::ClosedFloor(4)));
    let r = route(&m).unwrap();
    assert_eq!(r.audit.budgets.len(), 2);
    assert!(r.audit.budgets.iter().all(|b| b.relative_residual < 1e-14));
    assert_eq!(r.owner[1], Terminal::WaterBody(10));
    assert_eq!(r.owner[3], Terminal::ClosedFloor(4));
    assert!(r.audit.total_source_m3_s > 7.49);
    // Perturb the terrain, not an allowed receiver list: the central divide moves.
    m.bed[1] = 6.;
    let moved = route(&m).unwrap();
    assert_ne!(r.owner[2], moved.owner[2]);
}
#[test]
fn lakes_own_inflows_without_forcing_their_declared_outlet_to_ocean() {
    let flow = [1, 2, 3, 4, NONE, 2];
    let bodies = [NONE, NONE, 42, NONE, 7, NONE];
    let owners = terminal_owners(&flow, &bodies).unwrap();
    assert_eq!(owners[0], Terminal::WaterBody(42));
    assert_eq!(owners[5], Terminal::WaterBody(42));
    assert_eq!(owners[3], Terminal::WaterBody(7));
    assert!(terminal_owners(&[1, 0], &[NONE, NONE]).is_err());
}
#[test]
fn unsupported_components_and_heads_below_bed_fail_loudly() {
    let mut m = line(&[5., 6.], &[1., 1.]);
    assert!(route(&m).is_err());
    m.fixed[0] = Some((4., Terminal::WaterBody(0)));
    assert!(route(&m).is_err());
    m.fixed[0] = Some((5., Terminal::WaterBody(0)));
    m.edges[0] = vec![NONE; 8];
    m.edges[1] = vec![NONE; 8];
    assert!(route(&m).is_err());
}
#[test]
fn longitude_wrap_and_spherical_area_are_global() {
    let g = Grid::new(4, 8, 6_371_000.).unwrap();
    assert_eq!(g.neighbor(8, -1, 0), Some(15));
    assert_eq!(g.neighbor(0, 0, -1), None);
    let total = (0..g.len()).map(|i| g.area(i, 0., 1., 1.)).sum::<f64>();
    assert!(
        relative(
            total - 4. * std::f64::consts::PI * g.radius * g.radius,
            total
        ) < 1e-14
    );
}

#[test]
fn catchment_includes_upstream_lake_and_its_donors() {
    let flow = [1, 2, 3, 4, NONE];
    let bodies = [NONE, 8, NONE, 9, NONE];
    let owners = terminal_owners(&flow, &bodies).unwrap();
    let core = domain::upstream_closure(&flow, &owners, &Terminal::WaterBody(9)).unwrap();
    assert_eq!(core, vec![true, true, true, true, false]);
}

#[test]
fn parent_head_bounds_do_not_mutate_beds_and_account_actual_rerouting() {
    let mut m = line(&[0., 1., 0.5, 4.], &[0., 1., 2., 3.]);
    m.fixed[0] = Some((0., Terminal::WaterBody(0)));
    let old = m.bed.clone();
    let r = route_with_floor(&m, Some(&[0., 2., 0.5, 4.])).unwrap();
    assert_eq!(m.bed, old);
    assert_eq!(r.head, vec![0., 2., 2., 4.]);
    assert_eq!(r.discharge[0], 6.);
    assert_eq!(r.audit.uphill_heads, 0);
    assert!(r.audit.relative_residual < 1e-14);
    assert!(route_with_floor(&m, Some(&[1., 2., 0.5, 4.])).is_err());
}
