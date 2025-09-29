#[cfg(feature = "bind-threads")]
use hwlocality::{
    cpu::cpuset::CpuSet,
    object::types::ObjectType,
    topology::support::{DiscoverySupport, FeatureSupport},
    Topology,
};

#[cfg(not(feature = "bind-threads"))]
pub fn bind_threads() {}

#[cfg(feature = "bind-threads")]
pub fn bind_threads() {
    use std::sync::Arc;

    use hwlocality::{cpu::binding::CpuBindingFlags, Topology};
    use rayon::ThreadPoolBuilder;

    let builder = ThreadPoolBuilder::new();
    let topo = Arc::new(Topology::new().unwrap());
    let n_thread = std::env::var("RAYON_NUM_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);
    if let Some(cores) = get_proc_list(&topo) {
        let mut cores = cores.into_iter().cycle();
        builder
            .num_threads(n_thread)
            .spawn_handler(|t_builder| {
                let topo = topo.clone();
                let core = cores.next().expect("E: unreachable"); // due to cycle

                std::thread::spawn(move || {
                    // bind
                    let tid = hwlocality::current_thread_id();
                    topo.bind_thread_cpu(tid, &core, CpuBindingFlags::empty())
                        .unwrap();

                    // work
                    t_builder.run();
                });

                Ok(())
            })
            .build_global()
            .unwrap();
    } else {
        builder
            .num_threads(n_thread)
            .spawn_handler(|t_builder| {
                std::thread::spawn(move || {
                    // work
                    t_builder.run();
                });

                Ok(())
            })
            .build_global()
            .unwrap();
    };
}

#[cfg(feature = "bind-threads")]
#[derive(Debug, Default)]
enum BindingPolicy {
    /// Disable thread binding.
    /// Corresponding `RAYON_PROC_BIND_VALUE`: `false`.
    Disable,
    /// Enable thread binding & prioritize binding of PUs over cores?.
    /// Corresponding `RAYON_PROC_BIND_VALUE`: `close`.
    Close,
    /// Enable thread binding & prioritize binding across cores over filling PUs?.
    /// Corresponding `RAYON_PROC_BIND_VALUE`: `spread`. Default value.
    #[default]
    Spread,
}

/// Environment variable controlling the thread-binding policy.
///
/// The name of this variable and its possible values reflect the OpenMP equivalents.
#[cfg(feature = "bind-threads")]
const RAYON_PROC_BIND_VAR: &str = "RAYON_PROC_BIND";

#[cfg(feature = "bind-threads")]
impl BindingPolicy {
    fn from_env() -> Self {
        match std::env::var(RAYON_PROC_BIND_VAR) {
            Ok(val) => match val.to_lowercase().as_str() {
                "false" => Self::Disable,
                "close" => Self::Close,
                "spread" => Self::Spread,
                "" => Self::default(),
                _ => {
                    eprintln!("W: unrecognized RAYON_PROC_BIND value (!= false | close | spread)");
                    eprintln!("   continuing with default (spread)");
                    Self::default()
                }
            },
            Err(e) => {
                match e {
                    std::env::VarError::NotPresent => {}
                    std::env::VarError::NotUnicode(_) => {
                        eprintln!("W: non-unicode RAYON_PROC_BIND value");
                        eprintln!("   continuing with default (spread)");
                    }
                }
                Self::default()
            }
        }
    }
}

#[cfg(feature = "bind-threads")]
fn check_hwloc_support(topology: &Topology) -> Result<(), String> {
    if !topology.supports(FeatureSupport::discovery, DiscoverySupport::pu_count) {
        return Err("missing PU reporting support".to_string());
    }
    if !topology
        .feature_support()
        .cpu_binding()
        .is_some_and(|s| s.get_thread() && s.set_thread())
    {
        return Err("missing binding support".to_string());
    }

    Ok(())
}

/// Return a list of bind targets ordered according to the desired policy.
///
/// The desired policy is read from an environment variable (see [`RAYON_PROC_BIND_VAR`]). For details on each policy,
/// see [`BindingPolicy`].
///
/// The returned list is used by iterating over a `cycle`d version of it, which corresponds to a round robin logic.
#[cfg(feature = "bind-threads")]
fn get_proc_list(topology: &Topology) -> Option<Vec<CpuSet>> {
    if check_hwloc_support(topology).is_err() {
        return None;
    }

    let binding_policy = BindingPolicy::from_env();
    let core_depth = topology
        .depth_or_below_for_type(ObjectType::Core)
        .expect("E: unreachable");

    match binding_policy {
        BindingPolicy::Disable => None,
        BindingPolicy::Close => {
            let mut pu_set = Vec::with_capacity(256);
            topology.objects_at_depth(core_depth).for_each(|c| {
                let target = c.cpuset().unwrap().clone_target();
                let w = target.weight();
                if !(w == Some(1) || w == Some(2)) {
                    panic!()
                }
                target
                    .iter_set()
                    .map(CpuSet::from)
                    .for_each(|t| pu_set.push(t));
            });
            Some(pu_set)
        }
        BindingPolicy::Spread => {
            let mut first_pu_set = Vec::with_capacity(128);
            let mut second_pu_set = Vec::with_capacity(128);
            topology.objects_at_depth(core_depth).for_each(|c| {
                let target = c.cpuset().expect("E: unreachable").clone_target();
                // match required bc some modern CPUs have HT/SMT on only some of their cores :)
                match target.weight() {
                    Some(1) => {
                        // one PU per core -> no HT/SMT
                        first_pu_set.push(target);
                    }
                    Some(2) => {
                        // two PUs per core -> HT/SMT
                        let [first_pu, second_pu]: [CpuSet; 2] = target
                            .iter_set()
                            .map(CpuSet::from)
                            .collect::<Vec<_>>()
                            .try_into()
                            .expect("E: unreachable");
                        first_pu_set.push(first_pu);
                        second_pu_set.push(second_pu);
                    }
                    Some(_) | None => {
                        panic!("E: architecture too cursed")
                    }
                }
            });
            first_pu_set.append(&mut second_pu_set);
            Some(first_pu_set)
        }
    }
}
