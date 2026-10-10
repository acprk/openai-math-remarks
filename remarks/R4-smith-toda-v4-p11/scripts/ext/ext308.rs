//! Ext over the reduced-power algebra P (q_part = 0) at an odd prime p.
//! Usage: ext308 p s_max w_max  s1,w1 s2,w2 ...   (weight w = t/(2p-2))
use ext::chain_complex::{ChainComplex, FreeChainComplex};
use sseq::coordinates::Bidegree;

fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let p: i32 = a[1].parse()?;
    let s_max: i32 = a[2].parse()?;
    let w_max: i32 = a[3].parse()?;
    let q = 2 * p - 2;
    let json = serde_json::json!({
        "p": p, "algebra": ["milnor"],
        "profile": {"truncated": false, "q_part": 0, "p_part": []},
        "type": "finite dimensional module", "gens": {"x0": 0}, "actions": []
    });
    let save = std::env::var("EXT308_SAVE").ok().map(std::path::PathBuf::from);
    let res = ext::utils::construct((json, "milnor"), save)?;
    let t0 = std::time::Instant::now();
    res.compute_through_bidegree(Bidegree::s_t(s_max, q * w_max));
    eprintln!("resolved through s={s_max}, t={} in {:?}", q * w_max, t0.elapsed());
    // sanity: weight-1 and weight-p classes
    for (s, w) in [(1, 1), (1, p), (2, p), (2, p * p), (1, p + 1), (2, 2)] {
        println!("sanity Ext^({s},w={w}) = {}", res.number_of_gens_in_bidegree(Bidegree::s_t(s, q * w)));
    }
    for arg in &a[4..] {
        let v: Vec<i32> = arg.split(',').map(|x| x.parse().unwrap()).collect();
        let (s, w) = (v[0], v[1]);
        println!("Ext^(s={s}, w={w}) [t={}] dim = {}", q * w, res.number_of_gens_in_bidegree(Bidegree::s_t(s, q * w)));
    }
    Ok(())
}
