use ndarray::{prelude::s, Array, Array2, Array3, Axis, Zip};
use wasm_bindgen::prelude::*;

// macro_rules! log {
//     ( $( $t:tt )* ) => {
//         web_sys::console::log_1(&format!( $( $t )* ).into());
//     }
// }

// #[cfg(feature = "wee_alloc")]
// #[global_allocator]
// static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

const T: f64 = 298.15; // temperature in Kelvin
const K: f64 = 1.3806505E-23; // Boltzmann's constant in J/K
const KG_AMU: f64 = 1.660538921E-27; // kg/amu conversion factor
const PI: f64 = 3.141_592_653_589_793; // Pi
const E_CHARGE: f64 = 1.60217663E-19;
const REFINE_VOLTS: f64 = 100.0;
const SCALE: f64 = 1000.0;
const LARGE_NUMBER: f64 = 1844674407370955161.0;
const DT: f64 = 1E-9 * SCALE; // 1000 mm/m conversion. 1 ns step size
const EDGE_FACTOR: f64 = 3.0;
const CORNER_FACTOR: f64 = 4.5;

#[wasm_bindgen]
pub struct Environment {
    width: usize,
    height: usize,
    scale: usize,
    volts: Vec<f64>, // short vec of the volts of electrodes
    solids: Vec<i8>, // short vec of electrodes that are solid
    pixels: Vec<u8>,
    blank: Array2<i8>, // width * height of empty and magnetic
    emap: Array2<f64>,
    efmap: Array3<f64>,
    timeline: Array2<f64>,
    offsets: Vec<f64>,
    pressure: f64, // in Pa
    gas_mass: f64, // in amu / Da
    fdm_threshold: f64,
    magnets: Vec<f64>, // only 2 values possible - 0.0 is not magnet, 1.0 is magnet;
}

#[wasm_bindgen]
impl Environment {
    #[wasm_bindgen(constructor)]
    pub fn new(
        width: usize, height: usize, scale: usize, time: usize, volts: Vec<f64>, solids: Vec<i8>,
        pressure: f64, gas_mass: f64,
    ) -> Environment {
        let vlen = volts.len();
        let emap = Array2::<f64>::zeros((width, height));
        let blank = emap.mapv(|_| -1); // -1 = background; 0 = not splattable; 1 = splattable
        let mut efmap = Array3::zeros((width, height, vlen));
        for i in 0..(width) {
            for j in 0..height {
                efmap[[j, i, 1]] = (i + j) as f64;
            }
        }
        Environment {
            width,
            height,
            scale,
            pixels: vec![0; (width / scale) * (height / scale)],
            emap,
            solids,
            blank,
            efmap,
            volts,
            timeline: Array2::<f64>::ones((vlen, time)),
            offsets: vec![0.0; vlen],
            pressure,
            gas_mass,
            fdm_threshold: 10.0,
            magnets: vec![0.0, 0.0, 0.0, 0.0, 0.0],
        }
    }

    pub fn efpix(&self) -> Vec<f64> {
        (&self.efmap * Array::from_vec(self.volts.clone()))
            .sum_axis(Axis(2))
            .iter()
            .copied()
            .collect()
    }

    pub fn electrode_pixels(&self) -> Vec<u8> { self.pixels.clone() }

    pub fn timeline_length(&self) -> usize { self.timeline.ncols() }

    pub fn sum(&self) -> usize { self.pixels.iter().map(|&x| x as usize).sum() }

    pub fn width(&self) -> usize { self.width }

    pub fn height(&self) -> usize { self.height }

    pub fn scale(&self) -> usize { self.scale }

    pub fn update_voltages(&mut self, volts: Vec<f64>) { self.volts = volts; }

    pub fn update_solids(&mut self, solids: Vec<i8>) { self.solids = solids; }

    pub fn update_offsets(&mut self, offsets: Vec<f64>) { self.offsets = offsets; }

    pub fn update_magnets(&mut self, magnets: Vec<f64>) { self.magnets = magnets; }

    pub fn update_pressure(&mut self, pressure: f64) { self.pressure = pressure; }

    pub fn update_gas_mass(&mut self, gas_mass: f64) { self.gas_mass = gas_mass; }

    pub fn update_fdm_threshold(&mut self, threshold: f64) { self.fdm_threshold = threshold; }

    pub fn update_timeline(&mut self, dcs: Vec<f64>, fs: Vec<f64>, sts: Vec<f64>, eds: Vec<f64>) {
        for i in 0..self.timeline.nrows() {
            let frequency = fs[i] * SCALE * PI * 2.0 * 1E-9;
            let (dc_mode, pls_start, pls_end) = (dcs[i], sts[i] as usize, eds[i] as usize);
            if dc_mode > 0.0 {
                for j in 0..self.timeline.ncols() {
                    self.timeline[[i, j]] = (j < pls_end && j >= pls_start) as i32 as f64;
                }
            } else {
                for j in 0..self.timeline.ncols() {
                    self.timeline[[i, j]] = (frequency * j as f64).sin();
                }
            }
            // self.timeline = Array2::from_shape_vec((self.volts.len(), self.timeline.ncols()), new_time).unwrap()
        }
    }

    pub fn clear(&mut self) {
        self.pixels = vec![0; (self.width / self.scale) * (self.height / self.scale)];
        self.blank = self.blank.mapv(|_| -1);
    }

    pub fn getef(&self, x: usize, y: usize) -> f64 {
        self.volts.iter().enumerate().map(|(i, v)| self.efmap[[y, x, i]] * v).sum()
    }

    pub fn scalars(&self, time: usize) -> Vec<f64> {
        self.timeline
            .slice(s![.., time])
            .iter()
            .enumerate()
            .map(|(i, s)| (s * self.volts[i]) + self.offsets[i])
            .collect()
    }

    pub fn save_simion_pa(&self) -> Vec<u8> {
        let max: f64 = 100_000.0;
        let mut pa: Vec<u8> = vec![];
        pa.extend(&(-1i32).to_le_bytes());
        pa.extend(&(1i32).to_le_bytes());
        pa.extend(&(max).to_le_bytes());
        pa.extend(&(self.width as i32).to_le_bytes());
        pa.extend(&(self.height as i32).to_le_bytes());
        pa.extend(&(1i32).to_le_bytes());
        pa.extend(&(0i32).to_le_bytes());
        for y in (0..self.width).rev() {
            for x in 0..self.height {
                let v = self.emap[[y, x]] + if self.blank[[y, x]] == -1 { 0.0 } else { 2.0 * max };
                pa.extend(&v.to_le_bytes());
            }
        }
        pa
    }

    pub fn brush(&mut self, x: usize, y: usize, color: u8) {
        let index = (y * self.width / self.scale) + x;
        if index < self.pixels.len() {
            self.pixels[index] = color;
        }
    }

    pub fn generate_electrodes(&mut self) {
        let pix_row_width = self.width / self.scale;
        for (i, &value) in self.pixels.iter().enumerate() {
            let (x, y, v) = (i % pix_row_width * 2, i / pix_row_width * 2, value as usize);
            for &pixel in &[[y, x], [y, x + 1], [y + 1, x], [y + 1, x + 1]] {
                self.emap[pixel] = self.volts[v];
                self.blank[pixel] = if self.magnets[v] != 0.0 { -1 } else { self.solids[v] };
            }
        }
    }

    pub fn generate_electric_fields(&mut self) {
        self.generate_electrodes();
        self.efmap = Array3::zeros((self.width, self.height, self.volts.len()));
        for (i, &volts) in self.volts.iter().enumerate() {
            if i == 0 || !&self.emap.iter().any(|v| v == &volts) || self.magnets[i] != 0.0 {
                continue; // skip if voltage not found in the elec map or is blank/magnet
            }
            let mut e_field = Array2::<f64>::zeros((self.emap.ncols(), self.emap.nrows()));
            Zip::from(&mut e_field).and(&self.emap).and(&self.blank).for_each(|ef, &e, &c| {
                *ef = if (e == volts) && (c != -1) { REFINE_VOLTS } else { 0.0 };
            });
            let this = e_field.clone();
            let all = self.blank.mapv(|elec| if elec == 0 || elec == 1 { 1.0 } else { 0.0 });
            for &step in &[32, 16, 8, 4, 2, 1] {
                e_field = fdm_step(&e_field, &this, &all, step, REFINE_VOLTS, self.fdm_threshold);
                e_field = resize(&e_field, self.width, self.height);
            }
            let mut field = self.efmap.slice_mut(s![.., .., i]);
            field += &(e_field / REFINE_VOLTS);
        }
    }

    /// As passing values to and from rust/js is somewhat annoying, the "ion" is a single vec
    /// and the flight information is also a single vec.
    /// Ion fields: (1) x_position (2) y_position
    ///             (3) x_velocity (4) y_velocity
    ///             (5) mz         (6) rng seed
    /// The returned flight record will be time then x-position then y-position
    pub fn fly_ion(&self, ion: Vec<f64>) -> Vec<f64> {
        let (x, y) = (ion[0], ion[1]);
        let (xvel, yvel) = (ion[2] * SCALE, ion[3] * SCALE);
        let (mz, mass_kg) = (ion[4], ion[4] * KG_AMU); // use mz
        let (mut rng, ccs) = (Rng::new((ion[5] * LARGE_NUMBER) as u64), ccs(mz, self.gas_mass));
        let (mut pos, mut vel) = (Point { x, y }, Point { x: xvel, y: yvel });
        let (mut last_speed, mut mean_free_path) = (f64::MAX, f64::MAX);
        let (mut flight_record, mut collisions) = (vec![0.0, pos.x, pos.y], 0);

        for step in 1..self.timeline.ncols() {
            let e_vector = self.get_electric_vector(pos, step);
            let m_vector = self.get_magnetic_vector(pos, vel);
            let force = (e_vector + m_vector) * E_CHARGE;
            let acceleration = force / mass_kg;
            vel += acceleration * DT;
            pos += vel * DT;
            if self.pressure != 0.0 {
                let speed = vel.magnitude().max(1E-3);
                if (speed / last_speed - 1.0).abs() > 0.05 || step == 1 {
                    last_speed = speed;
                    mean_free_path = mfp(speed, ccs, self.pressure, self.gas_mass * KG_AMU);
                }
                if rng.f64() < (1.0 - (-speed * DT / mean_free_path).exp()) {
                    collisions += 1; // we hit something
                    vel = new_velocity(vel, mz, self.gas_mass, &mut rng);
                }
            }
            flight_record.extend_from_slice(&[DT * step as f64, pos.x, pos.y]);
            if pos.splatted((self.width - 1) as f64, (self.height - 1) as f64, &self.blank) {
                break;
            }
        }
        flight_record.push(collisions as f64);
        flight_record
    }

    fn get_electric_vector(&self, pos: Point, step: usize) -> Point {
        let scalars = self.scalars(step);
        let (xl, yl) = (pos.x.floor() as usize, pos.y.floor() as usize);
        let (xh, yh) = (pos.x.ceil().min(199.0) as usize, pos.y.ceil().min(199.0) as usize);
        let (xr, yr) = (pos.x.round().min(199.0) as usize, pos.y.round().min(199.0) as usize);
        let (mut vx, mut vy) = (0.0, 0.0);
        for (i, scalar) in scalars.iter().enumerate() {
            vy += (self.efmap[[yh, xr, i]] - self.efmap[[yl, xr, i]]) * scalar;
            vx += (self.efmap[[yr, xh, i]] - self.efmap[[yr, xl, i]]) * scalar;
        }
        -Point { x: vx, y: vy }
    }

    fn get_magnetic_vector(&self, pos: Point, vel: Point) -> Point {
        let pixel = [pos.y.floor() as usize, pos.x.floor() as usize];
        let mag = if self.blank[pixel] == -1 { self.emap[pixel] } else { 0.0 };
        Point { x: vel.y * mag, y: -1.0 * vel.x * mag } / SCALE
    }
}


fn fdm_step(
    e_field: &Array2<f64>, this: &Array2<f64>, all: &Array2<f64>, step: usize, volts: f64,
    threshold: f64,
) -> Array2<f64> {
    let mut e_now = e_field.to_owned().slice_mut(s![..;step, ..;step]).to_owned();
    let (width, height) = (e_field.ncols(), e_field.nrows());
    let (row, col) = (e_now.nrows(), e_now.ncols());
    let (mut hidden_this, mut hidden_all) = (Array2::zeros((row, col)), Array2::zeros((row, col)));
    for i in 0..row {
        for j in 0..col {
            let (x, y) = (j * step, i * step);
            let x_jump = if width % step == 0 || j < row - 1 { step } else { width % step };
            let y_jump = if height % step == 0 || i < col - 1 { step } else { height % step };
            if this.slice(s![y..y + y_jump, x..x + x_jump]).sum() != 0.0 {
                hidden_this[[i, j]] = volts;
            }
            if all.slice(s![y..y + y_jump, x..x + x_jump]).sum() == 0.0 {
                hidden_all[[i, j]] = 1.0; // so we have an elec somewhere inside.
            }
        }
    }
    let mut err = f64::INFINITY;
    while err > (threshold / (step as f64 * step as f64)) {
        let e_prev = e_now.clone();
        let mut fdm: Array2<f64> = Array2::zeros((row + 2, col + 2));
        let mut above = fdm.slice_mut(s![1..-1, 0..-2]);
        above += &e_now;
        let mut below = fdm.slice_mut(s![1..-1, 2..]);
        below += &e_now;
        let mut left = fdm.slice_mut(s![0..-2, 1..-1]);
        left += &e_now;
        let mut right = fdm.slice_mut(s![2.., 1..-1]);
        right += &e_now;
        // extract the center and update boundary conditions
        let mut boundaries = fdm.slice_mut(s![1..-1, 1..-1]);
        let mut center = boundaries.slice_mut(s![1..-1, 1..-1]);
        center /= 4.0;
        let mut first_row = boundaries.row_mut(0);
        first_row /= EDGE_FACTOR;
        let mut last_row = boundaries.row_mut(row - 1);
        last_row /= EDGE_FACTOR;
        let mut first_col = boundaries.column_mut(0);
        first_col /= EDGE_FACTOR;
        let mut last_col = boundaries.column_mut(col - 1);
        last_col /= EDGE_FACTOR;
        boundaries[[0, 0]] *= CORNER_FACTOR;
        boundaries[[0, row - 1]] *= CORNER_FACTOR;
        boundaries[[row - 1, 0]] *= CORNER_FACTOR;
        boundaries[[row - 1, row - 1]] *= CORNER_FACTOR;
        e_now = boundaries.to_owned();
        e_now *= &hidden_all;
        e_now += &hidden_this;
        err = (&e_now - &e_prev).mapv(|x| x.abs()).iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    }
    e_now
}


fn resize(e_field: &Array2<f64>, new_width: usize, new_height: usize) -> Array2<f64> {
    let (row, col) = (e_field.nrows() as f64 - 1.0, e_field.ncols() as f64 - 1.0);
    let (x_ratio, y_ratio) = (col / (new_width - 1) as f64, row / (new_height - 1) as f64);
    let mut resized = Array2::zeros((new_height, new_width));
    for i in 0..new_height {
        let y_l = (y_ratio * i as f64).floor() as usize;
        let y_h = (y_ratio * i as f64).ceil() as usize;
        let y_weight = (y_ratio * i as f64) - y_l as f64;
        for j in 0..new_width {
            let x_l = (x_ratio * j as f64).floor() as usize;
            let x_h = (x_ratio * j as f64).ceil() as usize;
            let x_weight = (x_ratio * j as f64) - x_l as f64;
            let a = e_field[[y_l, x_l]] * (1.0 - x_weight) * (1.0 - y_weight);
            let b = e_field[[y_l, x_h]] * x_weight * (1.0 - y_weight);
            let c = e_field[[y_h, x_l]] * y_weight * (1.0 - x_weight);
            let d = e_field[[y_h, x_h]] * x_weight * y_weight;
            resized[[i, j]] = a + b + c + d;
        }
    }
    resized
}

fn new_velocity(vel: Point, mass_amu: f64, gas_mass: f64, rng: &mut Rng) -> Point {
    let vel_stdev_gas = (K * T / (gas_mass * KG_AMU)).sqrt() / SCALE;
    let gas_vel = Point { x: rng.gaussian() * vel_stdev_gas, y: rng.gaussian() * vel_stdev_gas };
    // from https://en.wikipedia.org/wiki/Elastic_collision
    let v1 = (vel - gas_vel).magnitude(); // translate gas vel to zero & get magnitude
    let ion_angle = vel.y.atan2(vel.x);
    let (m1, m2, theta) = (mass_amu, gas_mass, 2.0 * PI * rng.f64() * 0.9999999999);
    let (costh, sinth, m1_2, m2_2) = (theta.cos(), theta.sin(), m1 * m1, m2 * m2);
    let delta_angle = (m2 * sinth).atan2(m1 + m2 * costh);
    let magnitude = v1 * (m1_2 + m2_2 + 2.0 * m1 * m2 * costh).sqrt() / (m1 + m2);
    // don't care about angle/magnitude for gas particle
    let angle = ion_angle + delta_angle;
    let vel_temp = Point { x: magnitude * angle.cos(), y: magnitude * angle.sin() };
    vel_temp + gas_vel // return frame of reference
}

/// speed in mm/us (m/ms), ccs in m^2, pressure in Pa, gas_mass in kg; λ = k * T / (√2 * π * d² * p)
fn mfp(speed: f64, ccs: f64, pressure: f64, gas_mass: f64) -> f64 {
    let c_bar = (8.0 * K * T / PI / gas_mass).sqrt() / SCALE;
    let s = speed / ((2.0 * K * T / gas_mass).sqrt() / SCALE);
    let c_bar_rel = c_bar * (s + 1.0 / (2.0 * s)) * 0.5 * PI.sqrt() * erf(s) + 0.5 * (-s * s).exp();
    K * T * (speed / c_bar_rel) / (pressure * ccs) * SCALE // mfp in mm
}

/// emperical formula that uses numbers from DOI: 10.3390/polym11040688
fn ccs(mass1: f64, mass2: f64) -> f64 {
    let ccs1 = 8.9 * (mass1 / 6.0).powf(2.0 / 3.0); // in angstroms squared
    let ccs2 = 8.9 * (mass2 / 6.0).powf(2.0 / 3.0); // in angstroms squared
    (((ccs1 / PI).sqrt() + (ccs2 / PI).sqrt()) * 1E-10).powi(2) * PI // ccs in meters squared
}

#[rustfmt::skip]
/// Error function (erf) stolen shamelessly from SIMION's collision_hs1.lua file
//  erf(z) = (2/sqrt(pi)) * integral[0..z] exp(-t^2) dt
fn erf(z: f64) -> f64 {
    let z2 = z.abs();
    let t = 1.0 / (1.0 + 0.32759109962 * z2);
    let mut res = (-1.061405429) * t;
    res = (res + 1.453152027) * t;
    res = (res - 1.421413741) * t;
    res = (res + 0.2844966736) * t;
    res = ((res - 0.254829592) * t) * (-z2 * z2).exp();
    res += 1.0;
    if z < 0.0 { -res } else { res }
}


#[derive(Copy, Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl core::ops::Add for Point {
    type Output = Self;

    fn add(self, rhs: Self) -> Self { Self { x: self.x + rhs.x, y: self.y + rhs.y } }
}

impl core::ops::AddAssign for Point {
    fn add_assign(&mut self, rhs: Self) { *self = Self { x: self.x + rhs.x, y: self.y + rhs.y }; }
}

impl core::ops::Sub for Point {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self { Self { x: self.x - rhs.x, y: self.y - rhs.y } }
}

impl core::ops::Mul for Point {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self { Self { x: self.x * rhs.x, y: self.y * rhs.y } }
}

impl core::ops::Neg for Point {
    type Output = Self;

    fn neg(self) -> Self { Self { x: -self.x, y: -self.y } }
}

impl core::ops::Div<f64> for Point {
    type Output = Self;

    fn div(self, rhs: f64) -> Self { Self { x: self.x / rhs, y: self.y / rhs } }
}

impl core::ops::Mul<f64> for Point {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self { Self { x: self.x * rhs, y: self.y * rhs } }
}

impl Point {
    fn splatted(self, x: f64, y: f64, blank: &ndarray::Array2<i8>) -> bool {
        self.x < 1.0 || self.y < 1.0 || self.x >= x || self.y >= y || blank[self.to_coords()] == 0
    }

    fn to_coords(self) -> [usize; 2] { [self.y as usize, self.x as usize] }

    fn magnitude(self) -> f64 { (self.x * self.x + self.y * self.y).sqrt() }
}


struct Rng(u64, u64);

impl Rng {
    const fn new(n: u64) -> Rng { Rng(n ^ 0xf4dbdf2183dcefb7, n ^ 0x1ad5be0d6dd28e9b) }

    fn f64(&mut self) -> f64 {
        let (mut x, y) = (self.0, self.1);
        x ^= x << 23;
        self.0 = y;
        self.1 = x ^ y ^ (x >> 17) ^ (y >> 26);
        (self.1.wrapping_add(y) >> 32) as f64 * 2.3283064365386963E-10
    }

    #[rustfmt::skip]
    fn gaussian(&mut self) -> f64 {
        let (v1, v2) = (2.0 * self.f64() - 1.0, 2.0 * self.f64() - 1.0);
        let s = v1 * v1 + v2 * v2;
        if s < 1.0 && s != 0.0 { v1 * (-2.0 * s.ln() / s).sqrt() } else { self.gaussian() }
    }
}
