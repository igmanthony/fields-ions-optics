use ndarray::{prelude::s, Array, Array2, Array3, Axis, Zip};
use wasm_bindgen::prelude::*;

macro_rules! log { ( $( $t:tt )* ) => { web_sys::console::log_1( &format!( $($t)* ).into()); } }

#[derive(Debug, Clone)]
struct SimulationError;
impl std::fmt::Display for SimulationError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result { write!(f, "simulation error") }
}
type R<T> = std::result::Result<T, SimulationError>;
type A2 = Array2<f64>;

const T: f64 = 298.15; // temperature in Kelvin
const K: f64 = 1.3806505E-23; // Boltzmann's constant in J/K
const KG_AMU: f64 = 1.660538921E-27; // kg/amu conversion factor
const PI: f64 = 3.141_592_653_589_793; // Pi
const E_CHARGE: f64 = 1.60217663E-19;
const REFINE_VOLTS: f64 = 10_000.0;
const SCALE: f64 = 1000.0; // 1000 mm/m conversion.
const LARGE_NUMBER: f64 = 1844674407370955161.0;
const DT: f64 = 1E-9; // 1 ns step size
const SCALAR_DISTANCE: f64 = 10_000.0 * SCALE; // I have no idea why this is so large
const EDGE_FACTOR: f64 = 3.0;
const CORNER_FACTOR: f64 = 4.5;


#[wasm_bindgen]
pub struct Environment {
    width: usize,
    height: usize,
    scale: usize,
    max_time: usize,
    pixels: Vec<u8>,
    solids: Vec<i8>,    // short vec of electrodes that are solid
    blank: Array2<i8>,  // width * height of empty and magnetic
    volts: Vec<f64>,    // short vec of the volts of electrodes
    emap: A2,           // 200 x 200
    efmap: Array3<f64>, // 200 x 200 x 5
    dc_mode: Vec<bool>,
    frequencies: Vec<f64>,
    pulse_starts: Vec<f64>,
    pulse_ends: Vec<f64>,
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
        width: usize, height: usize, scale: usize, max_time: usize, volts: Vec<f64>,
        solids: Vec<i8>, pressure: f64, gas_mass: f64,
    ) -> Environment {
        let (vlen, t) = (volts.len(), max_time as f64);
        let mut efmap = Array3::zeros((width, height, vlen));
        (0..width).for_each(|i| (0..height).for_each(|j| efmap[[j, i, 1]] = (i + j) as f64));
        Environment {
            width,
            height,
            scale,
            max_time,
            pixels: vec![0; (width / scale) * (height / scale)],
            solids,
            blank: Array2::<i8>::zeros((width, height)) - 1, // -1 = bg/mag; 0 = ghost; 1 = splat
            volts,
            emap: Array2::<f64>::zeros((width, height)),
            efmap,
            dc_mode: vec![true; vlen],
            frequencies: vec![0.0; vlen],
            pulse_starts: vec![0.0; vlen],
            pulse_ends: vec![0.0, t, t, t, t],
            offsets: vec![0.0; vlen],
            pressure,
            gas_mass,
            fdm_threshold: 10.0,
            magnets: vec![0.0; vlen],
        }
    }

    pub fn efpix(&self) -> Vec<f64> {
        (&self.efmap * Array::from_vec(self.volts.clone())).sum_axis(Axis(2)).into_iter().collect()
    }

    pub fn electrode_pixels(&self) -> Vec<u8> { self.pixels.clone() }

    pub fn timeline_length(&self) -> usize { self.max_time }

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

    pub fn update_max_time(&mut self, max_time_ns: usize) { self.max_time = max_time_ns; }

    pub fn update_timeline(&mut self, dcs: Vec<f64>, fs: Vec<f64>, sts: Vec<f64>, eds: Vec<f64>) {
        self.dc_mode = dcs.iter().map(|dc| dc > &0.0).collect();
        self.frequencies = fs.iter().map(|f| f * 1000.0).collect(); // convert freqs from kHz -> Hz
        self.pulse_starts = sts;
        self.pulse_ends = eds;
    }

    pub fn clear(&mut self) {
        self.pixels.iter_mut().for_each(|p| *p = 0);
        self.blank = self.blank.mapv(|_| -1);
    }

    pub fn getef(&self, x: usize, y: usize) -> f64 {
        self.volts.iter().enumerate().map(|(i, v)| self.efmap[[y, x, i]] * v).sum()
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
    /// The returned flight record will be time then x-position then y-position; finally collisions
    pub fn fly_ion(&self, ion: Vec<f64>) -> Vec<f64> {
        let (mut pos, mut vel) = (Point { x: ion[0], y: ion[1] }, Point { x: ion[2], y: ion[3] });
        let (mz, mass_kg) = (ion[4], ion[4] * KG_AMU); // use mz
        let (mut rng, ccs) = (Rng::new((ion[5] * LARGE_NUMBER) as u64), ccs(mz, self.gas_mass));
        let (mut last_speed, mut mean_free_path) = (f64::MAX, f64::MAX);
        let (mut flight_record, mut collisions) = (vec![0.0, pos.x, pos.y], 0);
        let (mut t, mut t_scalar) = (0.0, 10.0);
        log!("{}", self.max_time);
        while t < (self.max_time as f64 * DT) {
            let dt = DT * t_scalar;
            let (dp, dv) = match self.runge_kutta(pos, vel, t, dt, mass_kg) {
                Ok(value) => value,
                Err(_) => break,
            };
            pos += dp * SCALE;
            vel += dv;
            t += dt;
            let speed = vel.magnitude();
            if self.pressure != 0.0 {
                if (speed / last_speed - 1.0).abs() > 0.05 || t == 0.0 {
                    last_speed = speed;
                    mean_free_path = mfp(speed, ccs, self.pressure, self.gas_mass * KG_AMU);
                }
                if rng.f64() < (1.0 - (-speed * dt / mean_free_path).exp()) {
                    collisions += 1; // we hit something
                    vel = collide(vel, mz, self.gas_mass, &mut rng);
                }
            }
            flight_record.extend_from_slice(&[t, pos.x, pos.y]);
            t_scalar = match get_time_scalar(speed, dv.magnitude()) {
                Some(s) => s.min(100.0), // don't let scalar get too big
                None => break,
            };
            if pos.splatted((self.width - 3) as f64, (self.height - 3) as f64, &self.blank) {
                break;
            }
        }
        flight_record.push(collisions as f64);
        flight_record
    }

    
    /// fourth order Runge-Kutta method, adapted from
    /// https://gafferongames.com/post/integration_basics/
    fn runge_kutta(&self, pos: Point, vel: Point, t: f64, dt: f64, mass: f64) -> R<(Point, Point)> {
        let (rc1, rc2, rc3, rc4) = (0.0 * dt, 0.5 * dt, 0.5 * dt, 1.0 * dt);
        let (dp0, dv0) = (Point { x: 0.0, y: 0.0 }, Point { x: 0.0, y: 0.0 });
        let dp1 = vel + dv0 * rc1;
        let dv1 = self.acceleration(pos + dp0 * rc1, dp1, t + rc1, mass)?;
        let dp2 = vel + dv1 * rc2;
        let dv2 = self.acceleration(pos + dp1 * rc2, dp2, t + rc2, mass)?;
        let dp3 = vel + dv2 * rc3;
        let dv3 = self.acceleration(pos + dp2 * rc3, dp3, t + rc3, mass)?;
        let dp4 = vel + dv3 * rc4;
        let dv4 = self.acceleration(pos + dp3 * rc4, dp4, t + rc4, mass)?;
        let dp = (dp1 + (dp2 + dp3) * 2.0 + dp4) / 6.0 * dt;
        let dv = (dv1 + (dv2 + dv3) * 2.0 + dv4) / 6.0 * dt;
        Ok((dp, dv))
    }
    
    fn acceleration(&self, pos: Point, vel: Point, time: f64, mass: f64) -> R<Point> {
        if pos.x < 1.0 || pos.y < 1.0 || pos.x > 198.0 || pos.y > 198.0 {
            Err(SimulationError)
        } else {
            let e_vector = self.get_electric_vector(pos, time);
            let m_vector = self.get_magnetic_vector(pos, vel);
            let force = (e_vector + m_vector) * E_CHARGE;
            Ok(force / mass)
        }
    }

    /// four point spider (0.5 points in all direction):
    /// *     *---+-*     * || 00    01    02    03
    ///       |   | |       ||
    ///       |---O-|       ||
    /// *     *---+-*     * || 10    11    12    13
    ///        O--X--O      ||             
    ///           |         ||
    /// *     *   O *     * || 20    21    22    23
    ///                     ||
    ///                     ||
    /// *     *     *     * || 30    31    32    33
    /// Point x and y are between 1 and 2, so get fractional components
    fn get_electric_vector(&self, pos: Point, time: f64) -> Point {
        let (xf, x) = (pos.x as usize, pos.x.fract());
        let (yf, y) = (pos.y as usize, pos.y.fract());
        let xi = if x < 0.5 { xf - 1 } else { xf + 2 };
        let yi = if y < 0.5 { yf - 1 } else { yf + 2 };
        let (mut p0, mut p1, mut p2, mut p3) = (0.0, 0.0, 0.0, 0.0);
        let (mut p4, mut p5, mut p6, mut p7) = (0.0, 0.0, 0.0, 0.0);
        for i in 1..self.dc_mode.len() {
            let scalar = self.offsets[i] + self.volts[i] * self.timeline_scalar(i, time);
            p0 += self.efmap[[yf + 0, xf + 0, i]] * scalar;
            p1 += self.efmap[[yf + 0, xf + 1, i]] * scalar;
            p2 += self.efmap[[yf + 1, xf + 0, i]] * scalar;
            p3 += self.efmap[[yf + 1, xf + 1, i]] * scalar;
            p4 += self.efmap[[yf + 0, xi + 0, i]] * scalar;
            p5 += self.efmap[[yf + 1, xi + 0, i]] * scalar;
            p6 += self.efmap[[yi + 0, xf + 0, i]] * scalar;
            p7 += self.efmap[[yi + 0, xf + 1, i]] * scalar;
        }
        let ev_x = match x < 0.5 {
            true => lerp(x + 0.5, y, p4, p0, p5, p2) - lerp(x + 0.5, y, p0, p1, p2, p3),
            false => lerp(x - 0.5, y, p0, p1, p2, p3) - lerp(x - 0.5, y, p1, p4, p3, p5),
        };
        let ev_y = match y < 0.5 {
            true => lerp(x, y + 0.5, p6, p7, p0, p1) - lerp(x, y + 0.5, p0, p1, p2, p3),
            false => lerp(x, y - 0.5, p0, p1, p2, p3) - lerp(x, y - 0.5, p2, p3, p6, p7),
        };
        Point { x: ev_x, y: ev_y } * SCALE
    }

    fn get_magnetic_vector(&self, pos: Point, vel: Point) -> Point {
        let pixel = [pos.y as usize, pos.x as usize];
        let mag = if self.blank[pixel] == -1 { self.emap[pixel] } else { 0.0 };
        Point { x: vel.y * mag, y: -1.0 * vel.x * mag } / SCALE // div by scale due to flux
    }

    fn timeline_scalar(&self, i: usize, time: f64) -> f64 {
        match self.dc_mode[i] {
            true => (time < self.pulse_ends[i] && time > self.pulse_starts[i]) as i32 as f64,
            false => (self.frequencies[i] * PI * 2.0 * time).sin(),
        }
    }
}


fn get_time_scalar(speed: f64, accel: f64) -> Option<f64> {
    match (accel == 0.0, speed == 0.0) {
        (true, true) => None,
        (true, false) => Some(SCALAR_DISTANCE / speed),
        (false, true) => Some((2.0 * SCALAR_DISTANCE / accel).sqrt()),
        (false, false) => {
            let t_v = SCALAR_DISTANCE / speed;
            let t_a = (2.0 * SCALAR_DISTANCE / accel).sqrt();
            Some(t_v * t_a / (t_v + t_a))
        }
    }
}

fn lerp(x: f64, y: f64, q12: f64, q22: f64, q11: f64, q21: f64) -> f64 {
    (1.0 - y) * ((1.0 - x) * q11 + x * q21) + (y) * ((1.0 - x) * q12 + x * q22)
}

fn fdm_step(field: &A2, this: &A2, all: &A2, step: usize, volts: f64, threshold: f64) -> A2 {
    let mut e_now = field.to_owned().slice_mut(s![..;step, ..;step]).to_owned();
    let (width, height) = (field.ncols(), field.nrows());
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
    let mut fdm = Array2::zeros((row + 2, col + 2));
    while err > (threshold / (step as f64 * step as f64)) {
        fdm.fill(0.0);
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
        boundaries[[0, col - 1]] *= CORNER_FACTOR;
        boundaries[[row - 1, 0]] *= CORNER_FACTOR;
        boundaries[[row - 1, col - 1]] *= CORNER_FACTOR;
        let mut e_new = boundaries.to_owned();
        e_new *= &hidden_all;
        e_new += &hidden_this;
        err = (&e_new - &e_now).iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b.abs()));
        e_now = e_new;
    }
    e_now
}


/// linear interpolation resizing function for array -> used to double the size
fn resize(e_field: &A2, new_width: usize, new_height: usize) -> A2 {
    let (row, col) = (e_field.nrows() as f64 - 1.0, e_field.ncols() as f64 - 1.0);
    let (x_ratio, y_ratio) = (col / (new_width - 1) as f64, row / (new_height - 1) as f64);
    let mut resized = Array2::zeros((new_height, new_width));
    for i in 0..new_height {
        let y_l = (y_ratio * i as f64) as usize; // should floor
        let y_h = (y_ratio * i as f64).ceil() as usize;
        let y_weight = (y_ratio * i as f64) - y_l as f64;
        for j in 0..new_width {
            let x_l = (x_ratio * j as f64) as usize; // should floor
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

/// simulate collision and new trajectory of ion and background gas molecule
fn collide(vel: Point, mass_amu: f64, gas_mass: f64, rng: &mut Rng) -> Point {
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
    let c_bar = (8.0 * K * T / PI / gas_mass).sqrt();
    let s = speed / ((2.0 * K * T / gas_mass).sqrt());
    let c_bar_rel = c_bar * (s + 1.0 / (2.0 * s)) * 0.5 * PI.sqrt() * erf(s) + 0.5 * (-s * s).exp();
    K * T * (speed / c_bar_rel) / (pressure * ccs) // mfp in mm
}

/// emperical formula that uses numbers from DOI: 10.3390/polym11040688
fn ccs(mass1: f64, mass2: f64) -> f64 {
    let ccs1 = 8.9 * (mass1 / 6.0).powf(2.0 / 3.0); // in angstroms squared
    let ccs2 = 8.9 * (mass2 / 6.0).powf(2.0 / 3.0); // in angstroms squared
    (((ccs1 / PI).sqrt() + (ccs2 / PI).sqrt()) * 1E-10).powi(2) * PI // ccs in meters squared
}

#[rustfmt::skip]
/// Error function (erf) stolen shamelessly from SIMION's collision_hs1.lua file
///  erf(z) = (2/sqrt(pi)) * integral[0..z] exp(-t^2) dt
fn erf(z: f64) -> f64 {
    let z2 = z.abs();
    let t = 1.0 / (1.0 + 0.32759109962 * z2);
    let mut res = (-1.061405429) * t;
    res = (res + 1.453152027) * t;
    res = (res - 1.421413741) * t;
    res = (res + 0.2844966736) * t;
    res = ((res - 0.254829592) * t) * (-z2 * z2).exp() + 1.0;
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
        self.x < 2.0 || self.y < 2.0 || self.x >= x || self.y >= y || blank[self.to_coords()] == 0
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
