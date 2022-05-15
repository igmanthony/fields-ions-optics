use core::ops::{Add, AddAssign, Sub, Div, Mul, Neg};
use ndarray::{prelude::*, Array2, Array3, Axis, Zip};
use wasm_bindgen::prelude::*;


macro_rules! log {
    ( $( $t:tt )* ) => {
        web_sys::console::log_1(&format!( $( $t )* ).into());
    }
}

#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

const T: f64 = 298.15; // temperature in Kelvin
const K: f64 = 1.3806505E-23; // Boltzmann's constant in J/K
const KG_AMU: f64 = 1.660538921E-27; // kg/amu conversion factor
const EV_J: f64 = 6.2415095E+18; // eV/J conversion factor
const PI: f64 = 3.1415926535897931; // Pi

#[wasm_bindgen]
pub struct Environment {
    width: usize,
    height: usize,
    scale: usize,
    volts: Vec<f64>,      // short vec of the volts of electrodes
    splattables: Vec<i8>, // short vec of electrodes that are collidable
    pixels: Vec<u8>,
    collidables: Array2<i8>,
    emap: Array2<f64>,
    efmap: Array3<f64>,
    // for rf    
    timeline: Array2<f64>,
    offsets: Vec<f64>,
    // for gas collisions
    pressure: f64, // in Pa
    gas_mass: f64, // in amu / Da
    ccs: f64,
}

#[wasm_bindgen]
impl Environment {
    #[wasm_bindgen(constructor)]
    pub fn new(
        width: usize, height: usize, scale: usize, max_time: usize, volts: Vec<f64>,
        splattables: Vec<i8>, pressure: f64, gas_mass: f64
    ) -> Environment {
        let pixels = vec![0; (width / scale) * (height / scale)];
        let emap = Array2::<f64>::zeros((width, height));
        let collidables = emap.mapv(|_| -1); // -1 = background; 0 = not splattable; 1 = splattable
        let mut efmap = Array3::zeros((width, height, volts.len()));
        for i in 0..(width) {
            for j in 0..height {
                efmap[[j, i, 1]] = (i * j) as f64;
            }
        }
        let timeline = Array2::<f64>::ones((volts.len(), max_time));
        let offsets = vec![0.0; volts.len()];
        let ccs = 0.0;
        Environment {
            width,
            height,
            scale,
            pixels,
            emap,
            splattables,
            collidables,
            efmap,
            volts,
            timeline,
            offsets,
            pressure,
            gas_mass,
            ccs,
        }
    }

    pub fn electric_field_pixels(&self) -> Vec<u8> {
        let efsum = (&self.efmap * Array::from_vec(self.volts.clone())).sum_axis(Axis(2));
        let max = efsum.fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let min = efsum.fold(f64::INFINITY, |a, &b| a.min(b));
        efsum.iter().map(|&x| ((x - min) / (max - min) * 255.0) as u8).collect()
    }

    pub fn efpix(&self) -> Vec<f64> {
        (&self.efmap * Array::from_vec(self.volts.clone()))
            .sum_axis(Axis(2))
            .iter()
            .map(|&x| x)
            .collect()
    }

    /// Return the elec pixels for plotting
    pub fn electrode_pixels(&self) -> Vec<u8> { self.pixels.clone() }

    pub fn update_voltages(&mut self, volts: Vec<f64>) { self.volts = volts; }

    pub fn update_splattables(&mut self, splattables: Vec<i8>) { self.splattables = splattables;}

    pub fn update_timeline(&mut self, timeline: Vec<f64>) {
        self.timeline = Array2::from_shape_vec((self.volts.len(), self.timeline.ncols()), timeline)
            .expect("oops!");
    }

    pub fn update_offsets(&mut self, offsets: Vec<f64>) { self.offsets = offsets; }

    pub fn update_pressure(&mut self, pressure: f64) { self.pressure = pressure; }

    pub fn update_gas_mass(&mut self, gas_mass: f64) { self.gas_mass = gas_mass; }

    pub fn timeline_length(&self) -> usize { return self.timeline.ncols(); }

    pub fn sum(&self) -> usize { self.pixels.iter().map(|&x| x as usize).sum() }

    pub fn width(&self) -> usize { self.width }

    pub fn height(&self) -> usize { self.height }

    pub fn scale(&self) -> usize { self.scale }

    pub fn clear(&mut self) {
        self.pixels = vec![0; (self.width / self.scale) * (self.height / self.scale)];
        self.collidables = self.collidables.mapv(|_| -1);
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
        let max_voltage: f64 = 100_000.0;
        let mut pa: Vec<u8> = vec![];
        pa.extend(&(-1i32).to_le_bytes());
        pa.extend(&(1i32).to_le_bytes());
        pa.extend(&(max_voltage).to_le_bytes());
        pa.extend(&(self.width as i32).to_le_bytes());
        pa.extend(&(self.height as i32).to_le_bytes());
        pa.extend(&(1i32).to_le_bytes());
        pa.extend(&(0i32).to_le_bytes());
        for y in (0..self.width).rev() {
            for x in 0..self.height {
                let mut value = self.emap[[y, x]];
                value += if self.collidables[[y, x]] == 1 || self.collidables[[y, x]] == 0 {
                    2.0 * max_voltage
                } else {
                    0.0
                };
                pa.extend(&value.to_le_bytes());
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
                self.collidables[pixel] = self.splattables[v];
            }
        }
    }

    pub fn generate_electric_fields(&mut self) {
        self.generate_electrodes();
        self.efmap = Array3::zeros((self.width, self.height, self.volts.len()));
        for (i, &volts) in self.volts.iter().enumerate() {
            if (i == 0) || !&self.emap.iter().any(|v| v == &volts) {
                log!("Rust skipping #{} volts: {}", i, volts);
                continue; // continue if voltage not found in the elec map
            }
            log!("Rust refining #{} volts: {}", i, volts);
            let refine_volts = 100.0; // TODO: Replace this with 10_000.0
            // let mut e_field = self.emap.mapv(|x| if x == volts { refine_volts } else { 0.0 });
            let mut e_field = Array2::<f64>::zeros((self.emap.ncols(), self.emap.nrows()));
            Zip::from(&mut e_field).and(&self.emap).and(&self.collidables).for_each(
                |ef, &e, &c| {
                    *ef = if (e == volts) && (c != -1) { refine_volts } else { 0.0 };
                },
            );

            let this = e_field.clone();
            let all = self.collidables.mapv(|elec| if elec == 0 || elec == 1 { 1.0 } else { 0.0 });
            for &step in &[32, 16, 8, 4, 2, 1] {
                // "volts" is now 10_0000 volts for refining.
                e_field = fdm_step(&e_field, &this, &all, step, refine_volts);
                e_field = resize(&e_field, self.width, self.height);
            }
            let mut field = self.efmap.slice_mut(s![.., .., i]);
            e_field /= refine_volts;
            log!("e_field sum: {}", &e_field.sum());
            field += &e_field;
        }
    }

    /// As passing values to and from rust/js is somewhat annoying, the "ion" is a single vec
    /// and the flight information is also a single vec.
    /// Ion fields: (1) x_position (2) y_position
    ///             (3) x_velocity (4) y_velocity
    ///             (5) mz         (6) rng seed
    /// The returned flight record will be time then x-position then y-position
    pub fn fly_ion(&self, ion: Vec<f64>) -> Vec<f64> {
        if self.efmap.sum() == 0.0 {
            return vec![0.0, 0.0, 0.0]; // don't fly ions if there isn't an electric field
        }
        let max_time = self.timeline.ncols();
        let mut flight_record = vec![];
        let confusing_adjustment_factor = 1000.0; // I have no idea why this is needed
        let (x, y) = (ion[0], ion[1]);
        let (xvel, yvel) = (ion[2] * confusing_adjustment_factor, ion[3] * confusing_adjustment_factor);
        let mz = ion[4]; // use mz as mass in amu (essentially assume charge = 1)
        let mass_kg= mz * KG_AMU;
        let lims = Point { x: (self.width - 1) as f64, y: (self.height - 1) as f64 };
        let dt = 1e-9 * 1000.0; // 1000 mm/m conversion. 1 ns step size

        let mut time = 0.0;
        let mut pos = Point { x, y };
        let mut vel = Point { x: xvel, y: yvel };
        // for collisions
        let mut collisions = 0;
        let mut last_speed = vel.srss();
        let mut mean_free_path = f64::MAX;
        let ccs = if self.ccs <= 0.0 { calculate_ccs(mz, self.gas_mass) } else { self.ccs };
        let mut rng = Rng::new( (ion[5] * 1844674407370955161.0) as u64); // seed
        log!("velocity {}x {}y", xvel, yvel);
        flight_record.extend_from_slice(&[time, pos.x, pos.y]);
        let mut splat = false;
        // log!("Scalars: {:?}", self.scalars(0));
        for step in 1..max_time {
            let electric_vector = self.get_electric_vector(pos, step);
            let force = electric_vector * 1.602E-19;
            let acceleration = force / mass_kg;
            vel += acceleration * dt;
            pos += vel * dt;
            time = dt * step as f64;
            flight_record.extend_from_slice(&[time, pos.x, pos.y]);
            if self.pressure != 0.0 {
                let speed = vel.srss().max(1E-3); // calculate speed and 
                mean_free_path = if (speed / last_speed - 1.0).abs() > 0.05 || step == 1 { // changed
                    last_speed = speed; // update last speed to speed
                    calculate_mean_free_path(speed, ccs, self.pressure, self.gas_mass) // in mm
                } else {
                    mean_free_path
                };
                // NOTE: dt might need to be 1 ns rather than 1000x 1 ns...
                if rng.f64() < (1.0 - (-speed * dt / mean_free_path).exp()) { // test for collision
                    collisions += 1;
                    vel = calculate_new_velocity(vel, mz, self.gas_mass, &mut rng);
                }
            }
            if pos.splatted(lims, &self.collidables) {
                log!("SPLAT! {} us; x:{} y:{} s:{} collisions: {}", dt * step as f64, pos.x, pos.y, step, collisions);
                splat = true;
                break;
            }
        }
        if !splat {
            log!("Maximum time reached: {} us {} steps {} collisions", dt * max_time as f64, max_time, collisions);
        }
        flight_record
    }

    fn get_electric_vector(&self, pos: Point, step: usize) -> Point {
        let scalars = self.scalars(step);
        // log!("scalars {:?}", scalars);
        let (xl, yl) = (pos.x.floor() as usize, pos.y.floor() as usize);
        let (xh, yh) = (pos.x.ceil().min(199.0) as usize, pos.y.ceil().min(199.0) as usize);
        let (xr, yr) = (pos.x.round().min(199.0) as usize, pos.y.round().min(199.0) as usize);
        let mut vx = 0.0;
        let mut vy = 0.0;
        for (i, scalar) in scalars.iter().enumerate() {
            vy += (&self.efmap[[yh, xr, i]] - &self.efmap[[yl, xr, i]]) * scalar;
            vx += (&self.efmap[[yr, xh, i]] - &self.efmap[[yr, xl, i]]) * scalar;
        }
        -Point { x: vx, y: vy }
    }
}


fn fdm_step(
    e_field: &Array2<f64>, this: &Array2<f64>, all: &Array2<f64>, step: usize, volts: f64,
) -> Array2<f64> {
    let (width, height) = (e_field.ncols(), e_field.nrows());
    let (edge_divisor, corner_compensator) = (3.0, 4.5); // (3.04)
    let mut e_now = e_field.to_owned().slice_mut(s![..;step, ..;step]).to_owned();
    let (row, col) = (e_now.nrows(), e_now.ncols());
    let mut hidden_this = Array2::zeros((row, col));
    let mut hidden_all = Array2::zeros((row, col));
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
    while err > (10.0 / (step as f64 * step as f64)) {
        let e_prev = e_now.clone();
        // perform the finite difference method
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
        first_row /= edge_divisor;
        let mut last_row = boundaries.row_mut(row - 1);
        last_row /= edge_divisor;
        let mut first_col = boundaries.column_mut(0);
        first_col /= edge_divisor;
        let mut last_col = boundaries.column_mut(col - 1);
        last_col /= edge_divisor;
        boundaries[[0, 0]] *= corner_compensator;
        boundaries[[0, row - 1]] *= corner_compensator;
        boundaries[[row - 1, 0]] *= corner_compensator;
        boundaries[[row - 1, row - 1]] *= corner_compensator;
        e_now = boundaries.to_owned();
        e_now *= &hidden_all;
        e_now += &hidden_this;
        err = (&e_now - &e_prev).mapv(|x| x.abs()).iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    }
    e_now
}


#[derive(Copy, Clone, Debug)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Add for Point {
    type Output = Self;
    fn add(self, rhs: Self) -> Self { Self { x: self.x + rhs.x, y: self.y + rhs.y } }
}

impl AddAssign for Point {
    fn add_assign(&mut self, rhs: Self) { *self = Self { x: self.x + rhs.x, y: self.y + rhs.y }; }
}

impl Sub for Point {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self { Self { x: self.x - rhs.x, y: self.y - rhs.y} }
}

impl Mul for Point {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self { Self { x: self.x * rhs.x, y: self.y * rhs.y } }
}

impl Neg for Point {
    type Output = Self;
    fn neg(self) -> Self { Self { x: -self.x, y: -self.y } }
}

impl Div<f64> for Point {
    type Output = Self;
    fn div(self, rhs: f64) -> Self { Self { x: self.x / rhs, y: self.y / rhs } }
}

impl Mul<f64> for Point {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self { Self { x: self.x * rhs, y: self.y * rhs } }
}

impl Point {
    fn splatted(self, lim: Point, collidables: &Array2<i8>) -> bool {
        self.x < 1.0
            || self.y < 1.0
            || self.x >= lim.x
            || self.y >= lim.y
            || collidables[[self.y as usize, self.x as usize]] == 0
    }

    fn srss(self) -> f64 { (self.x * self.x + self.y * self.y).sqrt() }
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

/// speed in mm/us (m/ms), ccs in meters^2, pressure in Pa, gas_mass in amu
/// λ = k * T / (√2 * π * d² * p)
/// returns mfp in mm
fn calculate_mean_free_path(speed: f64, ccs: f64, pressure: f64, gas_mass: f64) -> f64 {
    let c_bar_gas = (8.0 * K * T / PI / (gas_mass * KG_AMU)).sqrt() / 1000.0;
    let c_star_gas = (2.0 * K * T / (gas_mass * KG_AMU)).sqrt() / 1000.0;
    let s = speed / c_star_gas;
    let c_bar_rel = c_bar_gas * (s + 1.0/(2.0*s)) * 0.5 * PI.sqrt() * erf(s) + 0.5 * (-s * s).exp();
    let mfp = K * T * (speed / c_bar_rel) / (pressure * ccs); // mfp
    mfp * 1000.0 // mfp in mm
}

/// calculates ccs in m^2
fn calculate_ccs(mass1: f64, mass2: f64) -> f64 {
    let ccs_rad1 = (mass_to_ccs(mass1) / PI).sqrt() * 1E-10;
    let ccs_rad2 = (mass_to_ccs(mass2) / PI).sqrt() * 1E-10;
    (ccs_rad1 + ccs_rad2).powi(2) * PI
}

/// emperical function that uses numbers from DOI: 10.3390/polym11040688
/// returns value in angstroms squared
fn mass_to_ccs(mass: f64) -> f64 {
    8.9 * (mass / 6.0).powf(2.0/3.0)
}

fn calculate_new_velocity(vel: Point, mass_amu: f64, gas_mass: f64, rng: &mut Rng) -> Point {
    let vel_stdev_gas = (K * T / (gas_mass * KG_AMU)).sqrt() / 1000.0;
    let gas_vel = Point{
        x: rng.gaussian() * vel_stdev_gas,
        y: rng.gaussian() * vel_stdev_gas,
    };
    // from https://en.wikipedia.org/wiki/Elastic_collision
    let v1 = (vel - gas_vel).srss(); // translate gas vel to zero & get magnitude
    let ion_angle = vel.y.atan2(vel.x);
    let (m1, m2, theta) = (mass_amu, gas_mass, 2.0 * PI * rng.f64() * 0.9999999999);
    let (costh, sinth, m1_2, m2_2) = (theta.cos(), theta.sin(), m1 * m1, m2 * m2);
    let delta_angle = (m2 * sinth).atan2(m1 + m2 * costh);
    let magnitude = v1 * (m1_2 + m2_2 + 2.0 * m1 * m2 * costh).sqrt() / (m1 + m2);
    // don't care about angle/magnitude for gas particle
    let angle = ion_angle + delta_angle;
    let vel_temp = Point{ x: magnitude * angle.cos(), y: magnitude * angle.sin() };
    let vel = vel_temp + gas_vel; // return frame of reference
    vel
}

/// Error function (erf) stolen shamelessly from SIMION's collision_hs1.lua file
//   erf(z) = (2/sqrt(pi)) * integral[0..z] exp(-t^2) dt
fn erf(z: f64) -> f64 {
    let z2 = z.abs();
    let t = 1.0 / (1.0 + 0.32759109962 * z2);
    let mut res = (    - 1.061405429 ) * t;
    res = (res + 1.453152027 ) * t;
    res = (res - 1.421413741 ) * t;
    res = (res + 0.2844966736) * t;
    res =((res - 0.254829592 ) * t) * (-z2*z2).exp();
    res += 1.0;
    if z < 0.0 {-res} else {res}
}


struct Rng(u64, u64);

impl Rng {
    const fn new(n: u64) -> Rng { Rng(n^0xf4dbdf2183dcefb7, n^0x1ad5be0d6dd28e9b) }

    #[inline]
    fn f64(&mut self) -> f64 {
        let (mut x, y) = (self.0, self.1);
        self.0 = y;
        x ^= x << 23;
        self.1 = x ^ y ^ (x >> 17) ^ (y >> 26);
        (self.1.wrapping_add(y) >> 32) as f64 * 2.3283064365386963E-10
    }

    #[inline]
    fn gaussian(&mut self) -> f64 {
        let (v1, v2) = (2.0 * self.f64() - 1.0, 2.0 * self.f64() - 1.0);
        let s = v1*v1 + v2*v2;
        if s < 1.0 && s != 0.0 { v1 * (-2.0 * s.ln() / s).sqrt() } else { self.gaussian() }
    }
}
