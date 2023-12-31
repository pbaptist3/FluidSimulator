use std::f64::consts::PI;
use sdl2::render::Canvas;
use sdl2::video::Window;
use sdl2::rect::Rect;
use sdl2::rect::Point;
use sdl2::pixels::Color;
use rayon::prelude::*;
use crate::particle::{Particle, Vec2};
use crate::{PARTICLE_SIZE, WALL_PARTICLE_COLOR, WINDOW_SIZE};

const PARTICLE_COUNT: u32 = 400;
const PARTICLE_SPACING: f64 = 15.0;
const GRAVITY: f64 = 250.0;
const DAMPENING: f64 = 1.2;
const SMOOTHING_RADIUS: f64 = 40.0;
const PARTICLE_MASS: f64 = 1.0;
const TARGET_DENSITY: f64 = 0.2;
const DEBUG_TARGET: usize = 10000000;
const GAS_CONSTANT: f64 = 0.004;
const VISCOSITY: f64 = 400000.0;
const WALL_FRICTION: f64 = 1.10;
const SURFACE_TENSION: f64 = 0.002;
const MIN_SURFACE_TENSION: f64 = 1.0;
const WALL_PARTICLE_STRENGTH: f64 = 3.0;
const WALL_PARTICLE_SPACING: f64 = 10.0;
const WALL_SMOOTHING_RADIUS: f64 = 10.0;

pub struct App {
    particles: Vec<Particle>,
    grid: Vec<Vec<Vec<usize>>>,
}

impl App {
    pub fn new() -> Self {
        let mut points = Self::generate_particle_array();
        let grid_hor_size = WINDOW_SIZE.0 as usize / SMOOTHING_RADIUS as usize + 2;
        let grid_vert_size = WINDOW_SIZE.1 as usize / SMOOTHING_RADIUS as usize + 2;
        let grid = vec![vec![Vec::new(); grid_hor_size]; grid_vert_size];
        Self::generate_wall_particles(&mut points);

        Self {
            particles: points,
            grid,
        }
    }

    // update physics states
    fn phys_update(&mut self, delta_time: f64) {
        for particle in self.particles.iter_mut().filter(|p| !p.is_wall()) {
            particle.position.x += particle.velocity.x * delta_time;
            particle.position.y += particle.velocity.y * delta_time;
            particle.velocity.y += GRAVITY * delta_time;


            // we need some bounds
            if particle.position.y > WINDOW_SIZE.1 as f64 || particle.position.y < 0.0
            {
                particle.position.y = if particle.position.y < 0.0 {
                    0.0
                } else {
                    WINDOW_SIZE.1 as f64
                };

                particle.velocity.y *= -1.0 / DAMPENING;
                particle.velocity.x /= WALL_FRICTION;
            }

            if particle.position.x > WINDOW_SIZE.0 as f64 || particle.position.x < 0.0 {
                particle.position.x = if particle.position.x < 0.0 {
                    0.0
                } else {
                    WINDOW_SIZE.0 as f64
                };
                particle.velocity.x *= -1.0 / DAMPENING;
                // EXPERIMENTAL maybe stops the wall climbing?
                particle.velocity.y /= WALL_FRICTION;
            }
        }

        // update grid of particles
        self.update_grid();

        // calculate densities
        let densities: Vec<f64> = self.particles.iter()
            .filter(|p| !p.is_wall())
            .map(|particle| self.calculate_density(particle))
            .collect();

        // get pressures and forces
        let pressures = self.calculate_pressures(&densities);
        let mut pressure_forces = self.calculate_pressure_forces(&pressures, &densities);

        // get viscosity forces
        let viscosity_forces = self.calculate_viscosities(&densities);

        // get surface tension forces
        let surface_tension_forces = self.calculate_surface_tensions(&densities);

        // get wall forces
        let wall_forces = self.particles.iter()
            .filter(|p| p.is_wall())
            .collect::<Vec<_>>()
            .par_iter()
            .map(|wall_particle| {
                self.iter_nearby_particles(wall_particle.position, false)
                    .iter()
                    .map(|(i, particle)| {
                        let delta = particle.position - wall_particle.position;
                        let dist = (delta.x*delta.x + delta.y*delta.y).sqrt();
                        if dist > WALL_SMOOTHING_RADIUS {
                            (*i, Vec2 {x: 0.0, y: 0.0})
                        } else {
                            (*i, delta * (WALL_PARTICLE_STRENGTH / dist * (WALL_SMOOTHING_RADIUS - dist).powi(2)))
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        for (i, force) in wall_forces.iter().flatten() {
            let particle = &mut self.particles[*i];
            particle.velocity = particle.velocity + *force;
        }

        // apply all the forces
        self.particles.iter_mut()
            .filter(|p| !p.is_wall())
            .collect::<Vec<_>>()
            .par_iter_mut()
            .zip(pressure_forces.par_iter())
            .zip(viscosity_forces.par_iter())
            .zip(surface_tension_forces.par_iter())
            .enumerate()
            .for_each(|(i, (((point, pressure), viscosity), tension))| {
                let p1 = point;
                p1.velocity = p1.velocity
                    + (*pressure + *viscosity + *tension) * (1.0 / densities[i]);

                if i == DEBUG_TARGET {
                    println!("gradient: {:?}\tdensity: {:?}", pressure, densities[i]);
                    println!("viscosity: {:?}", viscosity);
                    println!("surface tension: {:?}", tension);
                    println!();
                }
            });
    }

    fn calculate_pressures(&self, densities: &[f64]) -> Vec<f64> {
        densities.iter()
            .map(|density| GAS_CONSTANT * (density - TARGET_DENSITY))
            .collect()
    }

    fn calculate_pressure_forces(&self, pressures: &[f64], densities: &[f64]) -> Vec<Vec2> {
        self.particles.iter()
            .filter(|p| !p.is_wall())
            .collect::<Vec<_>>()
            .par_iter()
            .zip(pressures)
            .zip(densities)
            .map(|((particle, pressure), density)|
                self.calculate_gradient(particle, densities, |p2, i| (*pressure + pressures[i]) / 2.0) * (1.0 / density)
            )
            .collect()
    }

    fn calculate_viscosities(&self, densities: &[f64]) -> Vec<Vec2> {
        self.particles.iter()
            .filter(|p| !p.is_wall())
            .collect::<Vec<_>>()
            .par_iter()
            .zip(densities)
            .map(|(particle, density)|
                self.calculate_laplacian(particle, densities, |p2, i| (p2.velocity - particle.velocity) * 0.5)
                    * VISCOSITY
            )
            .collect()
    }

    fn calculate_surface_tensions(&self, densities: &[f64]) -> Vec<Vec2> {
        self.particles.iter()
            .filter(|p| !p.is_wall())
            .collect::<Vec<_>>()
            .par_iter()
            .zip(densities)
            .map(|(particle, density)| {
                let gradient = self.calculate_gradient(particle, densities, |_, _| 1.0);
                let laplacian = self.calculate_laplacian(particle, densities, |p2, _|  {
                    let delta = p2.position - particle.position;
                    let d = (delta.x*delta.x + delta.y*delta.y).sqrt();
                    if d < MIN_SURFACE_TENSION {
                        Vec2 {x: 0.0, y: 0.0}
                    } else {
                        delta
                    }
                });
                let grad_magnitude = (gradient.x*gradient.x + gradient.y*gradient.y).sqrt();
                let laplacian = laplacian.x + laplacian.y;
                gradient * (-1.0 * SURFACE_TENSION * laplacian * grad_magnitude)
            })
            .collect()
    }

    fn calculate_property<F>(&self, particle: &Particle, densities: &[f64], property: F) -> f64
    where
        F: Fn(&Particle, usize) -> f64 + std::marker::Sync + Send
    {
        self.iter_nearby_particles(particle.position, false)
            .par_iter()
            .zip(densities)
            .enumerate()
            .fold(|| 0.0f64, |sum: f64, (i, ((_, particle2), density))| {
                // skip if particles are the same
                if particle.position == particle2.position {
                    return sum;
                }
                let property_value = property(particle2, i);
                let result = smoothing_kernel(&particle, &particle2) * (property_value * PARTICLE_MASS / density);
                sum + result
            })
            .sum()
    }

    fn calculate_gradient<F>(&self, particle: &Particle, densities: &[f64], property: F) -> Vec2
        where
            F: Fn(&Particle, usize) -> f64 + std::marker::Sync
    {
        self.iter_nearby_particles(particle.position, false)
            .par_iter()
            .zip(densities)
            .enumerate()
            .fold(|| Vec2 {x: 0.0, y: 0.0}, |sum, (i, ((_, particle2), density))| {
                // skip if particles are the same
                if *particle == **particle2 {
                    return sum;
                }

                let property_value = property(particle2, i);
                let result = if *density == 0.0 {
                    Vec2 {x: 0.0, y: 0.0}
                } else if particle.position == particle2.position {
                    Vec2 {
                        x: rand::random::<f64>()*0.1 - 0.05,
                        y: rand::random::<f64>()*0.1 - 0.05,
                    }
                } else {
                    kernel_gradient(&particle, &particle2) * (property_value * PARTICLE_MASS / density)
                };
                sum + result
            })
            .reduce(|| Vec2 {x: 0.0, y: 0.0}, |a, b| a+b)
    }

    fn calculate_laplacian<F>(&self, particle: &Particle, densities: &[f64], property: F) -> Vec2
        where
            F: Fn(&Particle, usize) -> Vec2 + std::marker::Sync
    {
        self.iter_nearby_particles(particle.position, false)
            .par_iter()
            .zip(densities)
            .enumerate()
            .fold(|| Vec2 {x: 0.0, y: 0.0}, |sum, (i, ((_, particle2), density))| {
                // skip if particles are the same
                if particle.position == particle2.position {
                    return sum;
                }
                let property_value = property(particle2, i);
                let result = if *density == 0.0 {
                    Vec2 {x: 0.0, y: 0.0}
                } else {
                    property_value * (kernel_laplacian(&particle, &particle2) * PARTICLE_MASS / density)
                };
                sum + result
            })
            .reduce(|| Vec2 {x: 0.0, y: 0.0}, |a,b| a+b)
    }

    /// calculates density by summing smoothing kernel * mass over all particles
    fn calculate_density(&self, particle: &Particle) -> f64 {
        let mut sum = 0.0;
        for (_, particle2) in self.iter_nearby_particles(particle.position, false) {
            let result = PARTICLE_MASS * smoothing_kernel(&particle, &particle2);
            sum += result;
        }
        sum
    }

    // update physics then draw
    pub fn draw(&mut self, canvas: &mut Canvas<Window>, delta_time: f64, is_paused: bool) {
        // dont do physics updates if paused
        if !is_paused {
            self.phys_update(delta_time);
        }

        // draw stuff
        for (i, particle) in self.particles.iter().enumerate() {
            // set color based on particle velocity
            let color_r = (particle.velocity.x*particle.velocity.x * particle.velocity.y*particle.velocity.y).sqrt() * 0.5;
            let color_b = 255.0 - color_r;
            if particle.is_wall() {
                canvas.set_draw_color(WALL_PARTICLE_COLOR);
            } else {
                canvas.set_draw_color(Color::RGB(color_r as u8, 0, color_b as u8));
            }

            // if there is a debug target, do some special stuff for it
            if i == DEBUG_TARGET {
                canvas.set_draw_color(Color::RGB(0, 255, 0));

                // influence square (it should be circle)
                let influence = Rect::from_center(
                    Point::new(particle.position.x as i32, particle.position.y as i32), 2*SMOOTHING_RADIUS as u32, 2*SMOOTHING_RADIUS as u32
                );
                canvas.draw_rect(influence);
            }
            //canvas.set_draw_color(Color::RGB(200, 0, 0));
            let rect = Rect::from_center(
                Point::new(particle.position.x as i32, particle.position.y as i32), PARTICLE_SIZE, PARTICLE_SIZE
            );
            canvas.fill_rect(rect);
        }
    }

    // create a starting grid of particles
    fn generate_particle_array() -> Vec<Particle> {
        let mut particles = Vec::new();

        let midpoint = (WINDOW_SIZE.0 as f64 / 2.0, WINDOW_SIZE.1 as f64 / 4.0 * 3.0);
        // generate a relatively square grid of particles at start
        let height = PARTICLE_COUNT * PARTICLE_SPACING as u32 / (WINDOW_SIZE.0) + 1;
        //let height = (PARTICLE_COUNT as f64 * WINDOW_SIZE.0 as f64 / PARTICLE_SPACING - 10.0).floor() as u32;
        let width = PARTICLE_COUNT / height + 1;
        'particle_loop: for row in 0..height {
            for col in 0..width {
                // exit if we have enough particles
                if particles.len() >= PARTICLE_COUNT as usize {
                    break 'particle_loop;
                }

                let x = midpoint.0 - (width as f64 / 2.0 - col as f64) * PARTICLE_SPACING;
                let y = midpoint.1 - (height as f64 / 2.0 - row as f64) * PARTICLE_SPACING;
                let particle = Particle::rand_velocity((x, y));
                particles.push(particle);
            }
        }

        particles
    }

    fn generate_wall_particles(particles: &mut Vec<Particle>) {
        for x in (0..WINDOW_SIZE.0).step_by(WALL_PARTICLE_SPACING as usize) {
            let top_particle = Particle::at_rest((x as f64, 0.0), true);
            let bottom_particle = Particle::at_rest((x as f64, WINDOW_SIZE.1 as f64), true);
            particles.push(top_particle);
            particles.push(bottom_particle);
        }

        for y in (0..WINDOW_SIZE.1).step_by(WALL_PARTICLE_SPACING as usize) {
            let top_particle = Particle::at_rest((0.0, y as f64), true);
            let bottom_particle = Particle::at_rest((WINDOW_SIZE.0 as f64, y as f64), true);
            particles.push(top_particle);
            particles.push(bottom_particle);
        }
    }

    pub fn force_points_in_radius(&mut self, force_size: f64, position: Vec2, direction: Vec2, radius: f64) {
        for particle in self.particles.iter_mut() {
            let delta = particle.position - position;
            let distance = (delta.x*delta.x + delta.y*delta.y).sqrt();
            if distance > radius {
                continue;
            }
            particle.velocity = /*particle.velocity +*/ (direction * force_size);
        }
    }

    fn iter_nearby_particles(&self, position: Vec2, include_wall: bool) -> Vec<(usize, &Particle)> {
        let smoothing_radius = if include_wall {
            WALL_SMOOTHING_RADIUS
        } else {
            SMOOTHING_RADIUS
        };
        let grid_x = (position.x / smoothing_radius).floor() as isize;
        let grid_y = (position.y / smoothing_radius).floor() as isize;

        let mut indices: Vec<usize> = Vec::new();
        for y in isize::max(grid_y-1, 0)..=isize::min(grid_y+1, (self.grid.len() - 1) as isize) {
            for x in isize::max(grid_x-1, 0)..=isize::min(grid_x+1, (self.grid[0].len() - 1) as isize) {
                self.grid[y as usize][x as usize].iter()
                    .for_each(|index| indices.push(*index));
            }
        }
        indices.iter()
            .map(|i| (*i, &self.particles[*i]))
            .filter(|(_, p)| p.is_wall() == include_wall)
            .collect::<Vec<(usize, &Particle)>>()
    }

    fn update_grid(&mut self) {
        // clear old grid
        self.grid.iter_mut()
            .map(|row| row.iter_mut())
            .flatten()
            .for_each(|grid| grid.clear());

        // add each particle to correct grid
        for (i, particle) in self.particles.iter().enumerate() {
            let grid_x = (particle.position.x / SMOOTHING_RADIUS).floor() as usize;
            let grid_y = (particle.position.y / SMOOTHING_RADIUS).floor() as usize;

            self.grid[grid_y][grid_x].push(i);
        }
    }
}

fn smoothing_kernel(p1: &Particle, p2: &Particle) -> f64 {
    let delta = p2.position - p1.position;
    let r_sq = delta.x*delta.x + delta.y*delta.y;
    // discard particles that are too far away
    if r_sq > SMOOTHING_RADIUS.powi(2) {
        return 0.0;
    }

    315.0 / 64.0 / PI / SMOOTHING_RADIUS.powi(9) * (SMOOTHING_RADIUS.powi(2) - r_sq).powi(3) * 20000.0
}

fn kernel_gradient(p1: &Particle, p2: &Particle) -> Vec2 {
    let delta = p2.position - p1.position;
    let distance = f64::sqrt(delta.x*delta.x + delta.y*delta.y);
    // discard particles that are too far away
    if distance > SMOOTHING_RADIUS {
        return Vec2 {x: 0.0, y: 0.0};
    }

    Vec2 {
        x: -1.0 * 3.0 * delta.x * (SMOOTHING_RADIUS - distance).powi(2) / distance,
        y: -1.0 * 3.0 * delta.y * (SMOOTHING_RADIUS - distance).powi(2) / distance,
    }
}

fn kernel_laplacian(p1: &Particle, p2: &Particle) -> f64 {
    let delta = p2.position - p1.position;
    let r = f64::sqrt(delta.x*delta.x + delta.y*delta.y);
    // discard particles that are too far away
    if r > SMOOTHING_RADIUS {
        return 0.0;
    }

    45.0 / PI / SMOOTHING_RADIUS.powi(6) * (SMOOTHING_RADIUS - r)
}