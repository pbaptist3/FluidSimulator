extern crate core;

mod app;
mod particle;

use std::cmp::max;
use sdl2::pixels::Color;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use std::time::{Duration, Instant};
use rayon::ThreadPoolBuilder;
use crate::app::App;
use crate::particle::Vec2;

const BG_COLOR: Color = Color::RGB(0, 0, 0);
const FRAME_RATE_CAP: i64 = 2000;
pub const PARTICLE_COLOR: Color = Color::RGB(0, 150, 0);
pub const WALL_PARTICLE_COLOR: Color = Color::RGB(0, 100, 0);
pub const PARTICLE_SIZE: u32 = 6;
pub const WINDOW_SIZE: (u32, u32) = (900, 600);

const MOUSE_DRAG_FORCE: f64 = 8.0;
const MOUSE_DRAG_RADIUS: f64 = 150.0;

const THREAD_COUNT: usize = 20;

pub fn main() {
    // init threadpool
    ThreadPoolBuilder::new()
        .num_threads(THREAD_COUNT)
        .build_global()
        .unwrap();

    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("Fluid Simulator", WINDOW_SIZE.0, WINDOW_SIZE.1)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();

    let mut app = App::new();

    let mut previous_draw = Instant::now();
    let mut is_paused = true;
    let mut is_framestep = false;
    canvas.set_draw_color(BG_COLOR);
    canvas.clear();
    canvas.present();
    let mut event_pump = sdl_context.event_pump().unwrap();
    'running: loop {
        canvas.set_draw_color(BG_COLOR);
        canvas.clear();
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} |
                Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                    break 'running
                },
                Event::KeyDown { keycode: Some(Keycode::Q), .. } => {
                    is_paused = !is_paused
                },
                Event::KeyDown { keycode: Some(Keycode::E), .. } => {
                    is_framestep = true
                },
                Event::MouseMotion { mousestate, x, y, xrel, yrel, .. } => {
                    let xrel = xrel as f64;
                    let yrel = yrel as f64;
                    let dist = (xrel*xrel+yrel*yrel).sqrt();
                    if mousestate.left() || mousestate.right() {
                        let direction = if mousestate.left() {1.0} else {-1.0};
                        app.force_points_in_radius(
                            MOUSE_DRAG_FORCE * direction * dist,
                            Vec2 {x: x.into(), y: y.into()},
                            Vec2 {x: xrel, y: yrel} * (1.0 / dist),
                            MOUSE_DRAG_RADIUS,
                        );
                    }
                }
                _ => {}
            }
        }
        // main game loop stuffs
        let time_now = Instant::now();
        let delta_duration = time_now.duration_since(previous_draw);
        let do_physics = (!is_paused) || is_framestep;
        app.draw(&mut canvas, delta_duration.as_secs_f64(), do_physics);
        is_framestep = false;
        previous_draw = time_now;

        canvas.present();
        let sleep_duration = Duration::new(
            0,
            max(1_000_000_000i64 / FRAME_RATE_CAP - delta_duration.as_nanos() as i64, 0) as u32
        );
        std::thread::sleep(sleep_duration);
    }
}