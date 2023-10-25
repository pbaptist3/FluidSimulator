use std::ops::{Add, Mul, Sub};

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

impl Add for Vec2 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Vec2 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub for Vec2 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Vec2 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Mul<f64> for Vec2 {
    type Output = Vec2;

    fn mul(self, rhs: f64) -> Self::Output {
        Vec2 {
            x: self.x * rhs,
            y: self.y * rhs
        }
    }
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Particle {
    pub position: Vec2,
    pub velocity: Vec2,
    is_wall: bool
}

impl Particle {
    pub fn at_rest(position: (f64, f64), is_wall: bool) -> Self {
        Self {
            position: Vec2 {x: position.0, y: position.1},
            velocity: Vec2 {x: 0.0, y: 0.0},
            is_wall
        }
    }

    pub fn rand_velocity(position: (f64, f64)) -> Self {
        Self {
            position: Vec2 {x: position.0, y: position.1},
            velocity: Vec2 {x: rand::random::<f64>() * 400.0 - 200.0, y: rand::random::<f64>() * 200.0 - 100.0},
            is_wall: false,
        }
    }

    pub fn is_wall(&self) -> bool {
        self.is_wall
    }
}