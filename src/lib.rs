mod utils;

use wasm_bindgen::prelude::*;
use std::collections::VecDeque;

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

const F9TH: f64 = 4.0/9.0;
const O9TH: f64 = 1.0/9.0;
const O36TH: f64 = 1.0/36.0;

#[wasm_bindgen]
pub struct Simulator {
    width: usize,
    height: usize,
    omega: f64,
    u0: f64,
    n_0: VecDeque<VecDeque<f64>>,
    n_n: VecDeque<VecDeque<f64>>,
    n_s: VecDeque<VecDeque<f64>>,
    n_e: VecDeque<VecDeque<f64>>,
    n_w: VecDeque<VecDeque<f64>>,
    n_ne: VecDeque<VecDeque<f64>>,
    n_nw: VecDeque<VecDeque<f64>>,
    n_se: VecDeque<VecDeque<f64>>,
    n_sw: VecDeque<VecDeque<f64>>,
    rho: Vec<Vec<f64>>,
    ux: Vec<Vec<f64>>,
    uy: Vec<Vec<f64>>,
    curl: Vec<f64>,
    barrier: VecDeque<VecDeque<bool>>,
    flat_bar: Vec<u8>,
    fn0: Vec<f64>,
    fnn: Vec<f64>,
    fns: Vec<f64>,
    fne: Vec<f64>,
    fnw: Vec<f64>,
    fnne: Vec<f64>,
    fnnw: Vec<f64>,
    fnse: Vec<f64>,
    fnsw: Vec<f64>,
}

impl Simulator {
    fn stream(&mut self) {
        self.n_n.rotate_right(1);
        self.n_n[self.height - 1] = self.n_n[0].clone();
        self.n_s.rotate_left(1);
        self.n_s[0] = self.n_s[self.height - 1].clone();

        self.n_ne.rotate_right(1);
        self.n_ne[self.height - 1] = self.n_ne[0].clone();
        self.n_nw.rotate_right(1);
        self.n_nw[self.height - 1] = self.n_nw[0].clone();
        self.n_se.rotate_left(1);
        self.n_se[0] = self.n_se[self.height - 1].clone();
        self.n_sw.rotate_left(1);
        self.n_sw[0] = self.n_sw[self.height - 1].clone();
        for i in 0..self.height {
            self.n_e[i].rotate_right(1);
            self.n_e[i][self.width - 1] = self.n_e[i][0].clone();
            self.n_w[i].rotate_left(1);
            self.n_w[i][0] = self.n_w[i][self.width - 1].clone();
            self.n_ne[i].rotate_right(1);
            self.n_ne[i][self.width - 1] = self.n_ne[i][0].clone();
            self.n_se[i].rotate_right(1);
            self.n_se[i][self.width - 1] = self.n_se[i][0].clone();
            self.n_nw[i].rotate_left(1);
            self.n_nw[i][0] = self.n_nw[i][self.width - 1].clone();
            self.n_sw[i].rotate_left(1);
            self.n_sw[i][0] = self.n_sw[i][self.width - 1].clone();
        }

        for i in 0..self.height {
            for j in 0..self.width {
                if self.barrier[i][j] {
                    self.n_n[i+1][j] = self.n_s[i][j].clone();
                    self.n_s[i-1][j] = self.n_n[i][j].clone();
                    self.n_e[i][j+1] = self.n_w[i][j].clone();
                    self.n_w[i][j-1] = self.n_e[i][j].clone();
                    self.n_ne[i+1][j+1] = self.n_sw[i][j].clone();
                    self.n_nw[i+1][j-1] = self.n_se[i][j].clone();
                    self.n_se[i-1][j+1] = self.n_nw[i][j].clone();
                    self.n_sw[i-1][j-1] = self.n_ne[i][j].clone();
                }
            }
        }
    }

    fn collide(&mut self) {
        for i in 0..self.height {
            for j in 0..self.width {
                self.rho[i][j] = self.n_0[i][j] + self.n_e[i][j] + self.n_n[i][j]
                               + self.n_w[i][j] + self.n_s[i][j] + self.n_ne[i][j]
                               + self.n_nw[i][j] + self.n_sw[i][j] + self.n_se[i][j];
                self.ux[i][j] = (
                    self.n_e[i][j] + self.n_ne[i][j] + self.n_se[i][j]
                  - self.n_w[i][j] - self.n_nw[i][j] - self.n_sw[i][j]
                ) / self.rho[i][j];
                self.uy[i][j] = (
                    self.n_n[i][j] + self.n_ne[i][j] + self.n_nw[i][j]
                  - self.n_s[i][j] - self.n_se[i][j] - self.n_sw[i][j]
                ) / self.rho[i][j];

                let ux2 = self.ux[i][j].powi(2);
                let uy2 = self.uy[i][j].powi(2);
                let u2 = ux2 + uy2;
                let omu215 = 1.0 - 1.5 * u2;
                let uxuy = self.ux[i][j] * self.uy[i][j];

                self.n_0[i][j] = (1.0 - self.omega) * self.n_0[i][j]
                              + self.omega * F9TH * self.rho[i][j]
                              * omu215;
                self.n_n[i][j] = (1.0 - self.omega) * self.n_n[i][j]
                              + self.omega * O9TH * self.rho[i][j]
                              * (omu215 + 3.0*self.uy[i][j] + 4.5*uy2);
                self.n_s[i][j] = (1.0 - self.omega) * self.n_s[i][j]
                              + self.omega * O9TH * self.rho[i][j]
                              * (omu215 - 3.0*self.uy[i][j] + 4.5*uy2);
                self.n_e[i][j] = (1.0 - self.omega) * self.n_e[i][j]
                              + self.omega * O9TH * self.rho[i][j]
                              * (omu215 + 3.0*self.ux[i][j] + 4.5*ux2);
                self.n_w[i][j] = (1.0 - self.omega) * self.n_w[i][j]
                              + self.omega * O9TH * self.rho[i][j]
                              * (omu215 - 3.0*self.ux[i][j] + 4.5*ux2);
                self.n_ne[i][j] = (1.0 - self.omega) * self.n_ne[i][j]
                              + self.omega * O36TH * self.rho[i][j]
                            * (
                                omu215
                              + 3.0 * (self.ux[i][j] + self.uy[i][j])
                              + 4.5 * (u2 + 2.0*uxuy)
                              );
                self.n_nw[i][j] = (1.0 - self.omega) * self.n_nw[i][j]
                              + self.omega * O36TH * self.rho[i][j]
                            * (
                                omu215
                              + 3.0 * (-self.ux[i][j] + self.uy[i][j])
                              + 4.5 * (u2 - 2.0*uxuy)
                              );
                self.n_se[i][j] = (1.0 - self.omega) * self.n_se[i][j]
                              + self.omega * O36TH * self.rho[i][j]
                            * (
                                omu215
                              + 3.0 * (self.ux[i][j] - self.uy[i][j])
                              + 4.5 * (u2 - 2.0*uxuy)
                              );
                self.n_sw[i][j] = (1.0 - self.omega) * self.n_sw[i][j]
                              + self.omega * O36TH * self.rho[i][j]
                            * (
                                omu215
                              + 3.0 * (-self.ux[i][j] - self.uy[i][j])
                              + 4.5 * (u2 + 2.0*uxuy)
                              );
                              
                self.n_0[0][j] = F9TH * (1.0 - 1.5*self.u0.powi(2));
                self.n_n[0][j] = O9TH * (1.0 - 1.5*self.u0.powi(2));
                self.n_s[0][j] = O9TH * (1.0 - 1.5*self.u0.powi(2));
                self.n_e[0][j] = O9TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_w[0][j] = O9TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_ne[0][j] = O36TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_se[0][j] = O36TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_nw[0][j] = O36TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_sw[0][j] = O36TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_0[self.height - 1][j] = F9TH * (1.0 - 1.5*self.u0.powi(2));
                self.n_n[self.height - 1][j] = O9TH * (1.0 - 1.5*self.u0.powi(2));
                self.n_s[self.height - 1][j] = O9TH * (1.0 - 1.5*self.u0.powi(2));
                self.n_e[self.height - 1][j] = O9TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_w[self.height - 1][j] = O9TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_ne[self.height - 1][j] = O36TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_se[self.height - 1][j] = O36TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_nw[self.height - 1][j] = O36TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
                self.n_sw[self.height - 1][j] = O36TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
            }
            self.n_0[i][0] = F9TH * (1.0 - 1.5*self.u0.powi(2));
            self.n_n[i][0] = O9TH * (1.0 - 1.5*self.u0.powi(2));
            self.n_s[i][0] = O9TH * (1.0 - 1.5*self.u0.powi(2));
            self.n_e[i][0] = O9TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_w[i][0] = O9TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_ne[i][0] = O36TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_se[i][0] = O36TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_nw[i][0] = O36TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_sw[i][0] = O36TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_0[i][self.width - 1] = F9TH * (1.0 - 1.5*self.u0.powi(2));
            self.n_n[i][self.width - 1] = O9TH * (1.0 - 1.5*self.u0.powi(2));
            self.n_s[i][self.width - 1] = O9TH * (1.0 - 1.5*self.u0.powi(2));
            self.n_e[i][self.width - 1] = O9TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_w[i][self.width - 1] = O9TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_ne[i][self.width - 1] = O36TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_se[i][self.width - 1] = O36TH * (1.0 + 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_nw[i][self.width - 1] = O36TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
            self.n_sw[i][self.width - 1] = O36TH * (1.0 - 3.0*self.u0 + 3.0*self.u0.powi(2));
        }
    }

    fn get_curl(&self) -> Vec<f64> {
        let mut a = self.uy.clone();
        let mut b = self.uy.clone();
        let mut c = self.ux.clone();
        let mut d = self.ux.clone();
        let mut curl = self.curl.clone();
        for i in 0..self.height {
            a[i].rotate_left(1);
            b[i].rotate_right(1);
        }
        c.rotate_left(1);
        d.rotate_right(1);
        for i in 1..(self.height - 1) {
            for j in 1..(self.width - 1) {
                curl[i * self.width + j] = a[i][j] - b[i][j] - c[i][j] + d[i][j];
            }
        }

        curl
    }
}

#[wasm_bindgen]
impl Simulator {
    pub fn new(w: usize, h: usize, vis: f64, u_not: f64) -> Simulator {
        utils::set_panic_hook();

        let width = w;
        let height = h;
        let omega = 1.0 / (3.0 * vis + 0.5);
        let u0 = u_not;

        let n_0: VecDeque<VecDeque<f64>> = vec![
            vec![F9TH*(1.0-1.5*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        let n_n: VecDeque<VecDeque<f64>> = vec![
            vec![O9TH*(1.0-1.5*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        let n_s: VecDeque<VecDeque<f64>> = vec![
            vec![O9TH*(1.0-1.5*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        let n_e: VecDeque<VecDeque<f64>> = vec![
            vec![O9TH*(1.0+3.0*u0+3.0*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        let n_w: VecDeque<VecDeque<f64>> = vec![
            vec![O9TH*(1.0-3.0*u0+3.0*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        let n_ne: VecDeque<VecDeque<f64>> = vec![
            vec![O36TH*(1.0+3.0*u0+3.0*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        let n_se: VecDeque<VecDeque<f64>> = vec![
            vec![O36TH*(1.0+3.0*u0+3.0*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        let n_nw: VecDeque<VecDeque<f64>> = vec![
            vec![O36TH*(1.0-3.0*u0+3.0*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        let n_sw: VecDeque<VecDeque<f64>> = vec![
            vec![O36TH*(1.0-3.0*u0+3.0*u0.powi(2)); width].into_iter().collect(); height
        ].into_iter().collect();
        
        let fn0 = vec![F9TH*(1.0-1.5*u0.powi(2)); width * height];
        let fnn = vec![O9TH*(1.0-1.5*u0.powi(2)); width * height];
        let fns = vec![O9TH*(1.0-1.5*u0.powi(2)); width * height];
        let fne = vec![O9TH*(1.0+3.0*u0+3.0*u0.powi(2)); width * height];
        let fnw = vec![O9TH*(1.0-3.0*u0+3.0*u0.powi(2)); width * height];
        let fnne = vec![O36TH*(1.0+3.0*u0+3.0*u0.powi(2)); width * height];
        let fnnw = vec![O36TH*(1.0-3.0*u0+3.0*u0.powi(2)); width * height];
        let fnse = vec![O36TH*(1.0+3.0*u0+3.0*u0.powi(2)); width * height];
        let fnsw = vec![O36TH*(1.0-3.0*u0+3.0*u0.powi(2)); width * height];

        let mut rho: Vec<Vec<f64>> = vec![Vec::with_capacity(width); height];
        let mut ux: Vec<Vec<f64>> = vec![Vec::with_capacity(width); height];
        let mut uy: Vec<Vec<f64>> = vec![Vec::with_capacity(width); height];
        for i in 0..height {
            for j in 0..width {
                rho[i].push(n_0[i][j] + n_e[i][j] + n_n[i][j] + n_w[i][j] + n_s[i][j]
                    + n_ne[i][j] + n_nw[i][j] + n_sw[i][j] + n_se[i][j]);
                ux[i].push((n_e[i][j]+n_ne[i][j]+n_se[i][j]-n_w[i][j]-n_nw[i][j]-n_sw[i][j])
                    / rho[i][j]);
                uy[i].push((n_n[i][j]+n_ne[i][j]+n_nw[i][j]-n_s[i][j]-n_se[i][j]-n_sw[i][j])
                    / rho[i][j]);
            }
        }

        let mut barrier: VecDeque<VecDeque<bool>> = vec![
            vec![false; width].into_iter().collect(); height
        ].into_iter().collect();
        barrier[height / 2][width / 2] = true;
        let mut reference: Vec<Vec<bool>> = vec![vec![false; width]; height];
        reference[height / 2][width / 2] = true;
        let mut flat_bar: Vec<u8> = Vec::with_capacity(width * height);
        for i in 0..height {
            for j in 0..width {
                flat_bar.push(reference[i][j] as u8);
            }
        }

        let mut a = uy.clone();
        let mut b = uy.clone();
        let mut c = ux.clone();
        let mut d = ux.clone();
        let mut curl: Vec<f64> = Vec::with_capacity(width * height);
        for i in 0..height {
            a[i].rotate_left(1);
            b[i].rotate_right(1);
        }
        c.rotate_left(1);
        d.rotate_right(1);
        for i in 0..height {
            for j in 0..width {
                curl.push(a[i][j] - b[i][j] - c[i][j] + d[i][j]);
            }
        }

        Simulator {
            width,
            height,
            omega,
            u0,
            n_0,
            n_n,
            n_s,
            n_e,
            n_w,
            n_ne,
            n_nw,
            n_se,
            n_sw,
            rho,
            ux,
            uy,
            curl,
            barrier,
            flat_bar,
            fn0,
            fnn,
            fns,
            fne,
            fnw,
            fnne,
            fnnw,
            fnse,
            fnsw,
        }
    }

    pub fn tick(&mut self) {
        self.stream();
        self.collide();
        self.curl = self.get_curl();
    }

    pub fn n_0(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fn0[i * self.width + j] = self.n_0[i][j];
            }
        }
        self.fn0.as_ptr()
    }
    pub fn n_n(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fnn[i * self.width + j] = self.n_n[i][j];
            }
        }
        self.fnn.as_ptr()
    }
    pub fn n_s(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fns[i * self.width + j] = self.n_s[i][j];
            }
        }
        self.fns.as_ptr()
    }
    pub fn n_e(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fne[i * self.width + j] = self.n_e[i][j];
            }
        }
        self.fne.as_ptr()
    }
    pub fn n_w(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fnw[i * self.width + j] = self.n_w[i][j];
            }
        }
        self.fnw.as_ptr()
    }
    pub fn n_ne(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fnne[i * self.width + j] = self.n_ne[i][j];
            }
        }
        self.fnne.as_ptr()
    }
    pub fn n_nw(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fnnw[i * self.width + j] = self.n_nw[i][j];
            }
        }
        self.fnnw.as_ptr()
    }
    pub fn n_se(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fnse[i * self.width + j] = self.n_se[i][j];
            }
        }
        self.fnse.as_ptr()
    }
    pub fn n_sw(&mut self) -> *const f64 {
        for i in 0..self.height {
            for j in 0..self.width {
                self.fnsw[i * self.width + j] = self.n_sw[i][j];
            }
        }
        self.fnsw.as_ptr()
    }

    pub fn barrier(&self) -> *const u8 {
        self.flat_bar.as_ptr()
    }

    pub fn curl(&self) -> *const f64 {
        self.curl.as_ptr()
    }
}