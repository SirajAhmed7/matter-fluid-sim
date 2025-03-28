const canvas = document.getElementById('canvas');
const gl = canvas.getContext('webgl');
canvas.width = window.innerWidth - 20;
canvas.height = window.innerHeight - 20;

canvas.focus();

const simHeight = 3;
const cScale = canvas.height / simHeight;
const simWidth = canvas.width / cScale;

const U_FIELD = 0;
const V_FIELD = 1;

const FLUID_CELL = 0;
const AIR_CELL = 1;
const SOLID_CELL = 2;

let cnt = 0;

const clamp = (x, min, max) => Math.max(min, Math.min(max, x));

// ----------------- start of simulator ------------------------------

class FlipFluid {
  constructor(density, width, height, spacing, particleRadius, maxParticles) {
    // fluid
    this.density = density;
    this.fNumX = Math.floor(width / spacing) + 1;
    this.fNumY = Math.floor(height / spacing) + 1;
    this.h = Math.max(width / this.fNumX, height / this.fNumY);
    this.fInvSpacing = 1.0 / this.h;
    this.fNumCells = this.fNumX * this.fNumY;

    this.u = new Float32Array(this.fNumCells);
    this.v = new Float32Array(this.fNumCells);
    this.du = new Float32Array(this.fNumCells);
    this.dv = new Float32Array(this.fNumCells);
    this.prevU = new Float32Array(this.fNumCells);
    this.prevV = new Float32Array(this.fNumCells);
    this.p = new Float32Array(this.fNumCells);
    this.s = new Float32Array(this.fNumCells);
    this.cellType = new Int32Array(this.fNumCells);
    this.cellColor = new Float32Array(3 * this.fNumCells);

    // particles
    this.maxParticles = maxParticles;

    this.particlePos = new Float32Array(2 * this.maxParticles);
    this.particleColor = new Float32Array(3 * this.maxParticles);
    for (let i = 0; i < this.maxParticles; i++) {
      // this.particleColor[3 * i + 2] = 1.0;
      this.particleColor[3 * i] = 1.0;
      this.particleColor[3 * i + 1] = 1.0;
      this.particleColor[3 * i + 2] = 1.0;
    }

    this.particleVel = new Float32Array(2 * this.maxParticles);
    this.particleDensity = new Float32Array(this.fNumCells);
    this.particleRestDensity = 0.0;

    this.particleRadius = particleRadius;
    this.pInvSpacing = 1.0 / (2.2 * particleRadius);
    this.pNumX = Math.floor(width * this.pInvSpacing) + 1;
    this.pNumY = Math.floor(height * this.pInvSpacing) + 1;
    this.pNumCells = this.pNumX * this.pNumY;

    this.numCellParticles = new Int32Array(this.pNumCells);
    this.firstCellParticle = new Int32Array(this.pNumCells + 1);
    this.cellParticleIds = new Int32Array(maxParticles);

    this.numParticles = 0;
  }

  integrateParticles(dt, gravity) {
    for (let i = 0; i < this.numParticles; i++) {
      this.particleVel[2 * i + 1] += dt * gravity;
      this.particlePos[2 * i] += this.particleVel[2 * i] * dt;
      this.particlePos[2 * i + 1] += this.particleVel[2 * i + 1] * dt;
    }
  }

  pushParticlesApart(numIters) {
    const colorDiffusionCoeff = 0.001;
    const minDist = 2.0 * this.particleRadius;
    const minDist2 = minDist * minDist;

    // count particles per cell
    this.numCellParticles.fill(0);

    for (let i = 0; i < this.numParticles; i++) {
      const x = this.particlePos[2 * i];
      const y = this.particlePos[2 * i + 1];

      const xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
      const yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
      const cellNr = xi * this.pNumY + yi;
      this.numCellParticles[cellNr]++;
    }

    // partial sums
    let first = 0;
    for (let i = 0; i < this.pNumCells; i++) {
      first += this.numCellParticles[i];
      this.firstCellParticle[i] = first;
    }
    this.firstCellParticle[this.pNumCells] = first; // guard

    // fill particles into cells
    for (let i = 0; i < this.numParticles; i++) {
      const x = this.particlePos[2 * i];
      const y = this.particlePos[2 * i + 1];

      const xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
      const yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
      const cellNr = xi * this.pNumY + yi;
      this.firstCellParticle[cellNr]--;
      this.cellParticleIds[this.firstCellParticle[cellNr]] = i;
    }

    // push particles apart
    for (let iter = 0; iter < numIters; iter++) {
      for (let i = 0; i < this.numParticles; i++) {
        const px = this.particlePos[2 * i];
        const py = this.particlePos[2 * i + 1];

        const pxi = Math.floor(px * this.pInvSpacing);
        const pyi = Math.floor(py * this.pInvSpacing);
        const x0 = Math.max(pxi - 1, 0);
        const y0 = Math.max(pyi - 1, 0);
        const x1 = Math.min(pxi + 1, this.pNumX - 1);
        const y1 = Math.min(pyi + 1, this.pNumY - 1);

        for (let xi = x0; xi <= x1; xi++) {
          for (let yi = y0; yi <= y1; yi++) {
            const cellNr = xi * this.pNumY + yi;
            const first = this.firstCellParticle[cellNr];
            const last = this.firstCellParticle[cellNr + 1];

            for (let j = first; j < last; j++) {
              const id = this.cellParticleIds[j];
              if (id === i) continue;

              const qx = this.particlePos[2 * id];
              const qy = this.particlePos[2 * id + 1];

              let dx = qx - px;
              let dy = qy - py;
              const d2 = dx * dx + dy * dy;
              if (d2 > minDist2 || d2 === 0.0) continue;

              const d = Math.sqrt(d2);
              const s = (0.5 * (minDist - d)) / d;
              dx *= s;
              dy *= s;

              this.particlePos[2 * i] -= dx;
              this.particlePos[2 * i + 1] -= dy;
              this.particlePos[2 * id] += dx;
              this.particlePos[2 * id + 1] += dy;

              // diffuse colors
              for (let k = 0; k < 3; k++) {
                const color0 = this.particleColor[3 * i + k];
                const color1 = this.particleColor[3 * id + k];
                const color = (color0 + color1) * 0.5;
                this.particleColor[3 * i + k] =
                  color0 + (color - color0) * colorDiffusionCoeff;
                this.particleColor[3 * id + k] =
                  color1 + (color - color1) * colorDiffusionCoeff;
              }
            }
          }
        }
      }
    }
  }

  handleParticleCollisions(obstacleX, obstacleY, obstacleRadius) {
    const h = 1.0 / this.fInvSpacing;
    const r = this.particleRadius;
    const or = obstacleRadius;
    const or2 = or * or;
    const minDist = obstacleRadius + r;
    const minDist2 = minDist * minDist;

    const minX = h + r;
    const maxX = (this.fNumX - 1) * h - r;
    const minY = h + r;
    const maxY = (this.fNumY - 1) * h - r;

    for (let i = 0; i < this.numParticles; i++) {
      let x = this.particlePos[2 * i];
      let y = this.particlePos[2 * i + 1];

      const dx = x - obstacleX;
      const dy = y - obstacleY;
      const d2 = dx * dx + dy * dy;

      // obstacle collision
      if (d2 < minDist2) {
        this.particleVel[2 * i] = scene.obstacleVelX;
        this.particleVel[2 * i + 1] = scene.obstacleVelY;
      }

      // wall collisions
      if (x < minX) {
        x = minX;
        this.particleVel[2 * i] = 0.0;
      }
      if (x > maxX) {
        x = maxX;
        this.particleVel[2 * i] = 0.0;
      }
      if (y < minY) {
        y = minY;
        this.particleVel[2 * i + 1] = 0.0;
      }
      if (y > maxY) {
        y = maxY;
        this.particleVel[2 * i + 1] = 0.0;
      }

      this.particlePos[2 * i] = x;
      this.particlePos[2 * i + 1] = y;
    }
  }

  updateParticleDensity() {
    const n = this.fNumY;
    const h = this.h;
    const h1 = this.fInvSpacing;
    const h2 = 0.5 * h;

    const d = this.particleDensity;
    d.fill(0.0);

    for (let i = 0; i < this.numParticles; i++) {
      let x = this.particlePos[2 * i];
      let y = this.particlePos[2 * i + 1];

      x = clamp(x, h, (this.fNumX - 1) * h);
      y = clamp(y, h, (this.fNumY - 1) * h);

      const x0 = Math.floor((x - h2) * h1);
      const tx = (x - h2 - x0 * h) * h1;
      const x1 = Math.min(x0 + 1, this.fNumX - 2);

      const y0 = Math.floor((y - h2) * h1);
      const ty = (y - h2 - y0 * h) * h1;
      const y1 = Math.min(y0 + 1, this.fNumY - 2);

      const sx = 1.0 - tx;
      const sy = 1.0 - ty;

      if (x0 < this.fNumX && y0 < this.fNumY) d[x0 * n + y0] += sx * sy;
      if (x1 < this.fNumX && y0 < this.fNumY) d[x1 * n + y0] += tx * sy;
      if (x1 < this.fNumX && y1 < this.fNumY) d[x1 * n + y1] += tx * ty;
      if (x0 < this.fNumX && y1 < this.fNumY) d[x0 * n + y1] += sx * ty;
    }

    if (this.particleRestDensity === 0.0) {
      let sum = 0.0;
      let numFluidCells = 0;

      for (let i = 0; i < this.fNumCells; i++) {
        if (this.cellType[i] === FLUID_CELL) {
          sum += d[i];
          numFluidCells++;
        }
      }

      if (numFluidCells > 0) this.particleRestDensity = sum / numFluidCells;
    }
  }

  transferVelocities(toGrid, flipRatio) {
    const n = this.fNumY;
    const h = this.h;
    const h1 = this.fInvSpacing;
    const h2 = 0.5 * h;

    if (toGrid) {
      this.prevU.set(this.u);
      this.prevV.set(this.v);

      this.du.fill(0.0);
      this.dv.fill(0.0);
      this.u.fill(0.0);
      this.v.fill(0.0);

      for (let i = 0; i < this.fNumCells; i++) {
        this.cellType[i] = this.s[i] === 0.0 ? SOLID_CELL : AIR_CELL;
      }

      for (let i = 0; i < this.numParticles; i++) {
        const x = this.particlePos[2 * i];
        const y = this.particlePos[2 * i + 1];
        const xi = clamp(Math.floor(x * h1), 0, this.fNumX - 1);
        const yi = clamp(Math.floor(y * h1), 0, this.fNumY - 1);
        const cellNr = xi * n + yi;
        if (this.cellType[cellNr] === AIR_CELL) {
          this.cellType[cellNr] = FLUID_CELL;
        }
      }
    }

    for (let component = 0; component < 2; component++) {
      const dx = component === 0 ? 0.0 : h2;
      const dy = component === 0 ? h2 : 0.0;

      const f = component === 0 ? this.u : this.v;
      const prevF = component === 0 ? this.prevU : this.prevV;
      const d = component === 0 ? this.du : this.dv;

      for (let i = 0; i < this.numParticles; i++) {
        let x = this.particlePos[2 * i];
        let y = this.particlePos[2 * i + 1];

        x = clamp(x, h, (this.fNumX - 1) * h);
        y = clamp(y, h, (this.fNumY - 1) * h);

        const x0 = Math.min(Math.floor((x - dx) * h1), this.fNumX - 2);
        const tx = (x - dx - x0 * h) * h1;
        const x1 = Math.min(x0 + 1, this.fNumX - 2);

        const y0 = Math.min(Math.floor((y - dy) * h1), this.fNumY - 2);
        const ty = (y - dy - y0 * h) * h1;
        const y1 = Math.min(y0 + 1, this.fNumY - 2);

        const sx = 1.0 - tx;
        const sy = 1.0 - ty;

        const d0 = sx * sy;
        const d1 = tx * sy;
        const d2 = tx * ty;
        const d3 = sx * ty;

        const nr0 = x0 * n + y0;
        const nr1 = x1 * n + y0;
        const nr2 = x1 * n + y1;
        const nr3 = x0 * n + y1;

        if (toGrid) {
          const pv = this.particleVel[2 * i + component];
          f[nr0] += pv * d0;
          d[nr0] += d0;
          f[nr1] += pv * d1;
          d[nr1] += d1;
          f[nr2] += pv * d2;
          d[nr2] += d2;
          f[nr3] += pv * d3;
          d[nr3] += d3;
        } else {
          const offset = component === 0 ? n : 1;
          const valid0 =
            this.cellType[nr0] !== AIR_CELL ||
            this.cellType[nr0 - offset] !== AIR_CELL
              ? 1.0
              : 0.0;
          const valid1 =
            this.cellType[nr1] !== AIR_CELL ||
            this.cellType[nr1 - offset] !== AIR_CELL
              ? 1.0
              : 0.0;
          const valid2 =
            this.cellType[nr2] !== AIR_CELL ||
            this.cellType[nr2 - offset] !== AIR_CELL
              ? 1.0
              : 0.0;
          const valid3 =
            this.cellType[nr3] !== AIR_CELL ||
            this.cellType[nr3 - offset] !== AIR_CELL
              ? 1.0
              : 0.0;

          const v = this.particleVel[2 * i + component];
          const denom = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

          if (denom > 0.0) {
            const picV =
              (valid0 * d0 * f[nr0] +
                valid1 * d1 * f[nr1] +
                valid2 * d2 * f[nr2] +
                valid3 * d3 * f[nr3]) /
              denom;
            const corr =
              (valid0 * d0 * (f[nr0] - prevF[nr0]) +
                valid1 * d1 * (f[nr1] - prevF[nr1]) +
                valid2 * d2 * (f[nr2] - prevF[nr2]) +
                valid3 * d3 * (f[nr3] - prevF[nr3])) /
              denom;
            const flipV = v + corr;

            this.particleVel[2 * i + component] =
              (1.0 - flipRatio) * picV + flipRatio * flipV;
          }
        }
      }

      if (toGrid) {
        for (let i = 0; i < f.length; i++) {
          if (d[i] > 0.0) f[i] /= d[i];
        }

        // restore solid cells
        for (let i = 0; i < this.fNumX; i++) {
          for (let j = 0; j < this.fNumY; j++) {
            const solid = this.cellType[i * n + j] === SOLID_CELL;
            if (
              solid ||
              (i > 0 && this.cellType[(i - 1) * n + j] === SOLID_CELL)
            ) {
              this.u[i * n + j] = this.prevU[i * n + j];
            }
            if (
              solid ||
              (j > 0 && this.cellType[i * n + j - 1] === SOLID_CELL)
            ) {
              this.v[i * n + j] = this.prevV[i * n + j];
            }
          }
        }
      }
    }
  }

  solveIncompressibility(numIters, dt, overRelaxation, compensateDrift = true) {
    this.p.fill(0.0);
    this.prevU.set(this.u);
    this.prevV.set(this.v);

    const n = this.fNumY;
    const cp = (this.density * this.h) / dt;

    for (let iter = 0; iter < numIters; iter++) {
      for (let i = 1; i < this.fNumX - 1; i++) {
        for (let j = 1; j < this.fNumY - 1; j++) {
          if (this.cellType[i * n + j] !== FLUID_CELL) continue;

          const center = i * n + j;
          const left = (i - 1) * n + j;
          const right = (i + 1) * n + j;
          const bottom = i * n + j - 1;
          const top = i * n + j + 1;

          const s = this.s[center];
          const sx0 = this.s[left];
          const sx1 = this.s[right];
          const sy0 = this.s[bottom];
          const sy1 = this.s[top];
          const sumS = sx0 + sx1 + sy0 + sy1;
          if (sumS === 0.0) continue;

          let div =
            this.u[right] - this.u[center] + this.v[top] - this.v[center];

          if (this.particleRestDensity > 0.0 && compensateDrift) {
            const k = 1.0;
            const compression =
              this.particleDensity[center] - this.particleRestDensity;
            if (compression > 0.0) div = div - k * compression;
          }

          const p = (-div / sumS) * overRelaxation;
          this.p[center] += cp * p;

          this.u[center] -= sx0 * p;
          this.u[right] += sx1 * p;
          this.v[center] -= sy0 * p;
          this.v[top] += sy1 * p;
        }
      }
    }
  }

  updateParticleColors() {
    // const h1 = this.fInvSpacing;
    // for (let i = 0; i < this.numParticles; i++) {
    //   const s = 0.01;
    //   this.particleColor[3 * i] = clamp(
    //     this.particleColor[3 * i] - s,
    //     0.0,
    //     1.0
    //   );
    //   this.particleColor[3 * i + 1] = clamp(
    //     this.particleColor[3 * i + 1] - s,
    //     0.0,
    //     1.0
    //   );
    //   this.particleColor[3 * i + 2] = clamp(
    //     this.particleColor[3 * i + 2] + s,
    //     0.0,
    //     1.0
    //   );
    //   const x = this.particlePos[2 * i];
    //   const y = this.particlePos[2 * i + 1];
    //   const xi = clamp(Math.floor(x * h1), 1, this.fNumX - 1);
    //   const yi = clamp(Math.floor(y * h1), 1, this.fNumY - 1);
    //   const cellNr = xi * this.fNumY + yi;
    //   const d0 = this.particleRestDensity;
    //   if (d0 > 0.0) {
    //     const relDensity = this.particleDensity[cellNr] / d0;
    //     if (relDensity < 0.7) {
    //       const s = 0.8;
    //       this.particleColor[3 * i] = s;
    //       this.particleColor[3 * i + 1] = s;
    //       this.particleColor[3 * i + 2] = 1.0;
    //     }
    //   }
    // }
  }

  setSciColor(cellNr, val, minVal, maxVal) {
    val = Math.min(Math.max(val, minVal), maxVal - 0.0001);
    const d = maxVal - minVal;
    val = d === 0.0 ? 0.5 : (val - minVal) / d;
    const m = 0.25;
    const num = Math.floor(val / m);
    const s = (val - num * m) / m;
    let r, g, b;

    // r = 1.0;
    // g = 1.0;
    // b = 1.0;

    switch (num) {
      case 0:
        r = 0.0;
        g = s;
        b = 1.0;
        break;
      case 1:
        r = 0.0;
        g = 1.0;
        b = 1.0 - s;
        break;
      case 2:
        r = s;
        g = 1.0;
        b = 0.0;
        break;
      case 3:
        r = 1.0;
        g = 1.0 - s;
        b = 0.0;
        break;
    }

    this.cellColor[3 * cellNr] = r;
    this.cellColor[3 * cellNr + 1] = g;
    this.cellColor[3 * cellNr + 2] = b;
  }

  updateCellColors() {
    this.cellColor.fill(0.0);

    for (let i = 0; i < this.fNumCells; i++) {
      if (this.cellType[i] === SOLID_CELL) {
        this.cellColor[3 * i] = 0.5;
        this.cellColor[3 * i + 1] = 0.5;
        this.cellColor[3 * i + 2] = 0.5;
      } else if (this.cellType[i] === FLUID_CELL) {
        let d = this.particleDensity[i];
        if (this.particleRestDensity > 0.0) d /= this.particleRestDensity;
        this.setSciColor(i, d, 0.0, 2.0);
      }
    }
  }

  simulate(
    dt,
    gravity,
    flipRatio,
    numPressureIters,
    numParticleIters,
    overRelaxation,
    compensateDrift,
    separateParticles,
    obstacleX,
    obstacleY,
    obstacleRadius,
  ) {
    const numSubSteps = 1;
    const sdt = dt / numSubSteps;

    for (let step = 0; step < numSubSteps; step++) {
      this.integrateParticles(sdt, gravity);
      if (separateParticles) this.pushParticlesApart(numParticleIters);
      this.handleParticleCollisions(obstacleX, obstacleY, obstacleRadius);
      this.transferVelocities(true);
      this.updateParticleDensity();
      this.solveIncompressibility(
        numPressureIters,
        sdt,
        overRelaxation,
        compensateDrift,
      );
      this.transferVelocities(false, flipRatio);
    }

    this.updateParticleColors();
    this.updateCellColors();
  }
}

// ----------------- end of simulator ------------------------------

const scene = {
  gravity: -9.81,
  dt: 1.0 / 120.0,
  flipRatio: 0.9,
  numPressureIters: 100,
  numParticleIters: 2,
  frameNr: 0,
  overRelaxation: 1.9,
  compensateDrift: true,
  separateParticles: true,
  obstacleX: 0.0,
  obstacleY: 0.0,
  // obstacleRadius: 0.15,
  obstacleRadius: 0.1,
  paused: false,
  showObstacle: true,
  obstacleVelX: 0.0,
  obstacleVelY: 0.0,
  showParticles: true,
  showGrid: false,
  fluid: null,
};

function setupScene() {
  scene.obstacleRadius = 0.1;
  scene.overRelaxation = 1.9;
  scene.dt = 1.0 / 60.0;
  scene.numPressureIters = 50;
  scene.numParticleIters = 2;

  const res = 100;
  const tankHeight = 1.0 * simHeight;
  const tankWidth = 1.0 * simWidth;
  const h = (tankHeight / res) * 3.5; // Changed by me
  const density = 1000.0;
  const relWaterHeight = 0.9; // Change to change overall water height
  const relWaterWidth = 0.6;

  // compute number of particles
  const r = 0.35 * h; // particle radius w.r.t. cell size // Changed by me
  const dx = 2.0 * r;
  const dy = (Math.sqrt(3.0) / 2.0) * dx;

  const numX = Math.floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
  const numY = Math.floor(
    (relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy,
  );
  const maxParticles = numX * numY;
  // const maxParticles = 28000;

  // create fluid
  const f = (scene.fluid = new FlipFluid(
    density,
    tankWidth,
    tankHeight,
    h,
    r,
    maxParticles,
  ));

  // create particles
  f.numParticles = numX * numY;
  let p = 0;
  for (let i = 0; i < numX; i++) {
    for (let j = 0; j < numY; j++) {
      f.particlePos[p++] = h + r + dx * i + (j % 2 === 0 ? 0.0 : r);
      f.particlePos[p++] = h + r + dy * j;
    }
  }

  // setup grid cells for tank
  const n = f.fNumY;
  for (let i = 0; i < f.fNumX; i++) {
    for (let j = 0; j < f.fNumY; j++) {
      f.s[i * n + j] = i === 0 || i === f.fNumX - 1 || j === 0 ? 0.0 : 1.0;
    }
  }

  setObstacle(3.0, 2.0, true);
}

// draw -------------------------------------------------------

const pointVertexShader = `
  attribute vec2 attrPosition;
  attribute vec3 attrColor;
  uniform vec2 domainSize;
  uniform float pointSize;
  uniform float drawDisk;

  varying vec3 fragColor;
  varying float fragDrawDisk;

  void main() {
    vec4 screenTransform = 
      vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);
    gl_Position =
      vec4(attrPosition * screenTransform.xy + screenTransform.zw, 0.0, 1.0);

    gl_PointSize = pointSize;
    fragColor = attrColor;
    fragDrawDisk = drawDisk;
  }
`;

const pointFragmentShader = `
  precision mediump float;
  varying vec3 fragColor;
  varying float fragDrawDisk;

  void main() {
    if (fragDrawDisk == 1.0) {
      float rx = 0.5 - gl_PointCoord.x;
      float ry = 0.5 - gl_PointCoord.y;
      float r2 = rx * rx + ry * ry;
      if (r2 > 0.25)
        discard;
    }
    gl_FragColor = vec4(fragColor, 1.0);
  }
`;

const meshVertexShader = `
  attribute vec2 attrPosition;
  uniform vec2 domainSize;
  uniform vec3 color;
  uniform vec2 translation;
  uniform float scale;

  varying vec3 fragColor;

  void main() {
    vec2 v = translation + attrPosition * scale;
    vec4 screenTransform = 
      vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);
    gl_Position =
      vec4(v * screenTransform.xy + screenTransform.zw, 0.0, 1.0);

    fragColor = color;
  }
`;

const meshFragmentShader = `
  precision mediump float;
  varying vec3 fragColor;

  void main() {
    gl_FragColor = vec4(fragColor, 0.0);
  }
`;

const createShader = (gl, vsSource, fsSource) => {
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

  const vsShader = gl.createShader(gl.VERTEX_SHADER);
  gl.shaderSource(vsShader, vsSource);
  gl.compileShader(vsShader);
  if (!gl.getShaderParameter(vsShader, gl.COMPILE_STATUS)) {
    console.log(
      'vertex shader compile error: ' + gl.getShaderInfoLog(vsShader),
    );
  }

  const fsShader = gl.createShader(gl.FRAGMENT_SHADER);
  gl.shaderSource(fsShader, fsSource);
  gl.compileShader(fsShader);
  if (!gl.getShaderParameter(fsShader, gl.COMPILE_STATUS)) {
    console.log(
      'fragment shader compile error: ' + gl.getShaderInfoLog(fsShader),
    );
  }

  const shader = gl.createProgram();
  gl.attachShader(shader, vsShader);
  gl.attachShader(shader, fsShader);
  gl.linkProgram(shader);

  return shader;
};

let pointShader = null;
let meshShader = null;
let pointVertexBuffer = null;
let pointColorBuffer = null;
let gridVertBuffer = null;
let gridColorBuffer = null;
let diskVertBuffer = null;
let diskIdBuffer = null;

function draw() {
  gl.clearColor(0.0, 0.0, 0.0, 1.0);
  gl.clear(gl.COLOR_BUFFER_BIT);
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

  // prepare shaders
  if (pointShader === null) {
    pointShader = createShader(gl, pointVertexShader, pointFragmentShader);
  }
  if (meshShader === null) {
    meshShader = createShader(gl, meshVertexShader, meshFragmentShader);
  }

  // grid
  if (gridVertBuffer === null) {
    const f = scene.fluid;
    gridVertBuffer = gl.createBuffer();
    const cellCenters = new Float32Array(2 * f.fNumCells);
    let p = 0;

    for (let i = 0; i < f.fNumX; i++) {
      for (let j = 0; j < f.fNumY; j++) {
        cellCenters[p++] = (i + 0.5) * f.h;
        cellCenters[p++] = (j + 0.5) * f.h;
      }
    }
    gl.bindBuffer(gl.ARRAY_BUFFER, gridVertBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, cellCenters, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
  }

  if (gridColorBuffer === null) {
    gridColorBuffer = gl.createBuffer();
  }

  if (scene.showGrid) {
    const pointSize = ((0.9 * scene.fluid.h) / simWidth) * canvas.width;

    gl.useProgram(pointShader);
    gl.uniform2f(
      gl.getUniformLocation(pointShader, 'domainSize'),
      simWidth,
      simHeight,
    );
    gl.uniform1f(gl.getUniformLocation(pointShader, 'pointSize'), pointSize);
    gl.uniform1f(gl.getUniformLocation(pointShader, 'drawDisk'), 0.0);

    gl.bindBuffer(gl.ARRAY_BUFFER, gridVertBuffer);
    const posLoc = gl.getAttribLocation(pointShader, 'attrPosition');
    gl.enableVertexAttribArray(posLoc);
    gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, gridColorBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.cellColor, gl.DYNAMIC_DRAW);

    const colorLoc = gl.getAttribLocation(pointShader, 'attrColor');
    gl.enableVertexAttribArray(colorLoc);
    gl.vertexAttribPointer(colorLoc, 3, gl.FLOAT, false, 0, 0);

    gl.drawArrays(gl.POINTS, 0, scene.fluid.fNumCells);

    gl.disableVertexAttribArray(posLoc);
    gl.disableVertexAttribArray(colorLoc);
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
  }

  // water
  if (scene.showParticles) {
    gl.clear(gl.DEPTH_BUFFER_BIT);
    const pointSize =
      ((2.0 * scene.fluid.particleRadius) / simWidth) * canvas.width;

    gl.useProgram(pointShader);
    gl.uniform2f(
      gl.getUniformLocation(pointShader, 'domainSize'),
      simWidth,
      simHeight,
    );
    gl.uniform1f(gl.getUniformLocation(pointShader, 'pointSize'), pointSize);
    gl.uniform1f(gl.getUniformLocation(pointShader, 'drawDisk'), 1.0);

    if (pointVertexBuffer === null) pointVertexBuffer = gl.createBuffer();
    if (pointColorBuffer === null) pointColorBuffer = gl.createBuffer();

    gl.bindBuffer(gl.ARRAY_BUFFER, pointVertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.particlePos, gl.DYNAMIC_DRAW);

    const posLoc = gl.getAttribLocation(pointShader, 'attrPosition');
    gl.enableVertexAttribArray(posLoc);
    gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, pointColorBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.particleColor, gl.DYNAMIC_DRAW);

    const colorLoc = gl.getAttribLocation(pointShader, 'attrColor');
    gl.enableVertexAttribArray(colorLoc);
    gl.vertexAttribPointer(colorLoc, 3, gl.FLOAT, false, 0, 0);

    gl.drawArrays(gl.POINTS, 0, scene.fluid.numParticles);

    gl.disableVertexAttribArray(posLoc);
    gl.disableVertexAttribArray(colorLoc);
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
  }

  // disk
  const numSegs = 50;

  if (diskVertBuffer === null) {
    diskVertBuffer = gl.createBuffer();
    const dphi = (2.0 * Math.PI) / numSegs;
    const diskVerts = new Float32Array(2 * numSegs + 2);
    let p = 0;
    diskVerts[p++] = 0.0;
    diskVerts[p++] = 0.0;
    for (let i = 0; i < numSegs; i++) {
      diskVerts[p++] = Math.cos(i * dphi);
      diskVerts[p++] = Math.sin(i * dphi);
    }
    gl.bindBuffer(gl.ARRAY_BUFFER, diskVertBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, diskVerts, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, null);

    diskIdBuffer = gl.createBuffer();
    const diskIds = new Uint16Array(3 * numSegs);
    p = 0;
    for (let i = 0; i < numSegs; i++) {
      diskIds[p++] = 0;
      diskIds[p++] = 1 + i;
      diskIds[p++] = 1 + ((i + 1) % numSegs);
    }

    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, diskIdBuffer);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, diskIds, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
  }

  gl.clear(gl.DEPTH_BUFFER_BIT);
  const diskColor = [1.0, 1.0, 1.0];

  gl.useProgram(meshShader);
  gl.uniform2f(
    gl.getUniformLocation(meshShader, 'domainSize'),
    simWidth,
    simHeight,
  );
  gl.uniform3f(
    gl.getUniformLocation(meshShader, 'color'),
    diskColor[0],
    diskColor[1],
    diskColor[2],
  );
  gl.uniform2f(
    gl.getUniformLocation(meshShader, 'translation'),
    scene.obstacleX,
    scene.obstacleY,
  );
  gl.uniform1f(
    gl.getUniformLocation(meshShader, 'scale'),
    scene.obstacleRadius + scene.fluid.particleRadius,
  );

  const posLoc = gl.getAttribLocation(meshShader, 'attrPosition');
  gl.enableVertexAttribArray(posLoc);
  gl.bindBuffer(gl.ARRAY_BUFFER, diskVertBuffer);
  gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, diskIdBuffer);
  gl.drawElements(gl.TRIANGLES, 3 * numSegs, gl.UNSIGNED_SHORT, 0);

  gl.disableVertexAttribArray(posLoc);
}

function setObstacle(x, y, reset) {
  let vx = 0.0;
  let vy = 0.0;

  if (!reset) {
    vx = (x - scene.obstacleX) / scene.dt;
    vy = (y - scene.obstacleY) / scene.dt;
  }

  scene.obstacleX = x;
  scene.obstacleY = y;
  const r = scene.obstacleRadius;
  const f = scene.fluid;
  const n = f.fNumY;
  const cd = Math.sqrt(2) * f.h;

  for (let i = 1; i < f.fNumX - 2; i++) {
    for (let j = 1; j < f.fNumY - 2; j++) {
      f.s[i * n + j] = 1.0;

      const dx = (i + 0.5) * f.h - x;
      const dy = (j + 0.5) * f.h - y;

      if (dx * dx + dy * dy < r * r) {
        f.s[i * n + j] = 0.0;
        f.u[i * n + j] = vx;
        f.u[(i + 1) * n + j] = vx;
        f.v[i * n + j] = vy;
        f.v[i * n + j + 1] = vy;
      }
    }
  }

  scene.showObstacle = true;
  scene.obstacleVelX = vx;
  scene.obstacleVelY = vy;
}

// interaction -------------------------------------------------------

let mouseDown = true;

const startDrag = (x, y) => {
  const bounds = canvas.getBoundingClientRect();
  const mx = x - bounds.left - canvas.clientLeft;
  const my = y - bounds.top - canvas.clientTop;
  mouseDown = true;

  const scaledX = mx / cScale;
  const scaledY = (canvas.height - my) / cScale;

  setObstacle(scaledX, scaledY, true);
  scene.paused = false;
};

const drag = (x, y) => {
  if (mouseDown) {
    const bounds = canvas.getBoundingClientRect();
    const mx = x - bounds.left - canvas.clientLeft;
    const my = y - bounds.top - canvas.clientTop;
    const scaledX = mx / cScale;
    const scaledY = (canvas.height - my) / cScale;
    setObstacle(scaledX, scaledY, false);
  }
};

const endDrag = () => {
  // mouseDown = false;
  scene.obstacleVelX = 0.0;
  scene.obstacleVelY = 0.0;
};

canvas.addEventListener('mousedown', (event) => {
  startDrag(event.x, event.y);
});

canvas.addEventListener('mouseup', () => {
  endDrag();
});

canvas.addEventListener('mousemove', (event) => {
  drag(event.x, event.y);
});

canvas.addEventListener('touchstart', (event) => {
  startDrag(event.touches[0].clientX, event.touches[0].clientY);
});

canvas.addEventListener('touchend', () => {
  endDrag();
});

canvas.addEventListener(
  'touchmove',
  (event) => {
    event.preventDefault();
    event.stopImmediatePropagation();
    drag(event.touches[0].clientX, event.touches[0].clientY);
  },
  { passive: false },
);

document.addEventListener('keydown', (event) => {
  switch (event.key) {
    case 'p':
      scene.paused = !scene.paused;
      break;
    case 'm':
      scene.paused = false;
      simulate();
      scene.paused = true;
      break;
  }
});

const toggleStart = () => {
  const button = document.getElementById('startButton');
  button.innerHTML = scene.paused ? 'Stop' : 'Start';
  scene.paused = !scene.paused;
};

// main -------------------------------------------------------

const simulate = () => {
  if (!scene.paused) {
    scene.fluid.simulate(
      scene.dt,
      scene.gravity,
      scene.flipRatio,
      scene.numPressureIters,
      scene.numParticleIters,
      scene.overRelaxation,
      scene.compensateDrift,
      scene.separateParticles,
      scene.obstacleX,
      scene.obstacleY,
      scene.obstacleRadius,
    );
  }
  scene.frameNr++;
};

const update = () => {
  simulate();
  draw();
  requestAnimationFrame(update);
};

setupScene();
update();
