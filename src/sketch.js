// FLIP Fluid Simulation in p5.js
// Adapted from WebGL implementation

// Global variables
let simWidth, simHeight;
let cScale;
let f; // fluid simulator
let scene;
let mouseDown = false;

// Constants
const U_FIELD = 0;
const V_FIELD = 1;
const FLUID_CELL = 0;
const AIR_CELL = 1;
const SOLID_CELL = 2;

function clamp(x, min, max) {
  return Math.max(min, Math.min(max, x));
}

// ----------------- Fluid Simulator ------------------------------
class FlipFluid {
  constructor(density, width, height, spacing, particleRadius, maxParticles) {
    // Fluid grid properties
    this.density = density;
    this.fNumX = Math.floor(width / spacing) + 1;
    this.fNumY = Math.floor(height / spacing) + 1;
    this.h = Math.max(width / this.fNumX, height / this.fNumY);
    this.fInvSpacing = 1.0 / this.h;
    this.fNumCells = this.fNumX * this.fNumY;

    // Grid velocity fields
    this.u = new Float32Array(this.fNumCells);
    this.v = new Float32Array(this.fNumCells);
    this.du = new Float32Array(this.fNumCells);
    this.dv = new Float32Array(this.fNumCells);
    this.prevU = new Float32Array(this.fNumCells);
    this.prevV = new Float32Array(this.fNumCells);

    // Pressure and solid fields
    this.p = new Float32Array(this.fNumCells);
    this.s = new Float32Array(this.fNumCells);
    this.cellType = new Int32Array(this.fNumCells);
    this.cellColor = new Float32Array(3 * this.fNumCells);

    // Particle properties
    this.maxParticles = maxParticles;
    this.particlePos = new Float32Array(2 * this.maxParticles);
    this.particleColor = new Float32Array(3 * this.maxParticles);

    // Initialize particle colors to blue
    for (let i = 0; i < this.maxParticles; i++) {
      this.particleColor[3 * i + 2] = 1.0;
    }

    this.particleVel = new Float32Array(2 * this.maxParticles);
    this.particleDensity = new Float32Array(this.fNumCells);
    this.particleRestDensity = 0.0;

    // Particle spatial indexing
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

    // Count particles per cell
    this.numCellParticles.fill(0);

    for (let i = 0; i < this.numParticles; i++) {
      const x = this.particlePos[2 * i];
      const y = this.particlePos[2 * i + 1];

      const xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
      const yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
      const cellNr = xi * this.pNumY + yi;
      this.numCellParticles[cellNr]++;
    }

    // Calculate partial sums
    let first = 0;
    for (let i = 0; i < this.pNumCells; i++) {
      first += this.numCellParticles[i];
      this.firstCellParticle[i] = first;
    }
    this.firstCellParticle[this.pNumCells] = first; // guard

    // Fill particles into cells
    for (let i = 0; i < this.numParticles; i++) {
      const x = this.particlePos[2 * i];
      const y = this.particlePos[2 * i + 1];

      const xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
      const yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
      const cellNr = xi * this.pNumY + yi;
      this.firstCellParticle[cellNr]--;
      this.cellParticleIds[this.firstCellParticle[cellNr]] = i;
    }

    // Push particles apart
    const minDist = 2.0 * this.particleRadius;
    const minDist2 = minDist * minDist;

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
              if (id == i) continue;

              const qx = this.particlePos[2 * id];
              const qy = this.particlePos[2 * id + 1];

              const dx = qx - px;
              const dy = qy - py;
              const d2 = dx * dx + dy * dy;

              if (d2 > minDist2 || d2 == 0.0) continue;

              const d = Math.sqrt(d2);
              const s = (0.5 * (minDist - d)) / d;
              const moveX = dx * s;
              const moveY = dy * s;

              this.particlePos[2 * i] -= moveX;
              this.particlePos[2 * i + 1] -= moveY;
              this.particlePos[2 * id] += moveX;
              this.particlePos[2 * id + 1] += moveY;

              // Diffuse colors
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
    const minDist = obstacleRadius + r;
    const minDist2 = minDist * minDist;

    const minX = h + r;
    const maxX = (this.fNumX - 1) * h - r;
    const minY = h + r;
    const maxY = (this.fNumY - 1) * h - r;

    for (let i = 0; i < this.numParticles; i++) {
      let x = this.particlePos[2 * i];
      let y = this.particlePos[2 * i + 1];

      // Obstacle collision
      const dx = x - obstacleX;
      const dy = y - obstacleY;
      const d2 = dx * dx + dy * dy;

      if (d2 < minDist2) {
        // Transfer obstacle velocity to particles
        this.particleVel[2 * i] = scene.obstacleVelX;
        this.particleVel[2 * i + 1] = scene.obstacleVelY;
      }

      // Wall collisions
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

    // Calculate rest density if not set
    if (this.particleRestDensity == 0.0) {
      let sum = 0.0;
      let numFluidCells = 0;

      for (let i = 0; i < this.fNumCells; i++) {
        if (this.cellType[i] == FLUID_CELL) {
          sum += d[i];
          numFluidCells++;
        }
      }

      if (numFluidCells > 0) this.particleRestDensity = sum / numFluidCells;
    }
  }

  transferVelocities(toGrid, flipRatio = 0.9) {
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

      // Mark cells as solid or air initially
      for (let i = 0; i < this.fNumCells; i++) {
        this.cellType[i] = this.s[i] == 0.0 ? SOLID_CELL : AIR_CELL;
      }

      // Mark fluid cells based on particle positions
      for (let i = 0; i < this.numParticles; i++) {
        const x = this.particlePos[2 * i];
        const y = this.particlePos[2 * i + 1];
        const xi = clamp(Math.floor(x * h1), 0, this.fNumX - 1);
        const yi = clamp(Math.floor(y * h1), 0, this.fNumY - 1);
        const cellNr = xi * n + yi;
        if (this.cellType[cellNr] == AIR_CELL) {
          this.cellType[cellNr] = FLUID_CELL;
        }
      }
    }

    // Process u and v components
    for (let component = 0; component < 2; component++) {
      const dx = component == 0 ? 0.0 : h2;
      const dy = component == 0 ? h2 : 0.0;

      const f = component == 0 ? this.u : this.v;
      const prevF = component == 0 ? this.prevU : this.prevV;
      const d = component == 0 ? this.du : this.dv;

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
          // Transfer particle velocities to grid
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
          // Update particle velocities from grid
          const offset = component == 0 ? n : 1;
          const valid0 =
            this.cellType[nr0] != AIR_CELL ||
            this.cellType[nr0 - offset] != AIR_CELL
              ? 1.0
              : 0.0;
          const valid1 =
            this.cellType[nr1] != AIR_CELL ||
            this.cellType[nr1 - offset] != AIR_CELL
              ? 1.0
              : 0.0;
          const valid2 =
            this.cellType[nr2] != AIR_CELL ||
            this.cellType[nr2 - offset] != AIR_CELL
              ? 1.0
              : 0.0;
          const valid3 =
            this.cellType[nr3] != AIR_CELL ||
            this.cellType[nr3 - offset] != AIR_CELL
              ? 1.0
              : 0.0;

          const v = this.particleVel[2 * i + component];
          const weight = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

          if (weight > 0.0) {
            // PIC (Particle-in-Cell) velocity
            const picV =
              (valid0 * d0 * f[nr0] +
                valid1 * d1 * f[nr1] +
                valid2 * d2 * f[nr2] +
                valid3 * d3 * f[nr3]) /
              weight;

            // FLIP (Fluid-Implicit-Particle) correction
            const corr =
              (valid0 * d0 * (f[nr0] - prevF[nr0]) +
                valid1 * d1 * (f[nr1] - prevF[nr1]) +
                valid2 * d2 * (f[nr2] - prevF[nr2]) +
                valid3 * d3 * (f[nr3] - prevF[nr3])) /
              weight;

            const flipV = v + corr;

            // Use combination of PIC and FLIP
            this.particleVel[2 * i + component] =
              (1.0 - flipRatio) * picV + flipRatio * flipV;
          }
        }
      }

      if (toGrid) {
        // Normalize velocities
        for (let i = 0; i < f.length; i++) {
          if (d[i] > 0.0) f[i] /= d[i];
        }

        // Restore solid cell velocities
        for (let i = 0; i < this.fNumX; i++) {
          for (let j = 0; j < this.fNumY; j++) {
            const solid = this.cellType[i * n + j] == SOLID_CELL;
            if (
              solid ||
              (i > 0 && this.cellType[(i - 1) * n + j] == SOLID_CELL)
            ) {
              this.u[i * n + j] = this.prevU[i * n + j];
            }
            if (
              solid ||
              (j > 0 && this.cellType[i * n + j - 1] == SOLID_CELL)
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
          if (this.cellType[i * n + j] != FLUID_CELL) continue;

          const center = i * n + j;
          const left = (i - 1) * n + j;
          const right = (i + 1) * n + j;
          const bottom = i * n + j - 1;
          const top = i * n + j + 1;

          const sx0 = this.s[left];
          const sx1 = this.s[right];
          const sy0 = this.s[bottom];
          const sy1 = this.s[top];
          const s = sx0 + sx1 + sy0 + sy1;

          if (s == 0.0) continue;

          let div =
            this.u[right] - this.u[center] + this.v[top] - this.v[center];

          // Apply density correction if needed
          if (this.particleRestDensity > 0.0 && compensateDrift) {
            const k = 1.0;
            const compression =
              this.particleDensity[center] - this.particleRestDensity;
            if (compression > 0.0) div = div - k * compression;
          }

          let p = -div / s;
          p *= overRelaxation;
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
    const h1 = this.fInvSpacing;

    for (let i = 0; i < this.numParticles; i++) {
      const s = 0.01;

      // Gradually shift colors to blue
      this.particleColor[3 * i] = clamp(
        this.particleColor[3 * i] - s,
        0.0,
        1.0,
      );
      this.particleColor[3 * i + 1] = clamp(
        this.particleColor[3 * i + 1] - s,
        0.0,
        1.0,
      );
      this.particleColor[3 * i + 2] = clamp(
        this.particleColor[3 * i + 2] + s,
        0.0,
        1.0,
      );

      // Change color based on density
      const x = this.particlePos[2 * i];
      const y = this.particlePos[2 * i + 1];
      const xi = clamp(Math.floor(x * h1), 1, this.fNumX - 1);
      const yi = clamp(Math.floor(y * h1), 1, this.fNumY - 1);
      const cellNr = xi * this.fNumY + yi;

      const d0 = this.particleRestDensity;

      if (d0 > 0.0) {
        const relDensity = this.particleDensity[cellNr] / d0;
        if (relDensity < 0.7) {
          // Low density areas become lighter blue
          const s = 0.8;
          this.particleColor[3 * i] = s;
          this.particleColor[3 * i + 1] = s;
          this.particleColor[3 * i + 2] = 1.0;
        }
      }
    }
  }

  setSciColor(cellNr, val, minVal, maxVal) {
    val = Math.min(Math.max(val, minVal), maxVal - 0.0001);
    const d = maxVal - minVal;
    val = d == 0.0 ? 0.5 : (val - minVal) / d;
    const m = 0.25;
    const num = Math.floor(val / m);
    const s = (val - num * m) / m;
    let r, g, b;

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
      if (this.cellType[i] == SOLID_CELL) {
        this.cellColor[3 * i] = 0.5;
        this.cellColor[3 * i + 1] = 0.5;
        this.cellColor[3 * i + 2] = 0.5;
      } else if (this.cellType[i] == FLUID_CELL) {
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

// Setup scene with initial conditions
function setupScene() {
  scene.obstacleRadius = 0.15;
  scene.overRelaxation = 1.9;
  scene.dt = 1.0 / 60.0;
  scene.numPressureIters = 50;
  scene.numParticleIters = 2;

  const res = 100;
  const tankHeight = 1.0 * simHeight;
  const tankWidth = 1.0 * simWidth;
  const h = (tankHeight / res) * 3;
  const density = 1000.0;

  const relWaterHeight = 0.8;
  const relWaterWidth = 0.6;

  // Calculate particles
  const r = 0.35 * h; // particle radius w.r.t. cell size
  const dx = 2.0 * r;
  const dy = (Math.sqrt(3.0) / 2.0) * dx;

  const numX = Math.floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
  const numY = Math.floor(
    (relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy,
  );
  const maxParticles = numX * numY;

  // Create fluid system
  f = scene.fluid = new FlipFluid(
    density,
    tankWidth,
    tankHeight,
    h,
    r,
    maxParticles,
  );

  // Create initial particles
  f.numParticles = numX * numY;
  let p = 0;
  for (let i = 0; i < numX; i++) {
    for (let j = 0; j < numY; j++) {
      f.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
      f.particlePos[p++] = h + r + dy * j;
    }
  }

  // Setup grid cells for tank
  const n = f.fNumY;
  for (let i = 0; i < f.fNumX; i++) {
    for (let j = 0; j < f.fNumY; j++) {
      let s = 1.0; // fluid
      if (i == 0 || i == f.fNumX - 1 || j == 0) s = 0.0; // solid
      f.s[i * n + j] = s;
    }
  }

  setObstacle(3.0, 2.0, true);
}

// Set obstacle position and velocity
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

// p5.js setup and draw functions
function setup() {
  createCanvas(windowWidth - 20, windowHeight - 20);

  // Calculate simulation dimensions
  simHeight = 3.0;
  cScale = height / simHeight;
  simWidth = width / cScale;

  // Initialize scene object
  scene = {
    gravity: -9.81,
    dt: 1.0 / 60.0,
    flipRatio: 0.9,
    numPressureIters: 50,
    numParticleIters: 2,
    frameNr: 0,
    overRelaxation: 1.9,
    compensateDrift: true,
    separateParticles: true,
    obstacleX: 0.0,
    obstacleY: 0.0,
    obstacleRadius: 0.15,
    paused: true,
    showObstacle: true,
    obstacleVelX: 0.0,
    obstacleVelY: 0.0,
    showParticles: true,
    showGrid: false,
    fluid: null,
  };

  setupScene();

  // Create UI elements
  let button = createButton('Start/Stop');
  button.position(10, 10);
  button.mousePressed(togglePause);

  let gridToggle = createButton('Toggle Grid');
  gridToggle.position(100, 10);
  gridToggle.mousePressed(() => (scene.showGrid = !scene.showGrid));
}

function draw() {
  // Simulate physics if not paused
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
    scene.frameNr++;
  }

  // Render the simulation
  drawSimulation();
}

function drawSimulation() {
  background(0);

  // Draw grid cells (optional)
  if (scene.showGrid) {
    const f = scene.fluid;
    const cellSize = f.h * cScale;

    for (let i = 0; i < f.fNumX; i++) {
      for (let j = 0; j < f.fNumY; j++) {
        const cellNr = i * f.fNumY + j;
        const cellType = f.cellType[cellNr];

        if (cellType === FLUID_CELL || cellType === SOLID_CELL) {
          const r = f.cellColor[3 * cellNr] * 255;
          const g = f.cellColor[3 * cellNr + 1] * 255;
          const b = f.cellColor[3 * cellNr + 2] * 255;

          fill(r, g, b);
          noStroke();
          rect(i * cellSize, height - (j + 1) * cellSize, cellSize, cellSize);
        }
      }
    }
  }

  // Draw particles
  if (scene.showParticles) {
    const f = scene.fluid;
    const particleRadius = f.particleRadius * cScale;

    noStroke();
    for (let i = 0; i < f.numParticles; i++) {
      const x = f.particlePos[2 * i] * cScale;
      const y = height - f.particlePos[2 * i + 1] * cScale;

      const r = f.particleColor[3 * i] * 255;
      const g = f.particleColor[3 * i + 1] * 255;
      const b = f.particleColor[3 * i + 2] * 255;

      fill(r, g, b);
      circle(x, y, particleRadius * 2);
    }
  }

  // Draw obstacle
  if (scene.showObstacle) {
    const obstacleX = scene.obstacleX * cScale;
    const obstacleY = height - scene.obstacleY * cScale;
    const obstacleRadius = scene.obstacleRadius * cScale;

    fill(255, 0, 0);
    noStroke();
    circle(obstacleX, obstacleY, obstacleRadius * 2);
  }
}

// Toggle simulation pause state
function togglePause() {
  scene.paused = !scene.paused;
}

// Mouse and touch interaction functions
function mouseDragged() {
  if (mouseIsPressed) {
    let x = mouseX / cScale;
    let y = (height - mouseY) / cScale;
    setObstacle(x, y, false);
  }
  return false; // Prevent default
}

function mousePressed() {
  let x = mouseX / cScale;
  let y = (height - mouseY) / cScale;
  setObstacle(x, y, true);
  scene.paused = false;
  return false; // Prevent default
}

function mouseReleased() {
  scene.obstacleVelX = 0.0;
  scene.obstacleVelY = 0.0;
  return false; // Prevent default
}

// Handle window resize
function windowResized() {
  resizeCanvas(windowWidth - 20, windowHeight - 20);

  // Recalculate simulation dimensions
  simHeight = 3.0;
  cScale = height / simHeight;
  simWidth = width / cScale;
}
