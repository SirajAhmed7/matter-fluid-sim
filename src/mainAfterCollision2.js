import './style.css';

const wrapper = document.querySelector('.wrapper');
const canvas = document.getElementById('canvas');
const gl = canvas.getContext('webgl');
// const gl2 = canvas.getContext('2d');
// canvas.width = window.innerWidth - 20;
// canvas.height = window.innerHeight - 20;
// canvas.width = wrapper.clientWidth + wrapper.clientWidth * 0.055;
// canvas.height = wrapper.clientHeight + wrapper.clientHeight * 0.041;
canvas.width = wrapper.clientWidth + wrapper.clientWidth * 0.2;
canvas.height = wrapper.clientHeight + wrapper.clientHeight * 0.27;

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

let logoTexture = null;

function loadTexture(gl, url) {
  const texture = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, texture);

  // Fill with a temporary color while the image loads
  gl.texImage2D(
    gl.TEXTURE_2D,
    0,
    gl.RGBA,
    1,
    1,
    0,
    gl.RGBA,
    gl.UNSIGNED_BYTE,
    new Uint8Array([255, 255, 255, 255]),
  );

  const image = new Image();
  image.onload = function () {
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);

    // Only generate mips if the image is a power of 2 in both dimensions
    if (isPowerOf2(image.width) && isPowerOf2(image.height)) {
      gl.generateMipmap(gl.TEXTURE_2D);
    } else {
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    }
  };
  image.src = url;
  return texture;
}

function isPowerOf2(value) {
  return (value & (value - 1)) === 0;
}

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

  // Improved handling of diamond collisions
  // Enhanced particle ejection from diamond obstacle
  handleParticleCollisions(obstacleX, obstacleY, obstacleRadius) {
    const h = 1.0 / this.fInvSpacing;
    const r = this.particleRadius;

    // Moving circular obstacle properties
    const minDistCircle = obstacleRadius + r;
    const minDistCircle2 = minDistCircle * minDistCircle;

    // Fixed diamond obstacle properties
    const diamond = scene.fixedObstacle;

    // Domain boundaries
    const minX = h + r;
    const maxX = (this.fNumX - 1) * h - r;
    const minY = h + r;
    const maxY = (this.fNumY - 1) * h - r;

    for (let i = 0; i < this.numParticles; i++) {
      let x = this.particlePos[2 * i];
      let y = this.particlePos[2 * i + 1];

      // Check collision with moving circular obstacle (keeping this the same)
      const dxCircle = x - obstacleX;
      const dyCircle = y - obstacleY;
      const d2Circle = dxCircle * dxCircle + dyCircle * dyCircle;

      if (d2Circle < minDistCircle2 && d2Circle > 0.0) {
        const d = Math.sqrt(d2Circle);

        // Move particle to surface
        x = obstacleX + dxCircle * (minDistCircle / d);
        y = obstacleY + dyCircle * (minDistCircle / d);

        // Reflect velocity with damping
        const dot =
          dxCircle * this.particleVel[2 * i] +
          dyCircle * this.particleVel[2 * i + 1];
        const damping = 0.7;
        this.particleVel[2 * i] =
          (this.particleVel[2 * i] -
            ((1.0 + damping) * dot * dxCircle) / d2Circle) *
          damping;
        this.particleVel[2 * i + 1] =
          (this.particleVel[2 * i + 1] -
            ((1.0 + damping) * dot * dyCircle) / d2Circle) *
          damping;
      }

      // Enhanced diamond collision and trapped particle handling
      if (diamond && diamond.show) {
        // First, check if the particle is inside the diamond
        const insideDiamond = isPointInDiamond(x, y, diamond);

        // Inside handleParticleCollisions where we define the vertices for the diamond
        const vertices = [
          {
            x: diamond.x - 0.28 * diamond.size,
            y: diamond.y + 0.25 + 0.3 * diamond.size,
          },
          { x: diamond.x, y: diamond.y + 0.25 + 0.3 * diamond.size },
          {
            x: diamond.x + 0.28 * diamond.size,
            y: diamond.y + 0.25 + 0.3 * diamond.size,
          },
          { x: diamond.x + 0.4 * diamond.size, y: diamond.y + 0.25 },
          {
            x: diamond.x + 0.2 * diamond.size,
            y: diamond.y + 0.25 - 0.35 * diamond.size,
          },
          { x: diamond.x, y: diamond.y + 0.25 - 0.7 * diamond.size },
          {
            x: diamond.x - 0.2 * diamond.size,
            y: diamond.y + 0.25 - 0.35 * diamond.size,
          },
          { x: diamond.x - 0.4 * diamond.size, y: diamond.y + 0.25 },
          {
            x: diamond.x - 0.28 * diamond.size,
            y: diamond.y + 0.25 + 0.3 * diamond.size,
          },
        ];

        // const vertices = [
        //   {
        //     x: diamond.x - 0.58 * diamond.size,
        //     y: diamond.y + 0.3 + 0.3 * diamond.size,
        //   }, // Top-left
        //   { x: diamond.x, y: diamond.y + 0.3 + 0.3 * diamond.size }, // Top-middle
        //   {
        //     x: diamond.x + 0.58 * diamond.size,
        //     y: diamond.y + 0.3 + 0.3 * diamond.size,
        //   }, // Top-right
        //   { x: diamond.x + 0.77 * diamond.size, y: diamond.y + 0.32 }, // Right point
        //   {
        //     x: diamond.x + 0.33 * diamond.size,
        //     y: diamond.y + 0.32 - 0.35 * diamond.size,
        //   }, // Bottom-right
        //   { x: diamond.x, y: diamond.y + 0.35 - 0.7 * diamond.size }, // Bottom point
        //   {
        //     x: diamond.x - 0.33 * diamond.size,
        //     y: diamond.y + 0.32 - 0.35 * diamond.size,
        //   }, // Bottom-left
        //   { x: diamond.x - 0.77 * diamond.size, y: diamond.y + 0.32 }, // Left point
        //   {
        //     x: diamond.x - 0.58 * diamond.size,
        //     y: diamond.y + 0.3 + 0.3 * diamond.size,
        //   }, // Back to first
        // ];

        if (insideDiamond) {
          // Particle is trapped inside - we need to eject it

          // Find the closest edge to eject the particle
          const diamondCenter = { x: diamond.x, y: diamond.y + 0.25 }; // Adjusted for yTransform

          // Find the closest edge
          let minDist = Infinity;
          let closestPoint = { x: 0, y: 0 };

          for (let j = 0; j < vertices.length - 1; j++) {
            const v1 = vertices[j];
            const v2 = vertices[j + 1];

            // Calculate closest point on line segment
            const edgeX = v2.x - v1.x;
            const edgeY = v2.y - v1.y;
            const edgeLengthSquared = edgeX * edgeX + edgeY * edgeY;

            let t =
              ((x - v1.x) * edgeX + (y - v1.y) * edgeY) / edgeLengthSquared;
            t = Math.max(0, Math.min(1, t)); // Clamp to segment

            const pointOnEdge = {
              x: v1.x + t * edgeX,
              y: v1.y + t * edgeY,
            };

            // Distance to this point
            const dx = x - pointOnEdge.x;
            const dy = y - pointOnEdge.y;
            const distSquared = dx * dx + dy * dy;

            if (distSquared < minDist) {
              minDist = distSquared;
              closestPoint = pointOnEdge;
            }
          }

          // Safety buffer to ensure the particle is pushed outside
          const safety = this.particleRadius * 2.5;

          // Calculate direction from closest point to particle
          const dx = x - closestPoint.x;
          const dy = y - closestPoint.y;
          const dist = Math.sqrt(dx * dx + dy * dy);

          if (dist < 0.001) {
            // If particle is exactly on edge, push outward from center
            const dxFromCenter = x - diamond.x;
            const dyFromCenter = y - (diamond.y + 0.25);
            const distFromCenter = Math.sqrt(
              dxFromCenter * dxFromCenter + dyFromCenter * dyFromCenter,
            );

            // Place on edge and push out slightly
            if (distFromCenter > 0.001) {
              x = closestPoint.x + (dxFromCenter / distFromCenter) * safety;
              y = closestPoint.y + (dyFromCenter / distFromCenter) * safety;
            } else {
              // Random ejection if at center
              const angle = Math.random() * 2 * Math.PI;
              x = diamond.x + Math.cos(angle) * diamond.size * 0.6;
              y = diamond.y + 0.25 + Math.sin(angle) * diamond.size * 0.6;
            }
          } else {
            // Normal case - push along the vector from closest point to particle
            x = closestPoint.x + (dx / dist) * safety;
            y = closestPoint.y + (dy / dist) * safety;
          }

          // Reflect velocity outward
          const normalX = (x - closestPoint.x) / safety;
          const normalY = (y - closestPoint.y) / safety;
          const dotProd =
            this.particleVel[2 * i] * normalX +
            this.particleVel[2 * i + 1] * normalY;

          // Only reflect if moving into the surface
          if (dotProd < 0) {
            const damping = 0.5; // Reduced damping for trapped particles
            this.particleVel[2 * i] =
              (this.particleVel[2 * i] - 1.5 * dotProd * normalX) * damping;
            this.particleVel[2 * i + 1] =
              (this.particleVel[2 * i + 1] - 1.5 * dotProd * normalY) * damping;
          } else {
            // Add a small outward impulse if not already moving outward
            this.particleVel[2 * i] += normalX * 0.5;
            this.particleVel[2 * i + 1] += normalY * 0.5;
          }
        } else {
          // Instead of a binary in/out test, compute distance to closest edge
          const diamondCenter = { x: diamond.x, y: diamond.y + 0.25 }; // Adjusted for yTransform

          // Calculate distance to each edge and find the closest one
          let minDist = Infinity;
          let closestEdgeNormal = { x: 0, y: 0 };

          for (let j = 0; j < vertices.length - 1; j++) {
            const v1 = vertices[j];
            const v2 = vertices[j + 1];

            // Calculate distance to line segment
            const edgeX = v2.x - v1.x;
            const edgeY = v2.y - v1.y;
            const edgeLengthSquared = edgeX * edgeX + edgeY * edgeY;

            // Find closest point on line segment
            let t =
              ((x - v1.x) * edgeX + (y - v1.y) * edgeY) / edgeLengthSquared;
            t = Math.max(0, Math.min(1, t)); // Clamp to segment

            const closestX = v1.x + t * edgeX;
            const closestY = v1.y + t * edgeY;

            // Distance to closest point
            const dx = x - closestX;
            const dy = y - closestY;
            const distSquared = dx * dx + dy * dy;

            if (distSquared < minDist) {
              minDist = distSquared;

              // Calculate edge normal (perpendicular to edge)
              const edgeLength = Math.sqrt(edgeLengthSquared);
              closestEdgeNormal.x = -edgeY / edgeLength; // Rotate 90 degrees
              closestEdgeNormal.y = edgeX / edgeLength;

              // Make sure normal points outward
              const dotProduct =
                dx * closestEdgeNormal.x + dy * closestEdgeNormal.y;
              if (dotProduct < 0) {
                closestEdgeNormal.x = -closestEdgeNormal.x;
                closestEdgeNormal.y = -closestEdgeNormal.y;
              }
            }
          }

          // Check if we need to handle a collision
          const minDistance = Math.sqrt(minDist);
          if (minDistance < this.particleRadius * 1.5) {
            // Calculate penetration depth
            const penetrationDepth = this.particleRadius * 1.5 - minDistance;

            // Move particle along the normal by penetration depth
            x += closestEdgeNormal.x * penetrationDepth;
            y += closestEdgeNormal.y * penetrationDepth;

            // Reflect velocity with reduced damping
            const dot =
              this.particleVel[2 * i] * closestEdgeNormal.x +
              this.particleVel[2 * i + 1] * closestEdgeNormal.y;
            const damping = 0.6; // Even gentler damping for diamond

            this.particleVel[2 * i] =
              (this.particleVel[2 * i] -
                (1.0 + damping) * dot * closestEdgeNormal.x) *
              damping;
            this.particleVel[2 * i + 1] =
              (this.particleVel[2 * i + 1] -
                (1.0 + damping) * dot * closestEdgeNormal.y) *
              damping;
          }
        }
      }

      // Wall collisions (same as before)
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

      // for (let i = 0; i < this.fNumX; i++) {
      //   for (let j = 0; j < this.fNumY; j++) {
      //     const cellX = (i + 0.5) * this.h;
      //     const cellY = (j + 0.5) * this.h;

      //     // Check if cell is inside diamond
      //     if (
      //       scene.fixedObstacle &&
      //       scene.fixedObstacle.show &&
      //       isPointInDiamond(cellX, cellY, scene.fixedObstacle)
      //     ) {
      //       this.u[i * n + j] = 0.0;
      //       this.v[i * n + j] = 0.0;
      //       if (i > 0) this.u[(i - 1) * n + j] = 0.0;
      //       if (j > 0) this.v[i * n + (j - 1)] = 0.0;
      //     }
      //   }
      // }
    }

    // if (!toGrid) {
    //   // ... existing code ...

    //   if (denom > 0.0) {
    //     const picV = /* ... */;
    //     const corr = /* ... */;
    //     const flipV = v + corr;

    //     // Add small damping factor (0.99 means 1% energy loss per frame)
    //     const damping = 0.99;
    //     this.particleVel[2 * i + component] =
    //       ((1.0 - flipRatio) * picV + flipRatio * flipV) * damping;
    //   }
    // }

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

    // // Mark cells inside the diamond as solid
    // const f = scene.fluid;
    // // const n = f.fNumY;

    // // Before simulation, mark diamond cells as solid
    // for (let i = 1; i < f.fNumX - 2; i++) {
    //   for (let j = 1; j < f.fNumY - 2; j++) {
    //     const cellX = (i + 0.5) * f.h;
    //     const cellY = (j + 0.5) * f.h;

    //     if (isPointInDiamond(cellX, cellY, scene.fixedObstacle)) {
    //       f.cellType[i * n + j] = SOLID_CELL; // Mark as solid cell
    //     }
    //   }
    // }
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
  flipRatio: 0.45,
  numPressureIters: 100,
  numParticleIters: 2,
  frameNr: 0,
  overRelaxation: 1.9,
  compensateDrift: true,
  separateParticles: true,
  obstacleX: 0.0,
  obstacleY: 0.0,
  // obstacleRadius: 0.15,
  obstacleRadius: 0.35,
  paused: false,
  showObstacle: true,
  obstacleVelX: 0.0,
  obstacleVelY: 0.0,
  showParticles: true,
  showGrid: false,
  fluid: null,
  fixedObstacle: {
    // x: simWidth / 2,
    // y: simHeight / 2,
    x: 0.0,
    y: 0.0,
    size: 0.65, // Half-width/height of the diamond
    show: false,
    color: [1.0, 1.0, 1.0], // White color
    textureIntensity: 1.0,
  },
};

// Add to setupScene after creating the fluid
function cleanupTrappedParticles() {
  const f = scene.fluid;
  const diamond = scene.fixedObstacle;

  if (!diamond || !diamond.show) return;

  // Check all particles for being inside the diamond
  for (let i = 0; i < f.numParticles; i++) {
    const x = f.particlePos[2 * i];
    const y = f.particlePos[2 * i + 1];

    if (isPointInDiamond(x, y, diamond)) {
      // Move the particle to a safe position outside
      const dx = x - diamond.x;
      const dy = y - (diamond.y + 0.25);
      const dist = Math.sqrt(dx * dx + dy * dy);

      // Push out along the ray from center
      if (dist > 0.001) {
        const safeRadius = diamond.size * 1.1; // Safe distance outside
        f.particlePos[2 * i] = diamond.x + (dx / dist) * safeRadius;
        f.particlePos[2 * i + 1] = diamond.y + 0.25 + (dy / dist) * safeRadius;
      } else {
        // Random position if at center
        const angle = Math.random() * 2 * Math.PI;
        f.particlePos[2 * i] = diamond.x + Math.cos(angle) * diamond.size * 1.1;
        f.particlePos[2 * i + 1] =
          diamond.y + 0.25 + Math.sin(angle) * diamond.size * 1.1;
      }

      // Zero out velocity to prevent immediate re-entry
      f.particleVel[2 * i] = 0;
      f.particleVel[2 * i + 1] = 0;
    }
  }
}

function setupScene() {
  scene.obstacleRadius = 0.35;
  scene.overRelaxation = 1.9;
  scene.dt = 1.0 / 60.0;
  scene.numPressureIters = 50;
  scene.numParticleIters = 2;

  const res = 100;
  const tankHeight = 1.0 * simHeight;
  const tankWidth = 1.0 * simWidth;
  const h = (tankHeight / res) * 5.0; // Changed by me
  const density = 1000.0;
  const relWaterHeight = 0.75; // Change to change overall water height
  const relWaterWidth = 0.8;

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

  // function initTextures() {
  diamondTexture = loadTexture(gl, 'diamond_vector.svg'); // Path to your SVG
  // }

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

  // For fixed obstacle
  scene.fixedObstacle.x = simWidth / 2;
  scene.fixedObstacle.y = simHeight / 2;
  scene.fixedObstacle.size = 1;
  scene.fixedObstacle.show = true;

  setObstacle(3.0, 2.0, true);
  setupFixedObstacle();
  // Call this at the end of your setupScene function:
  cleanupTrappedParticles();
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

const texturedParticleVertexShader = `
  attribute vec2 attrPosition;
  attribute vec3 attrColor;
  uniform vec2 domainSize;
  uniform float pointSize;

  varying vec3 fragColor;

  void main() {
    vec4 screenTransform = 
      vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);
    gl_Position =
      vec4(attrPosition * screenTransform.xy + screenTransform.zw, 0.0, 1.0);

    gl_PointSize = pointSize;
    fragColor = attrColor;
  }
`;

const texturedParticleFragmentShader = `
  precision mediump float;
  varying vec3 fragColor;
  uniform sampler2D texture;
  uniform float logoIntensity;

  void main() {
    // Get texture color using gl_PointCoord (automatically provided for point sprites)
    vec4 texColor = texture2D(texture, gl_PointCoord);
    
    // Only show opaque parts of the logo
    if (texColor.a < 0.1) discard;
    
    // Mix particle color with logo texture
    vec3 finalColor = mix(fragColor, texColor.rgb, logoIntensity);
    
    gl_FragColor = vec4(finalColor, 1.0);
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

// Vertex shader program
const diamondVertexShader = `
attribute vec4 aVertexPosition;
attribute vec4 aVertexColor;

varying lowp vec4 vColor;

void main() {
    gl_Position = aVertexPosition;
    vColor = aVertexColor;
}
`;

// Fragment shader program
const diamondFratmentShader = `
varying lowp vec4 vColor;

void main() {
    gl_FragColor = vColor;
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
// Add these with your other WebGL buffer declarations
let diamondVertBuffer = null;
let diamondIdBuffer = null;
let diamondTexture = null;
let diamondTexCoordsBuffer = null;
let texturedDiamondShader = null;

function createDiamondShape() {
  // Create vertices for the diamond shape as shown in the image
  // Even narrower than before
  const vertices = [];
  const colors = [];

  const yTransform = 0.25;

  // Center vertex (for a fan-like shape)
  vertices.push(0.0, 0.0 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0); // White center

  // Define the diamond shape:
  // Top-left vertex - even narrower
  vertices.push(-0.28, 0.3 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  // Top-middle (flat top)
  vertices.push(0.0, 0.3 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  // Top-right vertex - even narrower
  vertices.push(0.28, 0.3 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  // Right point - even narrower
  vertices.push(0.4, 0.0 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  // Bottom-right diagonal - even narrower
  vertices.push(0.2, -0.35 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  // Bottom point - same height
  vertices.push(0.0, -0.7 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  // Bottom-left diagonal - even narrower
  vertices.push(-0.2, -0.35 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  // Left point - even narrower
  vertices.push(-0.4, 0.0 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  // Back to the first vertex to close the shape
  vertices.push(-0.28, 0.3 + yTransform, 0.0);
  colors.push(1.0, 1.0, 1.0, 1.0);

  return {
    vertices: vertices,
    colors: colors,
    vertexCount: 10, // Center + 8 points + closing point
  };
}

function isPointInDiamond(x, y, diamond) {
  // Transform coordinates relative to diamond center
  const dx = x - diamond.x;
  const dy = y - diamond.y - 0.25; // Account for yTransform in createDiamondShape

  // Use the EXACT same vertices as in createDiamondShape
  const vertices = [
    [-0.28, 0.3], // Top-left vertex
    [0.0, 0.3], // Top-middle
    [0.28, 0.3], // Top-right vertex
    [0.4, 0.0], // Right point
    [0.2, -0.35], // Bottom-right diagonal
    [0.0, -0.7], // Bottom point
    [-0.2, -0.35], // Bottom-left diagonal
    [-0.4, 0.0], // Left point
    [-0.28, 0.3], // Back to first
  ];

  // Scale vertices by the diamond size
  const scaledVertices = vertices.map((v) => [
    v[0] * diamond.size,
    v[1] * diamond.size,
  ]);

  // Ray casting algorithm for point-in-polygon test
  let inside = false;
  for (
    let i = 0, j = scaledVertices.length - 1;
    i < scaledVertices.length;
    j = i++
  ) {
    const xi = scaledVertices[i][0],
      yi = scaledVertices[i][1];
    const xj = scaledVertices[j][0],
      yj = scaledVertices[j][1];

    const intersect =
      yi > dy != yj > dy && dx < ((xj - xi) * (dy - yi)) / (yj - yi) + xi;
    if (intersect) inside = !inside;
  }

  return inside;
}

// function setupFixedObstacle() {
//   // The diamond properties are already set in scene.fixedObstacle
//   // Mark cells inside the diamond as solid in the fluid grid
//   const f = scene.fluid;
//   const n = f.fNumY;

//   for (let i = 1; i < f.fNumX - 2; i++) {
//     for (let j = 1; j < f.fNumY - 2; j++) {
//       const cellX = (i + 0.5) * f.h;
//       const cellY = (j + 0.5) * f.h;

//       if (isPointInDiamond(cellX, cellY, scene.fixedObstacle)) {
//         f.s[i * n + j] = 0.0; // Mark as solid cell
//       }
//     }
//   }
// }

function setupFixedObstacle() {
  const f = scene.fluid;
  const n = f.fNumY;
  const diamond = scene.fixedObstacle;

  // Define diamond vertices (same as in createDiamondShape)
  const vertices = [
    [-0.28, 0.3], // Top-left vertex
    [0.0, 0.3], // Top-middle
    [0.28, 0.3], // Top-right vertex
    [0.4, 0.0], // Right point
    [0.2, -0.35], // Bottom-right diagonal
    [0.0, -0.7], // Bottom point
    [-0.2, -0.35], // Bottom-left diagonal
    [-0.4, 0.0], // Left point
  ].map((v) => ({
    x: diamond.x + v[0] * diamond.size,
    y: diamond.y + 0.25 + v[1] * diamond.size,
  }));

  // Mark cells inside the diamond as solid
  for (let i = 1; i < f.fNumX - 2; i++) {
    for (let j = 1; j < f.fNumY - 2; j++) {
      const cellX = (i + 0.5) * f.h;
      const cellY = (j + 0.5) * f.h;

      if (isPointInDiamond(cellX, cellY, diamond)) {
        f.s[i * n + j] = 0.0; // Mark as solid cell
      }
    }
  }
}

function drawDiamondObstacle() {
  if (!scene.fixedObstacle.show || !diamondTexture) return;

  // Creates a shader of the given type, uploads the source and compiles it
  function loadShader(gl, type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);

    // Check if the compilation failed
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      alert(
        'An error occurred compiling the shaders: ' +
          gl.getShaderInfoLog(shader),
      );
      gl.deleteShader(shader);
      return null;
    }

    return shader;
  }

  // Initialize shader program
  function initShaderProgram(gl, vsSource, fsSource) {
    const vertexShader = loadShader(gl, gl.VERTEX_SHADER, vsSource);
    const fragmentShader = loadShader(gl, gl.FRAGMENT_SHADER, fsSource);

    // Create the shader program
    const shaderProgram = gl.createProgram();
    gl.attachShader(shaderProgram, vertexShader);
    gl.attachShader(shaderProgram, fragmentShader);
    gl.linkProgram(shaderProgram);

    // Check if creating the shader program failed
    if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {
      alert(
        'Unable to initialize the shader program: ' +
          gl.getProgramInfoLog(shaderProgram),
      );
      return null;
    }

    return shaderProgram;
  }

  // Initialize the shader program
  const shaderProgram = createShader(
    gl,
    diamondVertexShader,
    diamondFratmentShader,
  );

  // Create the shader data
  const programInfo = {
    program: shaderProgram,
    attribLocations: {
      vertexPosition: gl.getAttribLocation(shaderProgram, 'aVertexPosition'),
      vertexColor: gl.getAttribLocation(shaderProgram, 'aVertexColor'),
    },
  };

  // console.log(programInfo);

  // Create the diamond shape data
  const diamondData = createDiamondShape();

  // Clear the canvas
  // gl.clearColor(0.2, 0.2, 0.2, 1.0); // Dark gray background to match the image
  gl.clear(gl.COLOR_BUFFER_BIT);

  // Use the shader program
  gl.useProgram(programInfo.program);

  // Create buffers for the diamond vertices
  const positionBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
  gl.bufferData(
    gl.ARRAY_BUFFER,
    new Float32Array(diamondData.vertices),
    gl.STATIC_DRAW,
  );

  // Create buffers for the diamond colors
  const colorBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, colorBuffer);
  gl.bufferData(
    gl.ARRAY_BUFFER,
    new Float32Array(diamondData.colors),
    gl.STATIC_DRAW,
  );

  // Set up the position attribute
  gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
  gl.vertexAttribPointer(
    programInfo.attribLocations.vertexPosition,
    3, // 3 components per vertex (x,y,z)
    gl.FLOAT, // Type of the data
    false, // Don't normalize
    0, // Stride (0 = use type and numComponents)
    0, // Offset
  );
  gl.enableVertexAttribArray(programInfo.attribLocations.vertexPosition);

  // Set up the color attribute
  gl.bindBuffer(gl.ARRAY_BUFFER, colorBuffer);
  gl.vertexAttribPointer(
    programInfo.attribLocations.vertexColor,
    4, // 4 components per color (r,g,b,a)
    gl.FLOAT, // Type of the data
    false, // Don't normalize
    0, // Stride (0 = use type and numComponents)
    0, // Offset
  );
  gl.enableVertexAttribArray(programInfo.attribLocations.vertexColor);

  // Draw the diamond as a triangle fan
  gl.drawArrays(gl.TRIANGLE_FAN, 0, diamondData.vertexCount);
}

function drawDiamondBoundary() {
  if (!scene.fixedObstacle.show) return;

  const diamond = scene.fixedObstacle;
  const vertices = [
    [-0.28, 0.3], // Top-left vertex
    [0.0, 0.3], // Top-middle
    [0.28, 0.3], // Top-right vertex
    [0.4, 0.0], // Right point
    [0.2, -0.35], // Bottom-right diagonal
    [0.0, -0.7], // Bottom point
    [-0.2, -0.35], // Bottom-left diagonal
    [-0.4, 0.0], // Left point
  ].map((v) => ({
    x: diamond.x + v[0] * diamond.size,
    y: diamond.y + 0.25 + v[1] * diamond.size,
  }));

  // Draw lines connecting vertices (using WebGL)
  gl.useProgram(meshShader);
  gl.uniform2f(
    gl.getUniformLocation(meshShader, 'domainSize'),
    simWidth,
    simHeight,
  );
  gl.uniform3f(
    gl.getUniformLocation(meshShader, 'color'),
    1.0,
    0.0,
    0.0, // Red outline
  );

  // Draw each edge of the diamond
  for (let i = 0; i < vertices.length; i++) {
    const v1 = vertices[i];
    const v2 = vertices[(i + 1) % vertices.length];

    // Draw line from v1 to v2
    gl.lineWidth(3.0);
    gl.beginPath();
    gl.moveTo(v1.x, v1.y);
    gl.lineTo(v2.x, v2.y);
    gl.stroke();
  }
}

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

  // Load texture if not already loaded
  if (!logoTexture) {
    logoTexture = loadTexture(gl, '/rare_logo.png'); // Update this path
  }

  // // Draw fixed diamond obstacle

  drawDiamondObstacle();

  // drawDiamondBoundary();

  // water
  if (scene.showParticles) {
    gl.clear(gl.DEPTH_BUFFER_BIT);
    const pointSize =
      ((2.0 * scene.fluid.particleRadius) / simWidth) * canvas.width * 0.75;

    // Use textured shader
    const texturedShader = createShader(
      gl,
      texturedParticleVertexShader,
      texturedParticleFragmentShader,
    );
    gl.useProgram(texturedShader);

    // Set uniforms
    gl.uniform2f(
      gl.getUniformLocation(texturedShader, 'domainSize'),
      simWidth,
      simHeight,
    );
    gl.uniform1f(gl.getUniformLocation(texturedShader, 'pointSize'), pointSize);
    gl.uniform1f(gl.getUniformLocation(texturedShader, 'logoIntensity'), 0.8);

    // Bind texture
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, logoTexture);
    gl.uniform1i(gl.getUniformLocation(texturedShader, 'texture'), 0);

    // Set up buffers (same as before)
    if (pointVertexBuffer === null) pointVertexBuffer = gl.createBuffer();
    if (pointColorBuffer === null) pointColorBuffer = gl.createBuffer();

    gl.bindBuffer(gl.ARRAY_BUFFER, pointVertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.particlePos, gl.DYNAMIC_DRAW);

    const posLoc = gl.getAttribLocation(texturedShader, 'attrPosition');
    gl.enableVertexAttribArray(posLoc);
    gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, pointColorBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.particleColor, gl.DYNAMIC_DRAW);

    const colorLoc = gl.getAttribLocation(texturedShader, 'attrColor');
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
  // // Keep obstacle within bounds
  // const r = scene.obstacleRadius;
  // x = Math.max(r, Math.min(simWidth - r, x));
  // y = Math.max(r, Math.min(simHeight - r, y));

  let vx = 0.0;
  let vy = 0.0;

  if (!reset) {
    vx = (x - scene.obstacleX) / scene.dt;
    vy = (y - scene.obstacleY) / scene.dt;
  }

  // let vx = 0.0;
  // let vy = 0.0;

  // if (!reset) {
  //   vx = (x - scene.obstacleX) / scene.dt;
  //   vy = (y - scene.obstacleY) / scene.dt;
  // }

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

let prevMouseX = null;
let prevMouseY = null;
let prevTime = null;
let mouseVelocity = 0;
const maxRadius = 0.4; // Maximum radius when moving fast
const minRadius = 0.1; // Minimum radius when barely moving
const velocitySensitivity = 0.03; // How quickly radius responds to velocity

let mouseDown;

// const startDrag = (x, y) => {
//   const bounds = canvas.getBoundingClientRect();
//   const mx = x - bounds.left - canvas.clientLeft;
//   const my = y - bounds.top - canvas.clientTop;
//   mouseDown = true;

//   const scaledX = mx / cScale;
//   const scaledY = (canvas.height - my) / cScale;

//   setObstacle(scaledX, scaledY, true);
//   scene.paused = false;
// };

const startDrag = (x, y) => {
  const bounds = canvas.getBoundingClientRect();
  const mx = x - bounds.left - canvas.clientLeft;
  const my = y - bounds.top - canvas.clientTop;
  mouseDown = true;

  prevMouseX = x;
  prevMouseY = y;
  prevTime = performance.now();

  const scaledX = mx / cScale;
  const scaledY = (canvas.height - my) / cScale;
  setObstacle(scaledX, scaledY, true);
  scene.paused = false;
};

// const endDrag = () => {
//   mouseDown = false;

//   scene.obstacleVelX = 0.0;
//   scene.obstacleVelY = 0.0;
//   0;
// };

const endDrag = () => {
  mouseDown = false;
  prevMouseX = null;
  prevMouseY = null;
  prevTime = null;

  // Smoothly shrink radius when mouse stops
  const shrinkRadius = () => {
    scene.obstacleRadius *= 0.9;
    if (scene.obstacleRadius > 0.01) {
      requestAnimationFrame(shrinkRadius);
    } else {
      scene.obstacleRadius = 0;
    }
    setObstacle(scene.obstacleX, scene.obstacleY, false);
  };
  shrinkRadius();

  scene.obstacleVelX = 0.0;
  scene.obstacleVelY = 0.0;
};

// const drag = (x, y) => {
//   console.log(mouseDown);
//   if (mouseDown) {
//     const bounds = canvas.getBoundingClientRect();
//     const mx = x - bounds.left - canvas.clientLeft;
//     const my = y - bounds.top - canvas.clientTop;
//     const scaledX = mx / cScale;
//     const scaledY = (canvas.height - my) / cScale;
//     setObstacle(scaledX, scaledY, false);
//   } else startDrag(x, y);

//   // startDrag(x, y);
// };

const drag = (x, y) => {
  const now = performance.now();

  if (prevMouseX !== null && prevMouseY !== null && prevTime !== null) {
    const dt = (now - prevTime) / 1000; // Convert to seconds
    if (dt > 0) {
      const dx = (x - prevMouseX) / cScale;
      const dy = (y - prevMouseY) / cScale;
      const distance = Math.sqrt(dx * dx + dy * dy);
      mouseVelocity = distance / dt;

      // Calculate dynamic radius - scales with velocity but has smoothing
      const targetRadius = Math.min(
        maxRadius,
        minRadius + mouseVelocity * velocitySensitivity,
      );
      scene.obstacleRadius = scene.obstacleRadius * 0.7 + targetRadius * 0.3; // Smoothing
    }
  }

  prevMouseX = x;
  prevMouseY = y;
  prevTime = now;

  if (mouseDown) {
    const bounds = canvas.getBoundingClientRect();
    const mx = x - bounds.left - canvas.clientLeft;
    const my = y - bounds.top - canvas.clientTop;
    const scaledX = mx / cScale;
    const scaledY = (canvas.height - my) / cScale;
    setObstacle(scaledX, scaledY, false);
  } else startDrag(x, y);
};

canvas.addEventListener('mousedown', (event) => {
  startDrag(event.x, event.y);
});

canvas.addEventListener('mouseup', () => {
  endDrag();
});

// let prevMouse = { x: 0, y: 0 };

let mouseMoveTimer;

canvas.addEventListener('mousemove', (event) => {
  clearTimeout(mouseMoveTimer);
  drag(event.x, event.y);

  mouseMoveTimer = setTimeout(() => {
    // console.log('hello');
    // startDrag(event.x, event.y);
    endDrag();
  }, 100);
});

// Event handlers for scroll and mouse enter

canvas.addEventListener('mouseenter', (event) => {
  startDrag(event.x, event.y);
});

canvas.addEventListener('mouseleave', (event) => {
  // console.log(event.x);
  // endDrag();
  startDrag(event.x, event.y);
});

// Handle chaos on focus and blur
window.addEventListener('blur', () => {
  // console.log('Blur');
  // canvas.style.pointerEvents = 'none';
  endDrag();
  // console.log(canvas.style.pointerEvents);
});

window.addEventListener('focus', () => {
  // console.log('Focus');
  // canvas.style.pointerEvents = 'auto';
  endDrag();
  // console.log(canvas.style.pointerEvents);
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
    // scene.fluid.simulate(
    //   scene.dt,
    //   scene.gravity,
    //   scene.flipRatio,
    //   scene.numPressureIters,
    //   scene.numParticleIters,
    //   scene.overRelaxation,
    //   scene.compensateDrift,
    //   scene.separateParticles,
    //   scene.obstacleX,
    //   scene.obstacleY,
    //   scene.obstacleRadius,
    // );
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
  // Call this at the end of your setupScene function:
  cleanupTrappedParticles();
  requestAnimationFrame(update);
};

// Responsive
window.addEventListener('resize', () => {
  // console.log('resize');
  canvas.width = window.innerWidth - 20;
  canvas.height = window.innerHeight - 20;

  cScale = canvas.height / simHeight;
  simWidth = canvas.width / cScale;

  setupScene();
});

setupScene();
update();
