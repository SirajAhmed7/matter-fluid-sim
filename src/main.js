import './style.css';

import Matter from 'matter-js';

const matterContainer = document.querySelector('.matter');

// module aliases
const Engine = Matter.Engine,
  Render = Matter.Render,
  Runner = Matter.Runner,
  Bodies = Matter.Bodies,
  Body = Matter.Body,
  Events = Matter.Events,
  Constraint = Matter.Constraint,
  Composite = Matter.Composite;

// Matter.use('matter-attractors');

// create an engine
const engine = Engine.create();
engine.world.gravity.y = 0.7;

const THICKNESS = 60;

// create a renderer
const render = Render.create({
  element: matterContainer,
  engine: engine,
  options: {
    width: matterContainer.clientWidth,
    height: matterContainer.scrollHeight,
    background: 'transparent',
    wireframes: false,
    showAngleIndicator: true,
  },
});

const bigCircle = Bodies.circle(500, 500, 64, {
  // restitution: 0.1,
  render: {
    opacity: 0,
    visible: false,
  },
  // plugin: {
  //   attractors: [
  //     function (bodyA, bodyB) {
  //       return {
  //         x: (bodyA.position.x - bodyB.position.x) * 1e-8,
  //         y: (bodyA.position.y - bodyB.position.y) * 1e-8,
  //       };
  //     },
  //   ],
  // },
});

Composite.add(engine.world, bigCircle);

let mouse = Matter.Mouse.create(render.canvas);
let mouseConstraint = Matter.MouseConstraint.create(engine, {
  mouse,
  constraint: {
    stiffness: 0.2,
    render: {
      visible: false,
    },
  },
});
Composite.add(engine.world, mouseConstraint);

// mouseConstraint.pointA = mouse.position;
// mouseConstraint.bodyB = mouseConstraint.body = bigCircle;
// mouseConstraint.pointB = {
//   x: mouse.position.x - bigCircle.position.x,
//   y: mouse.position.y - bigCircle.position.y,
// };
// mouseConstraint.angleB = bigCircle.angle;

Matter.MouseConstraint.update = function (mouseConstraint, bodies) {
  const mcMouse = mouseConstraint.mouse,
    constraint = mouseConstraint.constraint,
    body = Composite.allBodies(engine.world)[0]; // this is where I find my target Body

  constraint.pointA = mcMouse.position;
  constraint.bodyB = mouseConstraint.body = body;
  // removing declaration of constraint.pointB worked
  constraint.angleB = body.angle;
};

// render.mouse = mouse;

// create two boxes and a ground
// const boxA = Bodies.rectangle(400, 200, 80, 80);
// const boxB = Bodies.rectangle(450, 50, 80, 80);

let x = matterContainer.clientWidth / 2;
const y = matterContainer.clientHeight * 0.2;

// function createLiquidParticles(x, y, particleCount, particleRadius) {
//   const particles = [];
//   const constraints = [];

//   // Create individual particles
//   for (let i = 0; i < particleCount; i++) {
//     const particle = Bodies.circle(
//       x + Math.random() * 200,
//       y + Math.random() * 200,
//       particleRadius,
//       {
//         friction: 0.05,
//         restitution: 0.0005,
//       },
//     );
//     particles.push(particle);
//   }

//   // Create constraints between nearby particles
//   for (let i = 0; i < particles.length; i++) {
//     for (let j = i + 1; j < particles.length; j++) {
//       const distance = Matter.Vector.magnitude(
//         Matter.Vector.sub(particles[i].position, particles[j].position),
//       );

//       if (distance < particleRadius * 2.5) {
//         const constraint = Constraint.create({
//           bodyA: particles[i],
//           bodyB: particles[j],
//           stiffness: 0.1,
//           damping: 0.05,
//         });
//         constraints.push(constraint);
//       }
//     }
//   }

//   return { particles, constraints };
// }

// const liquidSim = createLiquidParticles(400, 200, 1500, 10);
const particleCount = 800;

for (let i = 0; i < particleCount; i++) {
  x = (matterContainer.clientWidth / particleCount) * i;

  const collisionShape = Bodies.circle(x, y, 14, {
    isSensor: false,
    render: {
      opacity: 0,
      visible: false,
    },
  });

  const collisionMaterialOptions = {
    friction: 0.1,
    frictionAir: 0.001,
    restitution: 0.001,
    density: 0.00005,
    timeScale: 1,
    slop: 0.05,
    // plugin: {
    //   attractors: [
    //     function (bodyA, bodyB) {
    //       return {
    //         x: (bodyA.position.x - bodyB.position.x) * 0.000000001,
    //         y: (bodyA.position.y - bodyB.position.y) * 0.000000001,
    //       };
    //     },
    //   ],
    // },
    collisionFilter: {
      group: 0,
      category: 0x0001,
      mask: 0xffffffff,
    },
    // sleepThreshold: -1,
  };

  let shape;
  const shapeType = Math.floor(Math.random() * 3);

  // const materialOptions = {
  //   friction: 0.05,
  //   frictionAir: 0.0001,
  //   restitution: 1,
  //   density: 0.001,
  //   sleepThreshold: -1,
  // };

  const shapesColor = '#222';

  const size = 10;

  // if (i % 3 === 0) {
  if (shapeType === 0) {
    const horizontalRect = Bodies.rectangle(x, y, size, size / 4, {
      render: { fillStyle: shapesColor },
    });

    // Create a vertical rectangle
    const verticalRect = Bodies.rectangle(x, y, size / 4, size, {
      render: { fillStyle: shapesColor },
    });

    // Combine them into a compound body
    shape = Body.create({
      parts: [collisionShape, horizontalRect, verticalRect],
      ...collisionMaterialOptions,
    });

    // shape = Bodies.fromVertices(x, y, [plusVertices], materialOptions);
    // } else if (i === 0 || i % 2 === 0) {
  } else if (shapeType === 1) {
    const circle = Bodies.circle(x, y, size / 2, {
      render: { fillStyle: shapesColor },
    });

    shape = Body.create({
      parts: [collisionShape, circle],
      ...collisionMaterialOptions,
    });
    // } else if (i % 2 === 1) {
  } else if (shapeType === 2) {
    const square = Bodies.rectangle(x, y, size, size, {
      render: { fillStyle: shapesColor },
    });

    shape = Body.create({
      parts: [collisionShape, square],
      ...collisionMaterialOptions,
    });
  }

  // setTimeout(() => {
  //   Composite.add(engine.world, shape);
  // }, i * 8);
  Composite.add(engine.world, shape);
}

const ground = Bodies.rectangle(
  matterContainer.clientWidth / 2,
  matterContainer.clientHeight + THICKNESS / 2,
  // matterContainer.clientHeight,
  20000,
  THICKNESS,
  {
    isStatic: true,
  },
);

const ceil = Bodies.rectangle(
  matterContainer.clientWidth / 2,
  0 - THICKNESS / 2,
  20000,
  THICKNESS,
  {
    isStatic: true,
  },
);

const leftWall = Bodies.rectangle(
  0 - THICKNESS / 2,
  matterContainer.clientHeight / 2,
  // matterContainer.clientHeight,
  THICKNESS,
  matterContainer.clientHeight,
  {
    isStatic: true,
  },
);
const rightWall = Bodies.rectangle(
  matterContainer.clientWidth + THICKNESS / 2,
  matterContainer.clientHeight / 2,
  // matterContainer.clientHeight,
  THICKNESS,
  matterContainer.clientHeight,
  {
    isStatic: true,
  },
);

// add all of the bodies to the world
// Composite.add(engine.world, [boxA, boxB, ground, leftWall, rightWall]);
Composite.add(engine.world, [
  ground,
  ceil,
  leftWall,
  rightWall,
  // ...liquidSim.particles,
  // liquidSim.constraints,
]);

// // Track mouse position and velocity
// let lastMousePosition = { x: 0, y: 0 };
// let mouseVelocity = { x: 0, y: 0 };
// let lastTimestamp = 0;

// // Store original scale (1 is default)
// let currentScale = 0;

// function updateShapeScale(timestamp) {
//   const difference = Math.min(
//     Math.abs(mouse.position.x - lastMousePosition.x),
//     150
//   );

//   const scaleFactor = Math.max(0.5, Math.min(5, difference * 0.1));
//   currentScale = scaleFactor;
//   // console.log(scaleFactor);

//   Matter.Body.setPosition(
//     bigCircle,
//     Matter.Vector.create(mouse.position.x, mouse.position.y)
//   );

//   lastMousePosition.x = mouse.position.x;
//   lastMousePosition.y = mouse.position.y;

//   // if (!lastTimestamp) {
//   //   lastTimestamp = timestamp;
//   //   requestAnimationFrame(updateShapeScale);
//   //   return;
//   // }
//   // const deltaTime = (timestamp - lastTimestamp) / 1000; // Convert to seconds
//   // lastTimestamp = timestamp;
//   // if (deltaTime > 0) {
//   //   // Calculate mouse velocity (pixels per second)
//   //   mouseVelocity.x =
//   //     Math.abs(mouse.position.x - lastMousePosition.x) / deltaTime;
//   //   mouseVelocity.y = (mouse.position.y - lastMousePosition.y) / deltaTime;
//   //   lastMousePosition = { ...mouse.position };
//   //   // Calculate speed magnitude
//   //   const speed = mouseVelocity.x;
//   //   // Map speed to scale (adjust these values as needed)
//   //   const maxSpeed = 2000; // pixels per second
//   //   const minScale = 0.001;
//   //   const maxScale = 2;
//   //   // Smoothly transition the scale
//   //   const targetScale = Math.max(
//   //     minScale,
//   //     Math.min(maxScale, (speed / maxSpeed) * maxScale)
//   //   );
//   //   currentScale += (targetScale - currentScale) * 0.3; // Smoothing factor
//   //   // console.log(speed);
//   //   // Apply scale (set to 0 if speed is very low)
//   //   if (speed < 50) {
//   //     currentScale *= 0.9; // Quickly shrink when nearly stopped
//   //   }
//   //   console.log(currentScale);
//   //   // Update shape vertices with the new scale
//   //   Matter.Body.scale(bigCircle, currentScale, currentScale);
//   // }
//   // // console.log("hello");
//   requestAnimationFrame(updateShapeScale);
// }

// // Start the animation loop
// // requestAnimationFrame(updateShapeScale);
// updateShapeScale();

// run the renderer

Render.run(render);

// Events.on(mouseConstraint, "mousemove", (event) => {
//   // console.log(event);

//   // Store last mouse Position
//   lastMouseX = curMouseX;
//   lastMouseY = curMouseY;

//   // Store current mouse position
//   curMouseX = event.mouse.absolute.x;
//   curMouseY = event.mouse.absolute.y;

//   // Calculate mouse movement delta
//   const deltaX = curMouseX - lastMouseX;

//   // Scale factor based on mouse movement (adjust sensitivity as needed)
//   const scaleFactor = 1 + Math.abs(deltaX * 0.01);
//   lastScaleFactor = scaleFactor;
//   console.log(scaleFactor);

//   // Apply scaling to the rectangle
//   Body.scale(bigCircle, scaleFactor, scaleFactor);
//   // Only scale if the change is significant enough
// });

// const animate = function () {
//   console.log(deltaX);
//   // requestAnimationFrame(animate);
//   // if (Math.abs(deltaX) > 1) {
//   //   // Body.scale(bigCircle, scaleFactor, scaleFactor);
//   //   Matter.Body.scale(bigCircle, scaleFactor, scaleFactor);
//   // }
// };
// animate();

// create runner
const runner = Runner.create();

window.addEventListener('resize', () => {
  // Resize canvas
  render.canvas.width = matterContainer.clientWidth;
  render.canvas.height = matterContainer.clientHeight;

  // Reposition ground
  Matter.Body.setPosition(
    ground,
    Matter.Vector.create(
      matterContainer.clientWidth / 2,
      matterContainer.clientHeight + THICKNESS / 2,
    ),
  );

  Matter.Body.setPosition(
    rightWall,
    Matter.Vector.create(
      matterContainer.clientWidth + THICKNESS / 2,
      matterContainer.clientHeight / 2,
    ),
  );
});

// Allow scroll through canvas
mouse.element.removeEventListener('wheel', mouse.mousewheel);
// mouseConstraint.mouse.element.removeEventListener(
//   "mousewheel",
//   mouseConstraint.mouse.mousewheel
// );
// mouseConstraint.mouse.element.removeEventListener(
//   "DOMMouseScroll",
//   mouseConstraint.mouse.mousewheel
// );

// run the engine
Runner.run(runner, engine);
