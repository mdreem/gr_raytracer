# Rust Improvement Guide

This document provides detailed explanations and examples for the Rust idiom improvements suggested in the code review.

---

## 1. Using `.zip()` Iterators Instead of Indexed Loops

### The Problem

Indexed loops are common in languages like C/C++, but in Rust they're less idiomatic and more error-prone:

```rust
// Current pattern (found in schwarzschild.rs:543-557)
for i in 0..trajectory_a.len() {
    let step_a = &trajectory_a[i];
    let step_b = &trajectory_b[i];

    let step_a_cartesian = Point::new_spherical(
        step_a.x[0], step_a.x[1], step_a.x[2], step_a.x[3]
    ).get_spatial_vector_cartesian();

    let step_b_cartesian = Point::new_spherical(
        step_b.x[0], step_b.x[1], step_b.x[2], step_b.x[3]
    ).get_spatial_vector_cartesian();

    assert_abs_diff_eq!(step_a_cartesian[0], step_b_cartesian[0], epsilon = 1e-5);
    assert_abs_diff_eq!(step_a_cartesian[1], step_b_cartesian[2], epsilon = 1e-5);
    assert_abs_diff_eq!(step_a_cartesian[2], -step_b_cartesian[1], epsilon = 1e-5);
}
```

**Issues with this approach:**
1. **Panic risk:** If `trajectory_b.len() != trajectory_a.len()`, we'll panic with index out of bounds
2. **Less clear intent:** The reader must infer we're iterating both collections in parallel
3. **Manual indexing:** Unnecessary cognitive overhead tracking the index variable
4. **Not zero-cost:** The compiler may not optimize away bounds checks

### The Solution: Use `.zip()`

```rust
// Improved idiomatic version
for (step_a, step_b) in trajectory_a.iter().zip(trajectory_b.iter()) {
    let step_a_cartesian = Point::new_spherical(
        step_a.x[0], step_a.x[1], step_a.x[2], step_a.x[3]
    ).get_spatial_vector_cartesian();

    let step_b_cartesian = Point::new_spherical(
        step_b.x[0], step_b.x[1], step_b.x[2], step_b.x[3]
    ).get_spatial_vector_cartesian();

    assert_abs_diff_eq!(step_a_cartesian[0], step_b_cartesian[0], epsilon = 1e-5);
    assert_abs_diff_eq!(step_a_cartesian[1], step_b_cartesian[2], epsilon = 1e-5);
    assert_abs_diff_eq!(step_a_cartesian[2], -step_b_cartesian[1], epsilon = 1e-5);
}
```

**Benefits:**
1. **No panic on length mismatch:** `.zip()` stops at the shorter iterator's length
2. **Clear intent:** Immediately obvious we're processing pairs
3. **No manual indexing:** Compiler handles everything
4. **Better optimization:** No bounds checks needed

### When You Need the Index

If you also need the index, use `.enumerate()`:

```rust
// If you need both index and elements
for (i, (step_a, step_b)) in trajectory_a.iter().zip(trajectory_b.iter()).enumerate() {
    println!("Comparing step {}", i);
    // ... rest of logic
}
```

### Other Iterator Patterns in Your Codebase

#### Example 1: Windows for Adjacent Elements

Currently (in `scene.rs:96-104`):
```rust
for step_window in steps.steps.windows(2) {
    let last_step = &step_window[0];
    let step = &step_window[1];

    if let Some(intersection_color) =
        self.objects.intersects(last_step, step, observer_energy) {
        intersections.push(intersection_color);
    }
}
```

This is already good! But could be slightly cleaner:

```rust
// Alternative: explicitly pattern match the window
for window in steps.steps.windows(2) {
    let [last_step, step] = window else { unreachable!() };

    if let Some(intersection_color) =
        self.objects.intersects(last_step, step, observer_energy) {
        intersections.push(intersection_color);
    }
}

// Or even better with iterator methods:
let intersections: Vec<_> = steps.steps
    .windows(2)
    .filter_map(|window| {
        let [last_step, step] = window else { unreachable!() };
        self.objects.intersects(last_step, step, observer_energy)
    })
    .collect();
```

#### Example 2: Collecting Results

In tests (like `schwarzschild.rs:702-709`):
```rust
fn collect_points_step(steps: &IntegratedRay) -> Vec<MatchedPoint> {
    steps
        .iter()
        .map(|step| MatchedPoint {
            r: step.x[1],
            phi: step.x[3],
        })
        .collect()
}
```

This is already idiomatic! ✅ Good example of functional style.

#### Example 3: Comparing Rotated Vectors (kerr.rs:588-595)

Current:
```rust
for i in 0..trajectory_a.len() {
    let step_a = &trajectory_a[i];
    let step_b = &trajectory_b[i];

    assert_abs_diff_eq!(step_a.x[1], step_b.x[1], epsilon = 1e-5);
    assert_abs_diff_eq!(step_a.x[2], -step_b.x[3], epsilon = 1e-5);
    assert_abs_diff_eq!(step_a.x[3], step_b.x[2], epsilon = 1e-5);
}
```

Improved:
```rust
for (step_a, step_b) in trajectory_a.iter().zip(trajectory_b.iter()) {
    assert_abs_diff_eq!(step_a.x[1], step_b.x[1], epsilon = 1e-5);
    assert_abs_diff_eq!(step_a.x[2], -step_b.x[3], epsilon = 1e-5);
    assert_abs_diff_eq!(step_a.x[3], step_b.x[2], epsilon = 1e-5);
}
```

### Advanced: Assert Equal Length First

If you want to ensure both collections have the same length (which may be your intent):

```rust
assert_eq!(trajectory_a.len(), trajectory_b.len(),
           "Trajectories should have equal length");

for (step_a, step_b) in trajectory_a.iter().zip(trajectory_b.iter()) {
    // Now guaranteed to process all elements
    assert_abs_diff_eq!(step_a.x[1], step_b.x[1], epsilon = 1e-5);
    // ...
}
```

Or use a helper:

```rust
fn zip_exact<T, U>(a: &[T], b: &[U]) -> impl Iterator<Item = (&T, &U)> {
    assert_eq!(a.len(), b.len(), "Collections must have equal length");
    a.iter().zip(b.iter())
}

// Usage:
for (step_a, step_b) in zip_exact(&trajectory_a, &trajectory_b) {
    // ...
}
```

---

## 2. Compile-Time Coordinate System Safety with Phantom Types

### The Problem

Currently, coordinate systems are checked at runtime:

```rust
// From src/geometry/point.rs
pub struct Point {
    pub coordinate_system: CoordinateSystem,  // Runtime tag
    pub vector: Vector4<f64>,
}

impl Point {
    pub fn new_cartesian(t: f64, x: f64, y: f64, z: f64) -> Self {
        Point {
            coordinate_system: CoordinateSystem::Cartesian,
            vector: Vector4::new(t, x, y, z),
        }
    }

    pub fn new_spherical(t: f64, r: f64, theta: f64, phi: f64) -> Self {
        Point {
            coordinate_system: CoordinateSystem::Spherical,
            vector: Vector4::new(t, r, theta, phi),
        }
    }
}

// Usage requires runtime checks:
impl InnerProduct for Schwarzschild {
    fn inner_product(&self, position: &Point, v: &FourVector, w: &FourVector) -> f64 {
        debug_assert_eq!(position.coordinate_system, Spherical);  // Runtime check!
        // ... rest of implementation
    }
}
```

**Issues:**
1. **Runtime overhead:** Every operation must check the coordinate system
2. **Late error detection:** Mistakes only caught when code runs (or in debug mode)
3. **Easy to forget checks:** Developers might forget to validate
4. **Not zero-cost abstraction:** The tag occupies memory and requires branches

### The Solution: Phantom Types

Use Rust's type system to enforce coordinate systems at compile time:

```rust
use std::marker::PhantomData;

// Marker types for coordinate systems (zero-sized types)
pub struct Cartesian;
pub struct Spherical;

// Sealed trait to prevent users from adding new coordinate systems
mod sealed {
    pub trait Sealed {}
    impl Sealed for super::Cartesian {}
    impl Sealed for super::Spherical {}
}

pub trait CoordinateSystem: sealed::Sealed + Copy + 'static {
    const NAME: &'static str;
}

impl CoordinateSystem for Cartesian {
    const NAME: &'static str = "Cartesian";
}

impl CoordinateSystem for Spherical {
    const NAME: &'static str = "Spherical";
}

// Point is now parameterized by coordinate system
pub struct Point<CS: CoordinateSystem> {
    pub vector: Vector4<f64>,
    _phantom: PhantomData<CS>,  // Zero-sized, no runtime cost
}

impl Point<Cartesian> {
    pub fn new_cartesian(t: f64, x: f64, y: f64, z: f64) -> Self {
        Point {
            vector: Vector4::new(t, x, y, z),
            _phantom: PhantomData,
        }
    }

    // Conversion to spherical returns a different type
    pub fn to_spherical(&self) -> Point<Spherical> {
        let t = self.vector[0];
        let x = self.vector[1];
        let y = self.vector[2];
        let z = self.vector[3];

        let r = (x * x + y * y + z * z).sqrt();
        let theta = (z / r).acos();
        let phi = y.atan2(x);

        Point {
            vector: Vector4::new(t, r, theta, phi),
            _phantom: PhantomData,
        }
    }
}

impl Point<Spherical> {
    pub fn new_spherical(t: f64, r: f64, theta: f64, phi: f64) -> Self {
        Point {
            vector: Vector4::new(t, r, theta, phi),
            _phantom: PhantomData,
        }
    }

    pub fn to_cartesian(&self) -> Point<Cartesian> {
        let t = self.vector[0];
        let r = self.vector[1];
        let theta = self.vector[2];
        let phi = self.vector[3];

        let x = r * theta.sin() * phi.cos();
        let y = r * theta.sin() * phi.sin();
        let z = r * theta.cos();

        Point {
            vector: Vector4::new(t, x, y, z),
            _phantom: PhantomData,
        }
    }
}

// Common operations available regardless of coordinate system
impl<CS: CoordinateSystem> Point<CS> {
    pub fn radial_distance_spatial_part_squared(&self) -> f64 {
        let x = self.vector[1];
        let y = self.vector[2];
        let z = self.vector[3];
        x * x + y * y + z * z
    }
}

// Index trait for accessing components
impl<CS: CoordinateSystem> std::ops::Index<usize> for Point<CS> {
    type Output = f64;
    fn index(&self, index: usize) -> &f64 {
        &self.vector[index]
    }
}
```

### How Geometries Use This

Geometries are now parameterized by their coordinate system:

```rust
pub trait Geometry: InnerProduct + Signature + Clone + Sync {
    type Coords: CoordinateSystem;

    fn get_tetrad_at(&self, position: &Point<Self::Coords>) -> Tetrad<Self::Coords>;
    fn lorentz_transformation(
        &self,
        position: &Point<Self::Coords>,
        velocity: &FourVector<Self::Coords>
    ) -> Matrix4<f64>;
    // ... other methods
}

// Schwarzschild works in spherical coordinates
impl Geometry for Schwarzschild {
    type Coords = Spherical;

    fn inner_product(
        &self,
        position: &Point<Spherical>,  // Type signature enforces spherical!
        v: &FourVector<Spherical>,
        w: &FourVector<Spherical>
    ) -> f64 {
        // No runtime check needed - compiler guarantees it's spherical
        let r = position[1];
        let theta = position[2];

        let a = 1.0 - self.radius / r;

        a * v.vector[0] * w.vector[0]
            - v.vector[1] * w.vector[1] / a
            - r * r * v.vector[2] * w.vector[2]
            - r * r * theta.sin() * theta.sin() * v.vector[3] * w.vector[3]
    }
}

// Kerr works in Cartesian coordinates
impl Geometry for Kerr {
    type Coords = Cartesian;

    fn inner_product(
        &self,
        position: &Point<Cartesian>,  // Type signature enforces Cartesian!
        v: &FourVector<Cartesian>,
        w: &FourVector<Cartesian>
    ) -> f64 {
        let metric = self.metric(position[1], position[2], position[3]);
        (v.vector.transpose() * metric * w.vector).x
    }
}
```

### Benefits

**1. Compile-Time Safety:**
```rust
let p_cart = Point::new_cartesian(0.0, 1.0, 2.0, 3.0);
let p_sph = Point::new_spherical(0.0, 1.0, 0.5, 0.0);

let schwarzschild = Schwarzschild::new(2.0, 1e-4);

// This compiles:
schwarzschild.get_tetrad_at(&p_sph);  ✅

// This won't compile - caught at compile time!
schwarzschild.get_tetrad_at(&p_cart);  ❌
// error[E0308]: mismatched types
//   expected `Point<Spherical>`
//      found `Point<Cartesian>`
```

**2. Zero Runtime Cost:**
```rust
// The PhantomData<CS> has zero size
assert_eq!(
    std::mem::size_of::<Point<Cartesian>>(),
    std::mem::size_of::<Vector4<f64>>()
);

// No branching needed - the type determines the coordinate system
// Compiler can optimize aggressively
```

**3. Explicit Conversions:**
```rust
// Conversions are explicit and type-checked
let p_cart = Point::new_cartesian(0.0, 1.0, 0.0, 0.0);
let p_sph = p_cart.to_spherical();  // Returns Point<Spherical>

// Can't accidentally use wrong coordinate system
let schwarzschild = Schwarzschild::new(2.0, 1e-4);
schwarzschild.get_tetrad_at(&p_cart.to_spherical());  ✅
```

**4. Self-Documenting Code:**
```rust
// Function signatures make coordinate system requirements crystal clear
fn compute_schwarzschild_redshift(
    observer_pos: &Point<Spherical>,
    emitter_pos: &Point<Spherical>,
    geometry: &Schwarzschild,
) -> f64 {
    // No need to check - compiler guarantees spherical coordinates
    // ...
}

// Trying to pass Cartesian coordinates won't compile
let cart_pos = Point::new_cartesian(0.0, 1.0, 2.0, 3.0);
compute_schwarzschild_redshift(&cart_pos, &cart_pos, &geometry);  ❌
```

### Migration Strategy

Since this is a breaking change, here's how to migrate gradually:

#### Step 1: Add phantom types alongside existing runtime checks

```rust
pub struct Point<CS: CoordinateSystem = DynamicCS> {
    pub coordinate_system: CoordinateSystemTag,  // Keep for now
    pub vector: Vector4<f64>,
    _phantom: PhantomData<CS>,
}

// Special marker for dynamic dispatch
pub struct DynamicCS;
impl CoordinateSystem for DynamicCS {
    const NAME: &'static str = "Dynamic";
}
```

#### Step 2: Add type aliases for compatibility

```rust
// Backwards-compatible aliases
pub type DynamicPoint = Point<DynamicCS>;
pub type CartesianPoint = Point<Cartesian>;
pub type SphericalPoint = Point<Spherical>;
```

#### Step 3: Gradually convert code

Convert one module at a time, starting with geometry implementations:
1. Convert `Schwarzschild` to use `Point<Spherical>`
2. Convert `Kerr` to use `Point<Cartesian>`
3. Update tests
4. Remove runtime checks

#### Step 4: Remove dynamic dispatch

Once all code is converted, remove `DynamicCS` and the runtime tag.

### Advanced: Generic Over Coordinate Systems

Some functions work with any coordinate system:

```rust
// Generic function works with any coordinate system
fn compute_proper_time<CS: CoordinateSystem>(
    start: &Point<CS>,
    end: &Point<CS>,
) -> f64 {
    let dt = end[0] - start[0];
    // ... computation independent of coordinate system
    dt
}

// Can be called with either:
compute_proper_time(&point_cartesian, &other_cartesian);  ✅
compute_proper_time(&point_spherical, &other_spherical);  ✅

// But not mixed:
compute_proper_time(&point_cartesian, &point_spherical);  ❌
```

### Real-World Example: Camera

The camera could enforce its coordinate requirements:

```rust
impl Camera {
    pub fn new<G: Geometry>(
        position: Point<G::Coords>,  // Must match geometry's coordinates
        velocity: FourVector<G::Coords>,
        alpha: f64,
        rows: i64,
        columns: i64,
        phi: f64,
        theta: f64,
        psi: f64,
        geometry: &G,
    ) -> Result<Camera, CameraError> {
        // Type system guarantees position and velocity have correct coordinates
        let tetrad = geometry.get_tetrad_at(&position);
        // ...
    }
}

// Usage is type-safe:
let pos = Point::new_spherical(0.0, 10.0, PI/2.0, 0.0);
let vel = FourVector::new_spherical(1.0, 0.0, 0.0, 0.0);
let schwarzschild = Schwarzschild::new(2.0, 1e-4);

Camera::new(pos, vel, PI/2.0, 100, 100, 0.0, 0.0, 0.0, &schwarzschild);  ✅

// This won't compile if coordinates don't match:
let pos_cart = Point::new_cartesian(0.0, 10.0, 0.0, 0.0);
Camera::new(pos_cart, vel, PI/2.0, 100, 100, 0.0, 0.0, 0.0, &schwarzschild);  ❌
```

---

## Summary

### Iterator Pattern (.zip())

**When to use:**
- Processing multiple collections in parallel
- Comparing elements from different sequences
- Pairing related data

**Benefits:**
- Compiler-enforced safety
- Clearer intent
- Better optimization
- Less boilerplate

**Current locations to update:**
- `src/geometry/schwarzschild.rs:543-557`
- `src/geometry/kerr.rs:588-595`
- Any other indexed parallel iteration

### Phantom Types

**When to use:**
- Type-level enforcement of invariants
- Zero-runtime-cost abstractions
- Preventing invalid state at compile time

**Benefits:**
- Compile-time coordinate system checking
- Zero runtime overhead
- Self-documenting code
- Impossible to misuse API

**Trade-offs:**
- More complex type signatures
- Steeper learning curve
- Breaking change (requires migration)
- Generic code may need more type annotations

### Recommendation

1. **Iterator patterns:** Low-hanging fruit, apply immediately
2. **Phantom types:** Consider for v2.0 or major refactor

Both patterns make your code more idiomatic Rust and leverage the type system for correctness guarantees that would be impossible in languages like C++ or Python.
