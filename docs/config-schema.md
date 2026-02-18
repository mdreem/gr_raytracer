# Render Config Schema

Generated from Rust types in `src/configuration.rs`.

## RenderConfig

Top-level scene configuration loaded from TOML.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `celestial_temperature` | `number` | yes | Background emitter temperature in Kelvin. |
| `celestial_texture` | `TextureConfig` | yes | Texture mapped onto the celestial sphere/background. |
| `color_normalization` | `CIETristimulusNormalization` | yes | Normalization applied to final rendered color values. |
| `geometry_type` | `GeometryType` | yes | Geometry model used for ray integration. |
| `objects` | `array<ObjectsConfig>` | yes | Scene objects rendered in front of the celestial background. |

## CIETristimulusNormalization

Strategy used to normalize CIE XYZ tristimulus values before image output.

Variants:

### `NoNormalization`

Keep raw XYZ values unchanged.

### `Chromaticity`

Normalize XYZ values so X + Y + Z = 1 (chromaticity coordinates).

### `EqualLuminance`

Normalize XYZ values so Y = 1 (relative luminance normalization).


## GeometryType

Available spacetime geometry models.

Variants:

### `Euclidean`

Flat Euclidean geometry in Cartesian coordinates.

### `EuclideanSpherical`

Flat Euclidean geometry in spherical coordinates.

### `Schwarzschild`

Schwarzschild black hole geometry.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `horizon_epsilon` | `number` | yes | Small offset used to avoid singular behavior at the horizon. |
| `radius` | `number` | yes | Schwarzschild radius. |

### `Kerr`

Kerr rotating black hole geometry.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `a` | `number` | yes | Kerr spin parameter. |
| `horizon_epsilon` | `number` | yes | Small offset used to avoid singular behavior at the horizon. |
| `radius` | `number` | yes | Characteristic radius term. |


## ObjectsConfig

Scene object definitions.

Variants:

### `Sphere`

Spherical emitter or textured sphere.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `position` | `array` | yes | Center position in Cartesian coordinates (x, y, z). |
| `radius` | `number` | yes | Sphere radius. |
| `temperature` | `number` | yes | Emitter temperature in Kelvin. |
| `texture` | `TextureConfig` | yes | Surface texture configuration. |

### `Disc`

Thin disc emitter with inner and outer radius.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `inner_radius` | `number` | yes | Inner radius of the disc. |
| `outer_radius` | `number` | yes | Outer radius of the disc. |
| `temperature` | `number` | yes | Emitter temperature in Kelvin. |
| `texture` | `TextureConfig` | yes | Surface texture configuration. |

### `VolumetricDisc`

Participating-media disc with absorption/scattering.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `absorption` | `number` | yes | Absorption coefficient. |
| `axis` | `array | null` | no | Optional rotation axis for the disc volume (x, y, z). |
| `brightness_reference_temperature` | `number` | yes | Reference temperature used for brightness scaling. |
| `density_multiplier` | `number` | yes | Global multiplier applied to local density values. |
| `inner_radius` | `number` | yes | Inner radius of the disc volume. |
| `max_steps` | `integer` | yes | Maximum integration steps through the participating media. |
| `noise_offset` | `number` | yes | Scalar offset added to the sampled noise domain. |
| `noise_scale` | `array` | yes | Noise scaling factors (x, y, z). |
| `num_octaves` | `integer` | yes | Number of octaves used for procedural density noise. |
| `outer_radius` | `number` | yes | Outer radius of the disc volume. |
| `perlin_seed` | `integer | null` | no | Optional explicit seed for noise generation. |
| `scattering` | `number` | yes | Scattering coefficient. |
| `step_size` | `number` | yes | Integration step size through the media volume. |
| `temperature` | `number` | yes | Base emitter temperature in Kelvin. |
| `texture` | `TextureConfig` | yes | Texture sampled within the volume. |
| `thickness` | `number` | yes | Vertical thickness of the media disc. |


## TextureConfig

Texture definitions used by objects and background.

Variants:

### `Bitmap`

Sample colors from a bitmap image.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `beaming_exponent` | `number` | yes | Exponent used for relativistic beaming intensity scaling. |
| `color_normalization` | `CIETristimulusNormalization` | yes | Normalization applied to sampled texture colors. |
| `path` | `string` | yes | Filesystem path to the image file. |

### `Checker`

Procedural checkerboard texture.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `beaming_exponent` | `number` | yes | Exponent used for relativistic beaming intensity scaling. |
| `color1` | `array` | yes | RGB tuple for first checker color. |
| `color2` | `array` | yes | RGB tuple for second checker color. |
| `color_normalization` | `CIETristimulusNormalization` | yes | Normalization applied to generated texture colors. |
| `height` | `number` | yes | Height of each checker cell in UV space. |
| `width` | `number` | yes | Width of each checker cell in UV space. |

### `BlackBody`

Color from black-body spectrum at emitter temperature.

| Field | Type | Required | Description |
| --- | --- | --- | --- |
| `beaming_exponent` | `number` | yes | Exponent used for relativistic beaming intensity scaling. |
| `color_normalization` | `CIETristimulusNormalization` | yes | Normalization applied to generated black-body colors. |


