//! The layers of the globe: one per field of the world, in the world's
//! order, coloured exactly as the field's PNG map in the world directory
//! (`kuniumi_io::world_dir::map_colours`): class fields with the categorical
//! palette, `surface_height` with the elevation map at the sea level over the
//! axis of `hypsometry.png`, other quantities with viridis over their range.

use kuniumi_io::Rgb8;
use kuniumi_io::colormap::Colormap;
use kuniumi_io::world_dir::map_colours;
use kuniumi_sim::metrics::{ElevationMetrics, HYPSOMETRY_BIN_M, Hypsometry, MAX_HYPSOMETRY_BINS};
use kuniumi_sim::stages::crust::SURFACE_HEIGHT;
use kuniumi_world::{FieldKind, World};

/// The colour of every cell of a world without fields.
pub const NO_FIELD: Rgb8 = [170, 170, 170];

/// One field as a layer: how its values become colours.
#[derive(Clone, Debug, PartialEq)]
pub struct Layer {
    /// The field's name.
    pub name: String,
    /// The field's SI unit (`""` if dimensionless or a class).
    pub unit: String,
    /// Quantity or class.
    pub kind: FieldKind,
    /// The colour map.
    pub colormap: Colormap,
    /// The value drawn with the first colour of the map, in the field's unit.
    pub min: f64,
    /// The value drawn with the last colour of the map, in the field's unit.
    pub max: f64,
}

impl Layer {
    /// The layer's legend in one line, for the terminal.
    pub fn legend(&self) -> String {
        let (name, min, max) = (&self.name, self.min, self.max);
        let unit = if self.unit.is_empty() {
            "dimensionless"
        } else {
            &self.unit
        };
        match self.colormap {
            Colormap::Categorical => {
                format!("{name} (class ids {min} to {max}): categorical palette")
            }
            Colormap::Viridis => format!("{name} ({unit}): viridis from {min} to {max}"),
            Colormap::Diverging { centre } => {
                format!("{name} ({unit}): diverging around {centre}, from {min} to {max}")
            }
            Colormap::Elevation { sea_level } => format!(
                "{name} ({unit}): elevation, blue below the sea level at {sea_level} and green to \
                 brown above, from {min} to {max}"
            ),
        }
    }
}

/// The layers of a world and the one shown.
#[derive(Clone, Debug, PartialEq)]
pub struct Layers {
    layers: Vec<Layer>,
    current: usize,
}

impl Layers {
    /// The layers of `world`, showing `surface_height` if it has one, else
    /// its first field.
    pub fn new(world: &World) -> Self {
        let metrics = metrics(world);
        let layers = world
            .fields
            .iter()
            .map(|field| {
                let (colormap, range) = map_colours(field, world, metrics.as_ref());
                let (min, max) = range.resolve(&field.values);
                Layer {
                    name: field.name.clone(),
                    unit: field.unit.clone(),
                    kind: field.kind,
                    colormap,
                    min,
                    max,
                }
            })
            .collect();
        let current = world
            .fields
            .iter()
            .position(|field| field.name == SURFACE_HEIGHT)
            .unwrap_or(0);
        Self { layers, current }
    }

    /// The number of layers (the world's fields).
    pub fn len(&self) -> usize {
        self.layers.len()
    }

    /// The index of the layer shown.
    pub fn current_index(&self) -> usize {
        self.current
    }

    /// The layer shown, or `None` for a world without fields.
    pub fn current(&self) -> Option<&Layer> {
        self.layers.get(self.current)
    }

    /// Shows the next layer, after the last the first.
    pub fn next(&mut self) {
        if !self.layers.is_empty() {
            self.current = (self.current + 1) % self.layers.len();
        }
    }

    /// Shows the previous layer, before the first the last.
    pub fn previous(&mut self) {
        if !self.layers.is_empty() {
            self.current = (self.current + self.layers.len() - 1) % self.layers.len();
        }
    }

    /// Shows the layer of the field `name`.
    ///
    /// # Errors
    ///
    /// A message listing the fields, if the world has no field `name`.
    pub fn select(&mut self, name: &str) -> Result<(), String> {
        match self.layers.iter().position(|layer| layer.name == name) {
            Some(index) => {
                self.current = index;
                Ok(())
            }
            None => {
                let names: Vec<&str> = self
                    .layers
                    .iter()
                    .map(|layer| layer.name.as_str())
                    .collect();
                Err(format!(
                    "the world has no field `{name}`; its fields: {}",
                    names.join(", ")
                ))
            }
        }
    }

    /// The colour of every cell of `world` in the layer shown ([`NO_FIELD`]
    /// for a world without fields).
    pub fn colours(&self, world: &World) -> Vec<Rgb8> {
        match self.current() {
            Some(layer) => {
                let field = &world.fields[self.current];
                debug_assert_eq!(field.name, layer.name, "the layers of another world");
                let map = layer.colormap;
                field
                    .values
                    .iter()
                    .map(|&v| map.colour(v, layer.min, layer.max))
                    .collect()
            }
            None => vec![NO_FIELD; world.grid.len()],
        }
    }
}

/// The validation metrics of `world` as `world_dir::write` computes them for
/// its maps, or `None` where `write` would refuse the world (a non-finite
/// height, an elevation above the sea level that overflows, too many
/// histogram bins) or the world has no `surface_height`. A world directory
/// written by kuniumi always passes these checks.
pub fn metrics(world: &World) -> Option<ElevationMetrics> {
    let heights = &world.field(SURFACE_HEIGHT)?.values;
    let sea_level = world.sea_level_m;
    let finite = heights.iter().all(|h| (h - sea_level).is_finite());
    (finite
        && Hypsometry::bins_needed(heights, sea_level, HYPSOMETRY_BIN_M)
            <= MAX_HYPSOMETRY_BINS as f64)
        .then(|| ElevationMetrics::of_world(world))
        .flatten()
}

#[cfg(test)]
mod tests {
    use kuniumi_io::colormap::Range;
    use kuniumi_world::{PlanetParams, Settings};

    use super::*;

    fn world() -> World {
        let settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        kuniumi_sim::generate(&PlanetParams::earth(), &settings, 3, &mut |_| {}).unwrap()
    }

    #[test]
    fn every_field_is_a_layer_in_the_colours_of_its_map() {
        let world = world();
        assert!(world.fields.len() >= 10);
        let expected_metrics = ElevationMetrics::of_world(&world);
        assert!(expected_metrics.is_some());
        let mut layers = Layers::new(&world);
        assert_eq!(layers.len(), world.fields.len());
        // surface_height first.
        assert_eq!(layers.current().unwrap().name, SURFACE_HEIGHT);
        let mut seen_class = false;
        for field in &world.fields {
            layers.select(&field.name).unwrap();
            let layer = layers.current().unwrap();
            assert_eq!(
                (&layer.name, &layer.unit, layer.kind),
                (&field.name, &field.unit, field.kind)
            );
            let (map, range) = map_colours(field, &world, expected_metrics.as_ref());
            assert_eq!(layer.colormap, map, "{}", field.name);
            assert_eq!((layer.min, layer.max), range.resolve(&field.values));
            assert_eq!(
                layers.colours(&world),
                map.colours(&field.values, range),
                "{}",
                field.name
            );
            match field.kind {
                FieldKind::Class => {
                    seen_class = true;
                    assert_eq!(layer.colormap, Colormap::Categorical);
                }
                FieldKind::Quantity if field.name == SURFACE_HEIGHT => {
                    assert_eq!(
                        layer.colormap,
                        Colormap::Elevation {
                            sea_level: world.sea_level_m
                        }
                    );
                    // The axis of hypsometry.png: −11 km to +9 km above the sea level.
                    assert_eq!(
                        layer.min.to_bits(),
                        (world.sea_level_m - 11_000.0).to_bits()
                    );
                    assert_eq!(layer.max.to_bits(), (world.sea_level_m + 9_000.0).to_bits());
                }
                FieldKind::Quantity => assert_eq!(layer.colormap, Colormap::Viridis),
            }
        }
        assert!(seen_class);
    }

    #[test]
    fn cycling_visits_every_layer_and_wraps() {
        let world = world();
        let mut layers = Layers::new(&world);
        let n = layers.len();
        let start = layers.current_index();
        let mut visited = vec![false; n];
        for _ in 0..n {
            visited[layers.current_index()] = true;
            layers.next();
        }
        assert!(visited.iter().all(|&v| v));
        assert_eq!(layers.current_index(), start);
        layers.select(&world.fields[0].name).unwrap();
        layers.previous();
        assert_eq!(layers.current_index(), n - 1);
        layers.next();
        assert_eq!(layers.current_index(), 0);
        let error = layers.select("rainfall").unwrap_err();
        assert!(
            error.starts_with("the world has no field `rainfall`; its fields: cell_area, "),
            "{error}"
        );
        assert_eq!(layers.current_index(), 0);
    }

    #[test]
    fn worlds_without_relief_or_fields() {
        let mut world = world();
        world.fields.retain(|f| f.name != SURFACE_HEIGHT);
        let layers = Layers::new(&world);
        assert_eq!(layers.current().unwrap().name, world.fields[0].name);
        assert!(metrics(&world).is_none());
        world.fields.clear();
        let mut empty = Layers::new(&world);
        assert_eq!((empty.len(), empty.current()), (0, None));
        empty.next();
        empty.previous();
        assert_eq!(empty.colours(&world), vec![NO_FIELD; world.grid.len()]);
    }

    #[test]
    fn metrics_follow_the_checks_of_world_dir() {
        let mut world = world();
        assert_eq!(metrics(&world), ElevationMetrics::of_world(&world));
        world.field_mut(SURFACE_HEIGHT).unwrap().values[5] = f64::NAN;
        assert!(metrics(&world).is_none());
        world.field_mut(SURFACE_HEIGHT).unwrap().values[5] = 1e300;
        assert!(metrics(&world).is_none(), "too many bins");
        // Without metrics the elevation map spans the field's own range.
        world.field_mut(SURFACE_HEIGHT).unwrap().values[5] = f64::NAN;
        let layers = Layers::new(&world);
        let field = world.field(SURFACE_HEIGHT).unwrap();
        assert_eq!(
            (layers.current().unwrap().min, layers.current().unwrap().max),
            Range::Auto.resolve(&field.values)
        );
    }

    #[test]
    fn legends() {
        let layer = |colormap, min, max, kind, unit: &str| Layer {
            name: "f".into(),
            unit: unit.into(),
            kind,
            colormap,
            min,
            max,
        };
        assert_eq!(
            layer(Colormap::Viridis, -0.5, 2.25, FieldKind::Quantity, "m/s").legend(),
            "f (m/s): viridis from -0.5 to 2.25"
        );
        assert_eq!(
            layer(Colormap::Viridis, 0.0, 1.0, FieldKind::Quantity, "").legend(),
            "f (dimensionless): viridis from 0 to 1"
        );
        assert_eq!(
            layer(
                Colormap::Elevation { sea_level: -700.0 },
                -11_700.0,
                8_300.0,
                FieldKind::Quantity,
                "m"
            )
            .legend(),
            "f (m): elevation, blue below the sea level at -700 and green to brown above, from -11700 to 8300"
        );
        assert_eq!(
            layer(Colormap::Categorical, 0.0, 11.0, FieldKind::Class, "").legend(),
            "f (class ids 0 to 11): categorical palette"
        );
        assert_eq!(
            layer(
                Colormap::Diverging { centre: 0.0 },
                -1.0,
                3.0,
                FieldKind::Quantity,
                "K"
            )
            .legend(),
            "f (K): diverging around 0, from -1 to 3"
        );
    }
}
