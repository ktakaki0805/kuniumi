//! What the browser's page hands to the viewer (ticket MV-03): the paths of
//! the files the user chose, of which the world's are read, and the GPU
//! backend asked for in the page's address. Compiled for wasm32, where
//! `web` uses it, and for the native tests.

use kuniumi_io::world_dir::{GRID_FILE, METADATA_FILE};

/// The world among the files the user chose.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WorldFiles {
    /// The directory that holds `world.json`, as the browser gives it (`/`
    /// between the parts); empty if the files were chosen without one.
    pub dir: String,
    /// The files of that directory: each name with its index in the list
    /// chosen, sorted by name.
    pub files: Vec<(String, usize)>,
}

impl WorldFiles {
    /// The index in the list chosen of the file `name` of the world's
    /// directory.
    pub fn index(&self, name: &str) -> Option<usize> {
        self.files
            .binary_search_by(|(file, _)| file.as_str().cmp(name))
            .ok()
            .map(|k| self.files[k].1)
    }

    /// The name of the world's directory (its last part), or `world` if the
    /// files were chosen without one.
    pub fn name(&self) -> &str {
        match self.dir.rsplit('/').next() {
            Some(name) if !name.is_empty() => name,
            _ => "world",
        }
    }

    /// The indices of the files that reading the world may ask for, in the
    /// order of the list chosen: `world.json`, [`GRID_FILE`] and every
    /// `.kfield` file of the directory (the maps and the metrics are not
    /// read).
    pub fn wanted(&self) -> Vec<usize> {
        let mut wanted: Vec<usize> = self
            .files
            .iter()
            .filter(|(name, _)| {
                name == METADATA_FILE || name == GRID_FILE || name.ends_with(".kfield")
            })
            .map(|&(_, index)| index)
            .collect();
        wanted.sort_unstable();
        wanted
    }
}

/// `path` split into its directory (empty for a bare name) and its name.
fn split(path: &str) -> (&str, &str) {
    path.rsplit_once('/').unwrap_or(("", path))
}

/// The world among the files at `paths` (relative paths with `/`, as a
/// browser gives `webkitRelativePath`, or bare names): the one directory
/// that holds a `world.json`, and its files; files of other directories
/// are left out.
///
/// # Errors
///
/// If no path names a `world.json`, if several do (in different
/// directories: a directory holding several worlds was chosen), or if a
/// path is given twice.
pub fn world_files(paths: &[String]) -> Result<WorldFiles, String> {
    let mut sorted: Vec<&str> = paths.iter().map(String::as_str).collect();
    sorted.sort_unstable();
    if let Some(pair) = sorted.windows(2).find(|pair| pair[0] == pair[1]) {
        return Err(format!("{} was chosen twice", pair[0]));
    }
    let mut worlds: Vec<&str> = paths
        .iter()
        .map(|path| split(path))
        .filter(|&(_, name)| name == METADATA_FILE)
        .map(|(dir, _)| dir)
        .collect();
    worlds.sort_unstable();
    let dir = match worlds[..] {
        [] => {
            return Err(format!(
                "no {METADATA_FILE} among the {} files chosen: choose a world directory written \
                 by `kuniumi generate`",
                paths.len()
            ));
        }
        [dir] => dir,
        _ => {
            let shown: Vec<&str> = worlds
                .iter()
                .map(|dir| {
                    if dir.is_empty() {
                        "the files themselves"
                    } else {
                        dir
                    }
                })
                .collect();
            let (last, rest) = shown.split_last().expect("several");
            return Err(format!(
                "the files chosen hold {} worlds, in {} and {last}: choose one world directory",
                worlds.len(),
                rest.join(", ")
            ));
        }
    };
    let mut files: Vec<(String, usize)> = paths
        .iter()
        .enumerate()
        .filter_map(|(index, path)| {
            let (parent, name) = split(path);
            (parent == dir).then(|| (name.to_owned(), index))
        })
        .collect();
    files.sort_unstable();
    Ok(WorldFiles {
        dir: dir.to_owned(),
        files,
    })
}

/// The GPU backends for the backend `name` asked for in the page's address
/// (`?backend=…`): none, both (WebGPU where the browser has an adapter, else
/// WebGL2); `webgpu` or `webgl`, that one only.
///
/// # Errors
///
/// If `name` is none of these.
pub fn backends(name: &str) -> Result<wgpu::Backends, String> {
    match name {
        "" => Ok(wgpu::Backends::BROWSER_WEBGPU | wgpu::Backends::GL),
        "webgpu" => Ok(wgpu::Backends::BROWSER_WEBGPU),
        "webgl" => Ok(wgpu::Backends::GL),
        _ => Err(format!(
            "unknown backend `{name}` in the page's address: use backend=webgpu or backend=webgl"
        )),
    }
}

/// The times of loading a world in the page and the module's memory
/// afterwards, for the page's status line.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LoadTimes {
    /// From the user's choice until the viewer has the files: the browser
    /// reading them (ms).
    pub files_ms: f64,
    /// Reading the world from the files (ms).
    pub read_ms: f64,
    /// The scene: the mesh, the layers, the picker (ms).
    pub scene_ms: f64,
    /// The GPU: the adapter, the device, the buffers (ms).
    pub gpu_ms: f64,
    /// From the user's choice to the first frame presented, once it has
    /// been (ms).
    pub first_frame_ms: Option<f64>,
    /// The size of the module's memory after loading (bytes).
    pub memory_bytes: u64,
}

impl LoadTimes {
    /// The status line of the world `name` of `cells` cells and `fields`
    /// fields, drawn by `gpu` (backend and adapter).
    pub fn status(&self, name: &str, cells: usize, fields: usize, gpu: &str) -> String {
        let s = |ms: f64| format!("{:.2} s", ms / 1e3);
        let first = self.first_frame_ms.map_or(String::new(), |ms| {
            format!(", first frame {} after the choice", s(ms))
        });
        format!(
            "{name}: {cells} cells, {fields} fields; {gpu}; files {}, world {}, scene {}, GPU {}\
             {first}; memory {} MB",
            s(self.files_ms),
            s(self.read_ms),
            s(self.scene_ms),
            s(self.gpu_ms),
            self.memory_bytes / 1_000_000
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(paths: &[&str]) -> Vec<String> {
        paths.iter().map(|p| p.to_string()).collect()
    }

    #[test]
    fn a_picked_directory() {
        let paths = strings(&[
            "650k/plate_id.png",
            "650k/world.json",
            "650k/grid.kgrid",
            "650k/plate_id.kfield",
            "650k/metrics.txt",
            "650k/cell_area.kfield",
            "650k/hypsometry.png",
        ]);
        let world = world_files(&paths).unwrap();
        assert_eq!(world.dir, "650k");
        assert_eq!(world.name(), "650k");
        assert_eq!(
            world.files,
            [
                ("cell_area.kfield".to_string(), 5),
                ("grid.kgrid".to_string(), 2),
                ("hypsometry.png".to_string(), 6),
                ("metrics.txt".to_string(), 4),
                ("plate_id.kfield".to_string(), 3),
                ("plate_id.png".to_string(), 0),
                ("world.json".to_string(), 1),
            ]
        );
        assert_eq!(world.index("world.json"), Some(1));
        assert_eq!(world.index("plate_id.kfield"), Some(3));
        assert_eq!(world.index("rain.kfield"), None);
        assert_eq!(world.wanted(), [1, 2, 3, 5]);
    }

    #[test]
    fn a_nested_pick_and_bare_names() {
        // The parent of the world picked: a sibling directory without a
        // world, and files of the parent itself, are left out.
        let paths = strings(&[
            "MV-03/notes/world.txt",
            "MV-03/notes/grid.kgrid",
            "MV-03/650k/grid.kgrid",
            "MV-03/grid.kgrid",
            "MV-03/650k/world.json",
            "MV-03/650k/sub/x.kfield",
            "MV-03/650k/surface_height.kfield",
        ]);
        let world = world_files(&paths).unwrap();
        assert_eq!(world.dir, "MV-03/650k");
        assert_eq!(world.name(), "650k");
        assert_eq!(
            world.files,
            [
                ("grid.kgrid".to_string(), 2),
                ("surface_height.kfield".to_string(), 6),
                ("world.json".to_string(), 4),
            ]
        );
        assert_eq!(world.wanted(), [2, 4, 6]);
        // Files dropped without their directory.
        let bare = strings(&["grid.kgrid", "world.json", "a.kfield", "a.png"]);
        let world = world_files(&bare).unwrap();
        assert_eq!(world.dir, "");
        assert_eq!(world.name(), "world");
        assert_eq!(world.index("a.kfield"), Some(2));
        assert_eq!(world.wanted(), [0, 1, 2]);
    }

    #[test]
    fn no_world_several_worlds_and_repeats() {
        let none = strings(&["650k/grid.kgrid", "650k/plate_id.png"]);
        assert_eq!(
            world_files(&none).unwrap_err(),
            "no world.json among the 2 files chosen: choose a world directory written by \
             `kuniumi generate`"
        );
        assert_eq!(
            world_files(&[]).unwrap_err(),
            "no world.json among the 0 files chosen: choose a world directory written by \
             `kuniumi generate`"
        );
        let several = strings(&[
            "out/b/world.json",
            "out/a/world.json",
            "out/a/grid.kgrid",
            "out/world.json",
        ]);
        assert_eq!(
            world_files(&several).unwrap_err(),
            "the files chosen hold 3 worlds, in out, out/a and out/b: choose one world directory"
        );
        let twice = strings(&["w/world.json", "w/grid.kgrid", "w/grid.kgrid"]);
        assert_eq!(
            world_files(&twice).unwrap_err(),
            "w/grid.kgrid was chosen twice"
        );
    }

    #[test]
    fn backends_by_name() {
        use wgpu::Backends;
        assert_eq!(backends(""), Ok(Backends::BROWSER_WEBGPU | Backends::GL));
        assert_eq!(backends("webgpu"), Ok(Backends::BROWSER_WEBGPU));
        assert_eq!(backends("webgl"), Ok(Backends::GL));
        assert_eq!(
            backends("vulkan").unwrap_err(),
            "unknown backend `vulkan` in the page's address: use backend=webgpu or backend=webgl"
        );
    }

    #[test]
    fn status_by_hand() {
        let mut times = LoadTimes {
            files_ms: 850.0,
            read_ms: 1_204.0,
            scene_ms: 396.0,
            gpu_ms: 5.0,
            first_frame_ms: None,
            memory_bytes: 612_499_999,
        };
        assert_eq!(
            times.status("650k", 650_000, 13, "BrowserWebGpu on GPU"),
            "650k: 650000 cells, 13 fields; BrowserWebGpu on GPU; files 0.85 s, world 1.20 s, \
             scene 0.40 s, GPU 0.01 s; memory 612 MB"
        );
        times.first_frame_ms = Some(2_996.0);
        assert_eq!(
            times.status("650k", 650_000, 13, "Gl on GPU"),
            "650k: 650000 cells, 13 fields; Gl on GPU; files 0.85 s, world 1.20 s, scene 0.40 s, \
             GPU 0.01 s, first frame 3.00 s after the choice; memory 612 MB"
        );
    }
}
