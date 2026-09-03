# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Library Is

IPC Toolkit is a C++17 library for Incremental Potential Contact (IPC) — collision detection, distance computation, and contact potentials for physics simulation. It is **not** a full simulator; it provides reusable building blocks (barrier functions, CCD, distance queries, friction/adhesion potentials) for integration into external simulators. Python bindings are available via pybind11 as the `ipctk` package.

## Build Commands

All dependencies are automatically downloaded via CMake/CPM — no manual installation needed.

```bash
# Configure and build with tests (typical dev workflow)
cmake --preset=test
cmake --build --preset=test -j 8

# Release build
cmake --preset=release
cmake --build --preset=release -j 8

# Python bindings
cmake --preset=python
cmake --build --preset=python
# or: pip install .
```

## Running Tests

Tests use Catch2. Test files mirror the source structure: `tests/src/tests/` mirrors `src/ipc/`.

```bash
# Run all tests
cd build/test && ctest --verbose

# Run tests matching a pattern
cd build/test && ctest --verbose -R "test_name_pattern"

# Run test binary directly with Catch2 tag filter
./build/test/ipc_toolkit_tests "[tag]"
```

## Code Style

- **Formatter:** clang-format, WebKit-based style, **80-character column limit**. Pre-commit hooks enforce this — always run before committing:
  ```bash
  clang-format --style=file -i <file>
  ```
- **Linter:** clang-tidy with modernize-\*, performance-\*, readability-\*, bugprone-\* checks (treated as errors).
- **Naming:** Classes `CamelCase`, functions/methods `lower_case`, constants `UPPER_CASE`, member variables `lower_case`.

## Architecture

All library code lives under `src/ipc/` in the `ipc::` namespace.

### Core Types

- **`CollisionMesh`** (`collision_mesh.hpp`) — Central immutable data structure mapping a volumetric FE mesh to surface collision vertices/edges/faces. Nearly everything takes a `CollisionMesh` reference.
- **`ipc.hpp`** — Top-level API: `is_step_collision_free()`, `compute_collision_free_stepsize()`, `has_intersections()`.

### Collision Detection Pipeline

Data flows through three sequential stages:

1. **Broad phase** (`broad_phase/`) — Spatial pruning via AABB-based methods (LBVH, Sweep-and-Prune, Hash Grid, Brute Force; optionally CUDA-accelerated STQ).
2. **Candidates** (`candidates/`) — Typed potential collision pairs (vertex-vertex, edge-vertex, edge-edge, face-vertex, edge-face, face-face, plane-vertex) produced by the broad phase.
3. **CCD** (`ccd/`) — Narrow-phase continuous collision detection (Tight Inclusion, Additive CCD, Nonlinear CCD, Inexact CCD) returning exact collision times or boolean results.

### Potentials (`potentials/`)

Each potential provides `operator()`, `gradient()`, and `hessian()`:
- **`BarrierPotential`** — Main IPC contact barrier, parameterized by activation distance and stiffness.
- **`FrictionPotential`** — Lagged dissipative friction.
- **`NormalAdhesionPotential`** / **`TangentialAdhesionPotential`** — Surface adhesion variants.

### Supporting Modules

- `distance/` — Geometric distance functions (point-point, point-edge, edge-edge, point-triangle, etc.) with signed distance variants in `distance/signed/`.
- `tangent/` — Tangent space basis computation for friction/adhesion.
- `barrier/` — Barrier function math and derivatives.
- `collisions/` — Collision storage, split into `normal/` and `tangential/` subdirectories.

### Python Bindings (`python/`)

pybind11 bindings in `python/src/` with entry point `bindings.cpp` calling 50+ `define_*` functions that mirror the C++ API.

## Key CMake Options

| Option | Default | Purpose |
|--------|---------|---------|
| `IPC_TOOLKIT_BUILD_TESTS` | ON (top-level) | Build Catch2 tests |
| `IPC_TOOLKIT_BUILD_PYTHON` | OFF | Build pybind11 bindings |
| `IPC_TOOLKIT_WITH_CUDA` | OFF | GPU-accelerated CCD |
| `IPC_TOOLKIT_WITH_SIMD` | ON | SIMD optimizations |
| `IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION` | OFF | Exact rational arithmetic |

## Plotting Style

- Every benchmark figure (in a PR, doc, or paper artifact) is built with **Plotly** (`plotly.express`/`plotly.graph_objects`, not matplotlib).
- **Font** — leave the template's default sans-serif (Arial-family), regular weight throughout; don't switch families or italicize. Font size: 16pt for the figure title, 12pt for regular text (e.g., labels), 11pt for sublabels, 14pt for subtitles
- **Figure title** — one figure-level title via `fig.update_layout(title=dict(text=..., x=0.5, xanchor='center'))`, centered above the plot, naming the benchmark and what's measured
  - Include a subtitle on a second line noting whether smaller or bigger is better — the second line is smaller and greyed.
- **Footnote** — one small grey line below the axis stating methodology: how a path/config was selected per case, the number of runs, and the aggregation (e.g. "Median of three runs").

### Specific Styling

Which template applies depends on what is being compared:
- **Methods/algorithms, or the components of one algorithm** — follow the [Method Comparisons](#method-omparisons) template, the canonical reference for those; color encoding and layout below.
- **Before vs. after** — follow the [Before/After Comparisons](#beforeafter-comparisons); it is a different figure — horizontal, single-panel, no in-bar labels, with the speedup annotations colored by outcome.

#### Color encoding
Pick the palette to match what's being compared, and never reuse one plot's color meaning for another:
- **Methods/algorithms** (e.g. broad-phase methods) — `px.colors.qualitative.Set2`, one fixed `Set2` color per method; give the method the PR introduces a contrasting marker/line outline color (e.g. `crimson`, via `marker_line_color`/`marker_line_width`) so it doesn't just blend in.
- **Components of one algorithm** (e.g. pipeline stages in a stacked bar) — `px.colors.qualitative.T10` (Tableau10), taken in order, one color per component, consistent across every panel/scene so a color always means the same stage.
- **Before vs. after** — use specifically `#86b6ef` for before and `#2a78d6` for after; annotate each pair with its speedup factor rather than relying on the reader to eyeball bar lengths. See [Before/After Comparisons](#beforeafter-comparisons) for the full layout.

#### Method Comparisons

- **Template** — white background (`plot_bgcolor='white'`, `paper_bgcolor='white'`) with manual light-grey gridlines (`fig.update_yaxes(showgrid=True, gridcolor='#E0E0E0')`).
- **Subplot titles** — build panels with `make_subplots(subplot_titles=[...])` so each subplot gets its own smaller centered title naming the scene and its size (e.g. "Cloth-Funnel (Faces: 18.5K)") — never a shared/omitted per-panel title.
- **Y-axis title** — only the leftmost subplot in a row carries one (e.g. `fig.update_yaxes(title_text="Time (ms)", row=1, col=1)`); don't repeat it on every panel.
- **In-bar annotations** — each bar's raw numeric value via `text=`, `textposition='outside'`, horizontal (no `textangle`), just above the bar.
- **Legend** — one shared legend for the whole figure (`fig.update_layout(legend=dict(orientation='h', yanchor='bottom', y=-0.25, xanchor='center', x=0.5))`), centered below all subplots, laid out horizontally — never a separate legend per subplot.
- **Two-method improvement %** — when the comparison is exactly two methods, additionally annotate the faster method's percentage improvement via `fig.add_annotation(...)`, placed inside that bar near its top (e.g. `y = value * 0.95`), small font (~9pt), grey (`#333333`), no arrow.

#### Before/After Comparisons

A before/after figure measures one quantity across many cases, so it is **horizontal and single-panel** — not a grid of vertical-bar subplots. The rules above for panels, figure titles, and in-bar text do **not** apply here; these replace them.

- **Template** — `template='simple_white'` (Plotly's built-in template)
- **Bars** — horizontal grouped bars, one pair per case, `before` above `after` within the pair; cases run top to bottom in a meaningful order (ascending problem size, or grouped by family) — never alphabetical by accident.
- **X-axis** — log scale whenever the cases span more than ~20x, which they usually do; the axis title names the quantity, the unit, and the scale (e.g. `"gradient assembly, collision DOF (ms, log)"`).
- **In-bar annotations** — none, and no vertical text; the bars carry the shape.
- **Speedup annotation** — placed to the right of each pair, just past the longer bar, bold, and colored by outcome so the result reads while skimming, without comparing bar lengths:
  - **Green (`#2E7D32`)** — a real win, ratio ≥ 1.05x.
  - **Dark grey (`#444`)** — a wash, ratio 0.95x–1.04x; neutral by design — a 1.02x is not a result, and coloring it green oversells the change.
  - **Crimson** — a regression, ratio < 0.95x; never suppress these — a red row the PR explains is worth more than a chart that hides it.
  - **Format** — `1.43x`, two decimals, bucketed on the rounded value so the color always agrees with the printed number — a 0.997x must not render as a red `1.00x`, and a 1.0499x must not render green.
- **Category labels** — two lines, right-aligned against the axis: case name on the first line, its size and any variant/config on the second, separated by a middle dot (e.g. `"cloth-ball"` over `"3,012 · two-pass"`); the second line is smaller and greyed.
- **Legend** — inside the plot, top right, vertical, no frame (`legend=dict(x=0.99, y=0.99, xanchor='right', yanchor='top', bgcolor='rgba(0,0,0,0)')`), no legend title; entries name themselves (`"before (main)"`, `"after"`).

## Writing PR Descriptions

Start from `.github/pull_request_template.md` — keep its section skeleton (`# Description`, `## Type of change`, `# How Has This Been Tested?`, `# Checklist`), but write each section of the project's best PRs (e.g. [#246](https://github.com/ipc-sim/ipc-toolkit/pull/246), [#248](https://github.com/ipc-sim/ipc-toolkit/pull/248), [#230](https://github.com/ipc-sim/ipc-toolkit/pull/230), [#228](https://github.com/ipc-sim/ipc-toolkit/pull/228), [#222](https://github.com/ipc-sim/ipc-toolkit/pull/222), [#210](https://github.com/ipc-sim/ipc-toolkit/pull/210), [#207](https://github.com/ipc-sim/ipc-toolkit/pull/207), [#192](https://github.com/ipc-sim/ipc-toolkit/pull/192), [#188](https://github.com/ipc-sim/ipc-toolkit/pull/188)) instead of a one-line fill-in-the-blank.

- **Concise style** — write in bullet points, not paragraphs; one short sentence per bullet.
- **Title** — short, specific, states the concrete change. "Improve ..." for improving existing functionality. "Fix ..." for bug fixes.
- **Length** — a PR description is not a changelog. Draft it, then cut it roughly in half. Anything exhaustive (the full API surface, every bug fixed in passing, every deprecation) do not belong in the PR description; the PR carries only what a reviewer needs to review *this diff*.
- **`# Description`**
  - Start the description with the following:
    1. One sentence naming what was restructured and which areas it touches.
    2. A short list of the mechanism (the layers, the stages, the new types).
    3. A list of impacts this PR will have on the codebase and downstream users.
  - Attach each number to the capability it justifies, e.g. "Fixed-size optimization of shared 2D and 3D code. E.g., `point_edge_closest_point` at 8.7 ns on a matrix row vs 4.7 ns on a `Vector3d` — the same arithmetic, 1.9× apart, purely from losing the compile-time size."
  - **Subsections**
    - **API changes** — call out explicitly: new/changed function signatures, new flags and their default/fallback behavior, new public types, backward compatibility.
    - **Performance/behavior claims** — back with numbers (timings, percentages, speedup factors), never "faster"/"more efficient". Embed a plot or screenshot when one supports the point, with a bullet right after it saying what to look for (what's absent, what moved, why). When a benchmark spans multiple scenes, use a results table (scene, size, before, after, speedup) rather than prose or a single average. Any generated benchmark plot follows [Benchmark Plot Style](#benchmark-plot-style).
    - **Breaking changes** — give them their own visually distinct block (a subheader or bolded lead-in) separate from ordinary enhancements, and list deprecations explicitly rather than folding them into a general changes list.
    - **New dependencies** — list with license, fetch mechanism (e.g. CPM), and why they're safe to add (small/static footprint, no transitive bloat).
    - **References** — when the change implements or is based on a paper, cite it by author/year with a link (e.g. "Apetrei [2014]", "Chen et al. [2025]") right where the algorithm is introduced, not just in a bibliography.
    - A small set of emoji section markers (🚀 New Features, ⚡ Enhancements, 💥 Breaking Changes, 🐛 Bug Fixes) is acceptable for grouping a large multi-part PR — use them as `##`/`###` headers, not scattered inline, and only when the PR actually spans multiple such categories.
  - Don't rely on GitHub-generated diff-hunk permalinks (`[[1]](diffhunk://...)`) as the explanation — they render as noise off-GitHub; say what changed in the sentence itself.
  - Keep "Fixes # (issue)" only when an issue exists; delete the line otherwise.
- **Don't narrate the work** — no "three findings that cost time to get wrong", no "I probed", no "worth its own issue", no defending a decision already made. The reviewer is reading the diff, not the journey. State the result and stop.
- **Explain a speed-up by its hardware mechanism** where one exists — "rows of a column-major vertex matrix are now copied to the local stack for cache efficiency" tells a reviewer more than "recovers the compile-time dimension".
- **Use the project's own vocabulary** — these are `functions`, not "kernels". Never introduce a private term for something the codebase already names.
- **After each plot, one or two bullets of conclusion** — give the aggregate, then say *which class of thing moved*: "Dimension-agnostic functions benefit the most." Never restate the axes or narrate the picture.
- **`## Type of change`** — actually check the boxes that apply and delete the rest, per the template's own instructions.
- **`# How Has This Been Tested?`** — replace "Test A"/"Test B" placeholders with the real tests run, and when benchmarks are involved, state methodology here: hardware, compiler, build type, number of runs, observed variance.
- **`# Checklist`** — check only boxes you've actually verified; leave the rest unchecked for the reviewer rather than checking them by default.
- Never append a "🤖 Generated with Claude Code" (or similar AI-attribution) footer to the PR description, regardless of the [AI Co-Authorship](#ai-co-authorship-on-commits) rule for commits — that rule governs commit trailers only, not PR descriptions.

## PR Checklist

1. Code passes `clang-format` and introduces no new `clang-tidy` warnings.
2. Complex logic is commented.
3. Public API changes are reflected in documentation.
4. Unit tests added/updated (75% patch coverage target).

## AI Co-Authorship on Commits

When contributing to this project, AI co-authors may be included in the commit history to acknowledge their assistance in code generation, review, and optimization.
However, AI co-authors should not be listed for trivial changes such as formatting, renaming, or minor refactoring. The following guidelines help determine when to include an AI co-author:
- **Include an AI co-author** if the AI contributed to:
  - Adding new features or functionality.
  - Writing new functions or classes.
  - Refactoring or optimizing existing code.
- **Do not include an AI co-author** for:
  - Formatting changes (e.g., running clang-format).
  - Renaming variables or functions.
  - Minor refactoring (e.g., reordering code, adding comments).
  - Writing or improving unit tests.