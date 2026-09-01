#import "@preview/cetz:0.5.2"

// Export at print resolution with, e.g.:
//   typst compile --format png --ppi 300 fig1.typ fig1.png

#set page(width: auto, height: auto, margin: 1.5cm, fill: white)
#set text(size: 10.8pt)

#let INK = rgb("#111111")
#let ORANGE = rgb("#d98c3d")
#let BLUE = rgb("#4a90d9")
#let GREY = rgb("#777777")

#cetz.canvas({
  import cetz.draw: *

  // ---------------------------------------------------------------- helpers

  // Boxes are sized with generous headroom below the title/subtitle so a
  // PNG panel can be dropped in later without resizing anything.
  let node-box(x, y, w, h, title: none, subtitle: none, image-path: none, subtitle-color: GREY) = {
    rect((x, y), (x + w, y + h), stroke: INK + 1.4pt, radius: 4pt, fill: white)
    if title != none {
      content((x + w / 2, y + h - 0.4), text(size: 14.4pt, weight: "bold", fill: INK)[#title])
    }
    if subtitle != none {
      let sy = if title != none { y + h - 0.85 } else { y + h - 0.4 }
      content((x + w / 2, sy), text(size: 10.2pt, style: "italic", fill: subtitle-color)[#subtitle])
    }
    if image-path != none {
      let header = if title != none and subtitle != none { 1.15 }
        else if title != none or subtitle != none { 0.7 }
        else { 0 }
      let pad = 0.15
      let aw = w - 2 * pad
      let ah = h - header - pad
      content((x + w / 2, y + pad + ah / 2),
        image(image-path, width: aw * 1cm, height: ah * 1cm, fit: "contain"))
    }
  }

  // Box whose body is a short bold list of array/field names (e.g. what's
  // actually stored inside an npz), instead of a title + italic subtitle.
  // Wrapped to the box width so a long item can break onto its own line.
  let node-box-list(x, y, w, h, body) = {
    rect((x, y), (x + w, y + h), stroke: INK + 1.4pt, radius: 4pt, fill: white)
    content((x + w / 2, y + h / 2),
      box(width: (w - 0.2) * 1cm,
        align(center, text(size: 10.2pt, fill: INK)[#body])))
  }

  // Same as node-box, but drawn as a shallow stack of cards to signal
  // "there are usually several of these files".
  let node-box-stack(x, y, w, h, title: none, subtitle: none, image-path: none) = {
    let off = 0.16
    rect((x + 2 * off, y + 2 * off), (x + 2 * off + w, y + 2 * off + h), stroke: INK + 1.4pt, radius: 4pt, fill: white)
    rect((x + off, y + off), (x + off + w, y + off + h), stroke: INK + 1.4pt, radius: 4pt, fill: white)
    node-box(x, y, w, h, title: title, subtitle: subtitle, image-path: image-path)
  }

  let node-ellipse(cx, cy, w, h, txt, color) = {
    circle((cx, cy), radius: (w / 2, h / 2), stroke: color + 1.6pt, fill: white)
    content((cx, cy), align(center, text(size: 12pt, weight: "bold", fill: color)[#txt]))
  }

  // Rounded-rectangle process node (used for the "visualisation" tools).
  let node-rbox(x, y, w, h, txt, color) = {
    rect((x, y), (x + w, y + h), stroke: color + 1.6pt, radius: 8pt, fill: white)
    content((x + w / 2, y + h / 2), align(center, text(size: 12pt, weight: "bold", fill: color)[#txt]))
  }

  // Height is fixed (not content-driven) so callers can compute exact
  // top/bottom edges for arrows instead of guessing an offset from cy.
  let PILL_H = 0.78
  let node-pill(cx, cy, txt, color) = {
    content((cx, cy), box(
      height: PILL_H * 1cm, stroke: color + 1.2pt, radius: 12pt, inset: (x: 9pt),
      align(horizon, text(size: 10.8pt, fill: color)[#txt]),
    ))
  }

  // Plain centered label, no outline -- used for annotations that aren't
  // themselves a pipeline node (e.g. "this input is required here").
  let node-tag(cx, cy, txt, color) = {
    content((cx, cy), align(center, text(size: 10.8pt, fill: color)[#txt]))
  }

  let bracket(x1, x2, y, txt, flip: false) = {
    let t = if flip { -0.2 } else { 0.2 }
    line((x1, y), (x1, y + t), (x2, y + t), (x2, y),
      stroke: (paint: INK, thickness: 1.2pt, join: "round"))
    content(((x1 + x2) / 2, y + (if flip { -0.55 } else { 0.55 })),
      text(size: 16.8pt, fill: INK)[#txt])
  }

  // Vertical curly accolade grouping a span of rows on the right edge --
  // a real stretched "}" glyph (math.stretch, sized to an absolute length
  // via a zero-ratio relative) rather than a hand-drawn approximation.
  let vbracket(x, y1, y2, txt, color) = {
    let h = y2 - y1
    let ymid = (y1 + y2) / 2
    content((x, ymid),
      text(fill: color)[$ stretch(brace.r, size: #(h * 1cm)) $])
    content((x + 0.55, ymid), anchor: "west",
      text(size: 15.6pt, fill: color)[#txt])
  }

  // plain connector (no arrowhead) -- used for "merge" segments of a fan-in/fan-out
  let wire(p1, p2, color: INK) = {
    line(p1, p2, stroke: (paint: color, thickness: 1.2pt, join: "round"))
  }

  // Break a segment into a 45-degree piece followed by a horizontal/vertical
  // piece -- diagonal first, straight last, so every arrowhead approaches
  // its target head-on instead of at an angle (which is what read as a
  // "nick" where the bend landed right next to the mark).
  let elbow(x1, y1, x2, y2) = {
    let dx = x2 - x1
    let dy = y2 - y1
    let adx = calc.abs(dx)
    let ady = calc.abs(dy)
    if adx == 0 or ady == 0 or adx == ady {
      ((x1, y1), (x2, y2))
    } else {
      let sx = if dx > 0 { 1 } else { -1 }
      let sy = if dy > 0 { 1 } else { -1 }
      let diag = calc.min(adx, ady)
      ((x1, y1), (x1 + sx * diag, y1 + sy * diag), (x2, y2))
    }
  }

  // Arrows stop a hair short of their nominal target so the mark doesn't
  // visually fuse into the border it's pointing at, and corners are
  // rounded so short elbow bends don't render as a sharp notch.
  let arrow(p1, p2, color: INK, straight: false) = {
    let pts = if straight { (p1, p2) } else {
      elbow(p1.at(0), p1.at(1), p2.at(0), p2.at(1))
    }
    let n = pts.len()
    let a = pts.at(n - 2)
    let b = pts.at(n - 1)
    let dx = b.at(0) - a.at(0)
    let dy = b.at(1) - a.at(1)
    let len = calc.sqrt(dx * dx + dy * dy)
    let gap = 0.06
    let tip = if len > gap {
      (b.at(0) - dx / len * gap, b.at(1) - dy / len * gap)
    } else {
      b
    }
    let final-pts = pts.slice(0, n - 1) + (tip,)
    line(..final-pts, stroke: (paint: color, thickness: 1.4pt, join: "round"), mark: (end: ">"))
  }

  // ------------------------------------------------------ top section brackets

  bracket(0.1, 12.9, 17.3, [Filtering])
  bracket(13.0, 22.5, 17.3, [Normalization])

  // ---- Row 1: BAM/CRAM -> alignmentSieve -> BAM/CRAM filtered -> bamCoverage/Compare -> bigWig

  let bam = (0.1, 13.1, 3.6, 3.4)
  node-box-stack(..bam, title: [BAM/CRAM], subtitle: [aligned reads], image-path: "subfigs/bamfile.png")

  node-ellipse(6.37, 14.8, 3.6, 1.38, [alignmentSieve], ORANGE)

  let bam-f = (8.72, 13.1, 3.6, 3.4)
  node-box-stack(..bam-f, title: [BAM/CRAM], subtitle: [filtered reads], image-path: "subfigs/bamfile_sieve.png")

  node-ellipse(15.17, 14.8, 3.96, 2.28, [bamCoverage \ or \ bamCompare], ORANGE)

  let bigwig = (17.70, 13.1, 4.2, 3.4)
  node-box-stack(..bigwig, title: [bigWig/bedGraph], subtitle: [normalized signals or ratios], image-path: "subfigs/bigwig.png")

  arrow((bam.at(0) + bam.at(2), bam.at(1) + bam.at(3) / 2), (4.57, 14.8))
  arrow((8.17, 14.8), (bam-f.at(0), bam-f.at(1) + bam-f.at(3) / 2))
  arrow((bam-f.at(0) + bam-f.at(2), bam-f.at(1) + bam-f.at(3) / 2), (13.19, 14.8))
  arrow((17.15, 14.8), (bigwig.at(0), bigwig.at(1) + bigwig.at(3) / 2))

  // ---- Row 2: multiBamSummary / computeMatrix

  node-ellipse(10.52, 11.35, 4.32, 1.44, [multiBamSummary], ORANGE)

  node-ellipse(19.80, 11.35, 3.6, 1.44, [computeMatrix], ORANGE)

  arrow((bam-f.at(0) + bam-f.at(2) / 2, bam-f.at(1)), (10.52, 12.07))
  arrow((bigwig.at(0) + bigwig.at(2) / 2, bigwig.at(1)), (19.80, 12.07))

  // computeMatrix always needs a regions file, regardless of mode
  node-tag(15.17, 11.35, [regions], GREY)
  arrow((15.82, 11.35), (18.0, 11.35))

  // ---- Row 2.5: output-mode pills for both matrix builders (same pattern)
  // pill pairs are centered under their parent ellipse (multiBamSummary / computeMatrix)

  // only the BED-file mode of multiBamSummary needs a regions file -- bins does not
  node-tag(6.52, 8.85, [regions], GREY)
  arrow((7.22, 8.85), (8.22, 8.85))

  node-pill(9.12, 8.85, [BED-file], ORANGE)
  node-pill(11.92, 8.85, [bins], ORANGE)
  wire((10.52, 10.63), (10.52, 9.935))
  wire((9.12, 9.935), (11.92, 9.935))
  arrow((9.12, 9.935), (9.12, 9.24))
  arrow((11.92, 9.935), (11.92, 9.24))

  node-pill(18.4, 8.85, [scale-regions], ORANGE)
  node-pill(21.2, 8.85, [reference-point], ORANGE)
  wire((19.80, 10.63), (19.80, 9.935))
  wire((18.4, 9.935), (21.2, 9.935))
  arrow((18.4, 9.935), (18.4, 9.24))
  arrow((21.2, 9.935), (21.2, 9.24))

  // ---- Row 3: matrices (npz), one per mode

  let npz-bedfile = (7.82, 4.75, 2.6, 2.7)
  node-box-list(..npz-bedfile, [npz \ counts \ scaling factors \ BED-file])

  let npz-bins = (10.62, 4.75, 2.6, 2.7)
  node-box-list(..npz-bins, [npz \ counts \ scaling factors])

  let npz-scale = (17.1, 4.75, 2.6, 2.7)
  node-box-list(..npz-scale, [matrix.gz])

  let npz-ref = (19.9, 4.75, 2.6, 2.7)
  node-box-list(..npz-ref, [matrix.gz])

  arrow((9.12, 8.46), (9.12, 7.45))
  arrow((11.92, 8.46), (11.92, 7.45))
  arrow((18.4, 8.46), (18.4, 7.45))
  arrow((21.2, 8.46), (21.2, 7.45))

  // ---- Row 3 -> 4: both multiBamSummary matrices feed all three QC tools

  wire((9.12, 4.75), (9.12, 4.05))
  wire((11.92, 4.75), (11.92, 4.05))
  wire((6.02, 4.05), (15.02, 4.05))
  arrow((6.02, 4.05), (6.02, 3.45))
  arrow((11.02, 4.05), (11.02, 3.45))
  arrow((15.02, 4.05), (15.02, 3.45))

  arrow((18.4, 4.75), (18.4, 3.45), straight: true)
  arrow((21.2, 4.75), (21.2, 3.45), straight: true)

  // ---- Row 4: plotPCA / plotCorrelation / plotEnrichment / plotHeatmap
  // centers match their row-5 children exactly, so every connector below is vertical;
  // the PCA/Correlation/Enrichment trio is centered as a group under multiBamSummary

  node-rbox(4.52, 2.05, 3.0, 1.4, [plotPCA], BLUE)
  node-rbox(9.22, 2.05, 3.6, 1.4, [plotCorrelation], BLUE)
  node-rbox(13.22, 2.05, 3.6, 1.4, [plotEnrichment], BLUE)
  node-rbox(18.1, 2.05, 3.4, 1.4, [plotHeatmap], BLUE)

  wire((6.02, 2.05), (6.02, 1.55))
  wire((4.62, 1.55), (7.42, 1.55))
  arrow((4.62, 1.55), (4.62, 1.1))
  arrow((7.42, 1.55), (7.42, 1.1))

  arrow((11.02, 2.05), (11.02, 1.1))
  arrow((15.02, 2.05), (15.02, 1.1))

  arrow((19.8, 2.05), (19.8, 1.1))

  // ---- Row 5: final outputs -- figure spots left empty for now

  let scatter = (3.32, -2.1, 2.6, 3.2)
  node-box(..scatter, subtitle: [PCA scatter], image-path: "subfigs/pca_scatter.png", subtitle-color: INK)

  let scree = (6.12, -2.1, 2.6, 3.2)
  node-box(..scree, subtitle: [scree plot], image-path: "subfigs/scree.png", subtitle-color: INK)

  let corr-heatmap = (9.42, -2.3, 3.2, 3.4)
  node-box(..corr-heatmap, subtitle: [correlation heatmap], image-path: "subfigs/corr_heatmap.png", subtitle-color: INK)

  let bar = (13.42, -2.3, 3.2, 3.4)
  node-box(..bar, subtitle: [bar chart], image-path: "subfigs/enrichment_bar.png", subtitle-color: INK)

  let heatmap = (18.1, -2.3, 3.4, 3.4)
  node-box(..heatmap, subtitle: [heatmap], image-path: "subfigs/heatmap_example.png", subtitle-color: INK)

  // --------------------------------------------------- bottom section brackets

  bracket(2.12, 21.8, -2.8, [QC / Visualization], flip: true)

  // --------------------------------------------------- headline accolades

  vbracket(23.0, 8.2, 16.9, [new Rust backend], ORANGE)
  vbracket(23.0, -3.5, 3.3, [improved visualizations], BLUE)
})
