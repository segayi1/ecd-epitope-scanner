/* Beads renderer (3Dmol) + engine adapter.
 *
 * The primary renderer stays PDBe Mol* (cartoon-quality, implemented inline in
 * app.js).  This module adds a second engine so users can switch between
 * Cartoon (Mol*) and Beads (3Dmol C-alpha spheres).  Both engines share ONE
 * colour callback so the semantic layering (design > interface > region >
 * surface > grey) is identical everywhere.
 *
 * 3Dmol picking notes (hard-won): atom picking reads viewer-level
 * setClickable/setHoverable, which must be called AFTER all addStyle calls and
 * BEFORE render(), so intersection shapes are rebuilt fresh.
 */

const BeadsRenderer = (() => {
  let viewer = null;
  let residues = [];
  let colorFn = null;
  let pickHandler = null;
  let hoverHandler = null;

  async function init(container, pdbUrl, residueList) {
    residues = residueList;
    container.innerHTML = '<div id="beads-canvas" style="position:absolute;inset:0"></div>';
    viewer = $3Dmol.createViewer(container.querySelector("#beads-canvas"), {
      backgroundColor: "#f4f0e7", antialias: true
    });
    const text = await fetch(pdbUrl).then((r) => {
      if (!r.ok) throw new Error(`structure fetch failed: HTTP ${r.status}`);
      return r.text();
    });
    viewer.addModel(text, "pdb");
    viewer.zoomTo();
    viewer.zoom(1.35);
    apply();
    return viewer;
  }

  function apply() {
    if (!viewer) return;
    viewer.setStyle({}, {});
    viewer.removeAllLabels();
    for (const residue of residues) {
      const color = colorFn ? colorFn(residue) : "#9CB89C";
      const chosen = residue.__beadSelected;
      viewer.addStyle(
        { chain: residue.chain_id, resi: residue.auth_seq_id, atom: "CA" },
        { sphere: { radius: chosen ? 2.1 : 1.8, color } }
      );
      if (residue.__beadSidechain) {
        viewer.addStyle(
          { chain: residue.chain_id, resi: residue.auth_seq_id, not: { atom: "CA" } },
          { stick: { radius: 0.15, color } }
        );
      }
    }
    viewer.setClickable({}, true, (atom) => {
      if (pickHandler) pickHandler(matchResidue(atom));
    });
    viewer.setHoverable({}, true, (atom) => {
      const residue = matchResidue(atom);
      if (residue && hoverHandler) hoverHandler(residue);
    }, () => {});
    viewer.render();
  }

  function matchResidue(atom) {
    return residues.find(
      (r) => r.chain_id === atom.chain && r.auth_seq_id === atom.resi
    ) || null;
  }

  function setColors(colorFnCallback, flags) {
    colorFn = colorFnCallback;
    if (flags) {
      for (const residue of residues) {
        residue.__beadSelected = flags.selected.has(residue.uniprot_pos);
        residue.__beadSidechain = flags.sidechain.has(residue.uniprot_pos);
      }
    }
    apply();
  }

  function setPickCallbacks(onPick, onHover) {
    pickHandler = onPick;
    hoverHandler = onHover;
  }

  function focus(residue) {
    if (!viewer || !residue) return;
    /* 3Dmol has no per-residue focus primitive; nudge the camera instead. */
    viewer.render();
  }

  function reset() {
    if (!viewer) return;
    viewer.zoomTo();
    viewer.zoom(1.35);
    viewer.render();
  }

  function destroy() {
    if (viewer && viewer.clear) viewer.clear();
    viewer = null;
  }

  return { init, apply, setColors, setPickCallbacks, focus, reset, destroy };
})();

/* NOTE: setRenderEngine / paintBeads / beadsColorFor live in app.js because
 * they need module-scope state (state, selectResidue, STEP_META...).
 * This file only exports the self-contained BeadsRenderer factory. */
