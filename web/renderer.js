/* ECD Epitope Workbench — renderer layer (rewrite v3)
 *
 * Two engines behind one tiny contract:
 *   MolstarView  : PDBe Mol* plugin (cartoon quality, built-in orbit/zoom/pan)
 *   BeadsView    : 3Dmol C-alpha spheres (semantic colouring)
 *
 * Neither module touches app state: app.js passes everything in.
 *
 * Mol* picking notes: residue picking comes through DOM events
 *   "PDB.molstar.click" / "PDB.molstar.mouseover" on the host element;
 *   the payload is event.eventData / event.detail with auth_seq_id.
 * 3Dmol picking notes: setClickable/setHoverable MUST be called after all
 *   addStyle calls and BEFORE render(), or the pickers stay inert.
 */

const MolstarView = (() => {
  let plugin = null;
  let ready = false;
  let initOptions = null;

  /* plugin.render() can resolve BEFORE the internal molstar Plugin has built
     its managers (canvas controls make this window obvious).  Painting before
     managers exist throws "Cannot read properties of undefined (reading
     'managers')", so every caller waits for this. */
  async function waitReady(timeoutMs = 15000) {
    const t0 = Date.now();
    while (Date.now() - t0 < timeoutMs) {
      if (plugin && plugin.plugin && plugin.plugin.managers
          && plugin.assemblyRef
          && plugin.plugin.managers.structure.hierarchy.current.structures.length > 0) return;
      await new Promise((resolve) => setTimeout(resolve, 100));
    }
    throw new Error("Mol* structure loading timed out");
  }

  /* Mol* defaults ship with distance fog (intensity 15) and far-plane
     clipping — zooming in makes the molecule literally fade/cut into the
     background colour, which users report as "结构被背景遮挡". Both off. */
  function applyViewProps() {
    const c3d = plugin?.plugin?.canvas3d;
    if (!c3d) return;
    c3d.setProps({
      cameraFog: { name: "off", params: {} },
      cameraClipping: { far: false, radius: 0, minNear: 1 },
    });
  }

  const alive = () => Boolean(plugin && plugin.plugin && plugin.plugin.managers);

  async function init(host, structureUrl, representation) {
    ready = false;
    host.innerHTML = "";
    plugin = new window.PDBeMolstarPlugin();
    initOptions = {
      customData: { url: structureUrl, format: "cif", binary: false },
      defaultPreset: "default",
      visualStyle: { polymer: representation || "cartoon", het: "ball-and-stick", carbs: "carbohydrate" },
      hideStructure: ["water"],
      bgColor: "#f4f0e7",
      highlightColor: "#FFD166",
      selectColor: "#FFD166",
      granularity: "residue",
      hideControls: true,          // left panel only; canvas controls stay ON
      sequencePanel: false,
      leftPanel: false,
      rightPanel: false,
      logPanel: false,
      pdbeLink: false,
      loadingOverlay: true,
      reactive: true,
      subscribeEvents: false,
    };
    await plugin.render(host, initOptions);
    await waitReady();
    applyViewProps();
    ready = true;
    window.__ecdMolstar = plugin;   // debug/inspection handle
    return plugin;
  }

  async function update(structureUrl) {
    if (!plugin) return;
    if (!alive()) return;
    await plugin.visual.update(
      { customData: { url: structureUrl, format: "cif", binary: false }, alphafoldView: false },
      true
    );
  }

  async function setRepresentation(representation) {
    if (!plugin) return;
    if (!alive()) return;
    const visualStyle = { polymer: representation, het: "ball-and-stick", carbs: "carbohydrate" };
    /* PDBe Mol* 3.10.1 `visual.update({ visualStyle }, false)` only updates
       init params; it does not rebuild existing representations. The pinned
       build exposes the same applyVisualStyles method used during load. Keep
       a full-reload fallback for a future build that hides that method. */
    if (typeof plugin.applyVisualStyles === "function" && plugin.assemblyRef) {
      await plugin.applyVisualStyles(plugin.assemblyRef, visualStyle);
    } else {
      initOptions = { ...initOptions, visualStyle };
      await plugin.visual.update(initOptions, true);
      await waitReady();
    }
    await plugin.canvas.applySettings({
      color: "#f4f0e7",
      lighting: representation === "molecular-surface" ? "matte" : "flat",
    });
    applyViewProps();
  }

  async function color(queryRecords, nonSelectedColor) {
    if (!plugin) return;
    if (!alive()) return;
    await plugin.visual.select({ data: queryRecords, nonSelectedColor, keepRepresentations: false });
  }

  async function highlight(record) {
    if (!plugin) return;
    if (!alive()) return;
    await plugin.visual.highlight({ data: [record], color: "#FFD166", focus: false });
  }

  async function clearHighlight() {
    if (!plugin) return;
    if (!alive()) return;
    await plugin.visual.clearHighlight();
  }

  async function tooltips(records) {
    if (!plugin) return;
    if (!alive()) return;
    await plugin.visual.tooltips({ data: records });
  }

  async function focus(records) {
    if (!plugin) return;
    if (!alive()) return;
    await plugin.visual.focus(records);
  }

  async function reset() {
    if (!plugin) return;
    if (!alive()) return;
    await plugin.visual.reset({ camera: true });
  }

  function resize() {
    /* Mol* listens for window resize; the canvas parent can change size
       without a window event (e.g. the sequence dock expanding), so we
       re-dispatch one. */
    window.dispatchEvent(new Event("resize"));
  }

  async function spin(on) {
    if (!plugin) return;
    if (!alive()) return;
    await plugin.visual.toggleSpin(on, false);
  }

  function destroy() {
    try { if (plugin && plugin.plugin) plugin.plugin.dispose(); } catch (e) { /* noop */ }
    plugin = null;
    initOptions = null;
    ready = false;
  }

  return { init, update, waitReady, setRepresentation, color, highlight, clearHighlight, tooltips, focus, reset, resize, spin, destroy, isReady: () => ready && alive() };
})();

const BeadsView = (() => {
  let viewer = null;
  let residues = [];
  let pickCb = null;
  let hoverCb = null;

  async function init(host, pdbUrl, residueList) {
    residues = residueList;
    host.innerHTML = '<div class="beads-canvas" style="position:absolute;inset:0"></div>';
    viewer = $3Dmol.createViewer(host.querySelector(".beads-canvas"), {
      backgroundColor: "#f4f0e7", antialias: true,
    });
    const text = await fetch(pdbUrl).then((r) => {
      if (!r.ok) throw new Error(`structure fetch failed: HTTP ${r.status}`);
      return r.text();
    });
    viewer.addModel(text, "pdb");
    viewer.zoomTo();
    viewer.zoom(1.35);
    apply(() => "#9CB89C", { selected: new Set(), sidechain: new Set() });
    return viewer;
  }

  function apply(colorFn, flags) {
    if (!viewer) return;
    viewer.setStyle({}, {});
    viewer.removeAllLabels();
    for (const residue of residues) {
      const color = colorFn(residue);
      const selected = flags.selected.has(residue.uniprot_pos);
      const spec = selected
        ? { sphere: { radius: 2.1, color } }
        : { sphere: { radius: 1.8, color } };
      viewer.addStyle({ chain: residue.chain_id, resi: residue.auth_seq_id, atom: "CA" }, spec);
      if (flags.sidechain.has(residue.uniprot_pos)) {
        viewer.addStyle(
          { chain: residue.chain_id, resi: residue.auth_seq_id, not: { atom: "CA" } },
          { stick: { radius: 0.15, color } }
        );
      }
    }
    /* arming order matters: styles first, picking second, render last */
    viewer.setClickable({}, true, (atom) => {
      const residue = match(atom);
      if (residue && pickCb) pickCb(residue);
    });
    viewer.setHoverable({}, true, (atom) => {
      const residue = match(atom);
      if (residue && hoverCb) hoverCb(residue);
    }, () => {
      if (hoverCb) hoverCb(null);
    });
    viewer.render();
  }

  function match(atom) {
    return residues.find(
      (r) => r.chain_id === atom.chain && r.auth_seq_id === atom.resi
    ) || null;
  }

  function setColors(colorFn, flags) { apply(colorFn, flags); }
  function setCallbacks(onPick, onHover) { pickCb = onPick; hoverCb = onHover; }
  function focusResidue() { /* 3Dmol: no per-residue camera primitive; no-op */ }
  function reset() { if (!viewer) return; viewer.zoomTo(); viewer.zoom(1.35); viewer.render(); }
  function resize() { if (!viewer) return; viewer.resize(); viewer.render(); }
  function destroy() { if (viewer && viewer.clear) viewer.clear(); viewer = null; }

  return { init, setColors, setCallbacks, focusResidue, reset, resize, destroy };
})();
