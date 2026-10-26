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
  let chainView = null;
  const disposed = new WeakSet();

  function dispose(instance) {
    if (!instance?.plugin || disposed.has(instance)) return;
    disposed.add(instance);
    instance.plugin.dispose();
  }

  /* plugin.render() can resolve BEFORE the internal molstar Plugin has built
     its managers (canvas controls make this window obvious).  Painting before
     managers exist throws "Cannot read properties of undefined (reading
     'managers')", so every caller waits for this. */
  async function waitReady(timeoutMs = 15000, instance = plugin) {
    const t0 = Date.now();
    while (Date.now() - t0 < timeoutMs) {
      if (instance !== plugin) return false;
      if (instance?.plugin?.managers && instance.assemblyRef
          && instance.plugin.managers.structure.hierarchy.current.structures.length > 0) return true;
      await new Promise((resolve) => setTimeout(resolve, 100));
    }
    throw Object.assign(new Error("Mol* structure loading timed out"), { stage: "structure" });
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

  async function init(host, structureUrl, representation, inlineStructureText = null) {
    destroy();
    host.innerHTML = "";
    let instance;
    try { instance = new window.PDBeMolstarPlugin(); }
    catch (error) { error.stage = "webgl"; throw error; }
    plugin = instance;
    initOptions = {
      customData: { url: structureUrl, format: "cif", binary: false },
      defaultPreset: "default",
      visualStyle: { polymer: representation || "cartoon", het: "ball-and-stick", carbs: "carbohydrate" },
      hideStructure: ["water"],
      bgColor: "#f1f5f7",
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
    let resourceUrl = null;
    let stage = "fetch";
    try {
      let data = inlineStructureText;
      if (data == null) {
        const response = await fetch(structureUrl);
        if (!response.ok) throw new Error(`HTTP ${response.status}`);
        data = await response.text();
      }
      stage = "parse";
      const format = /^\s*(?:#[^\n]*\n\s*)*data_/i.test(data) ? "cif"
        : /^(?:HEADER|ATOM\s|HETATM|MODEL\s|COMPND|TITLE\s|REMARK)/m.test(data) ? "pdb" : null;
      if (!format) throw new Error("Expected mmCIF or PDB structure data");
      if (instance !== plugin) return null;
      resourceUrl = URL.createObjectURL(new Blob([data], { type: "text/plain" }));
      const options = { ...initOptions, customData: { ...initOptions.customData, url: resourceUrl, format } };
      stage = "render";
      await instance.render(host, options);
      if (!await waitReady(15000, instance)) return null;
      applyViewProps();
      ready = true;
      window.__ecdMolstar = instance;
      return instance;
    } catch (error) {
      error.stage ||= /webgl|context|gpu/i.test(error.message) ? "webgl" : stage;
      throw error;
    } finally {
      if (resourceUrl) URL.revokeObjectURL(resourceUrl);
      if (instance !== plugin) dispose(instance);
    }
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
    const instance = plugin;
    const visualStyle = { polymer: representation, het: "ball-and-stick", carbs: "carbohydrate" };
    /* PDBe Mol* 3.10.1 `visual.update({ visualStyle }, false)` only updates
       init params; it does not rebuild existing representations. The pinned
       build exposes the same applyVisualStyles method used during load. Keep
       a full-reload fallback for a future build that hides that method. */
    if (typeof instance.applyVisualStyles === "function" && instance.assemblyRef) {
      await instance.applyVisualStyles(instance.assemblyRef, visualStyle);
    } else {
      const previousChains = chainView;
      initOptions = { ...initOptions, visualStyle };
      await instance.visual.update(initOptions, true);
      if (!await waitReady(15000, instance)) return;
      if (previousChains) {
        await configureChains(previousChains.targetChains, previousChains.partnerChains, representation, previousChains.targetRecords);
        chainView.visible = previousChains.visible;
      }
    }
    if (instance !== plugin) return;
    applyChainVisibility();
    await instance.canvas.applySettings({
      color: "#f1f5f7",
      lighting: representation === "molecular-surface" ? "matte" : "flat",
    });
    if (instance === plugin) applyViewProps();
  }

  async function color(queryRecords, nonSelectedColor, { keepRepresentations = false } = {}) {
    if (!plugin) return;
    if (!alive()) return;
    const instance = plugin;
    await instance.visual.select({ data: queryRecords, nonSelectedColor, keepRepresentations });
    if (instance === plugin) applyChainVisibility();
  }

  async function configureChains(targetChains, partnerChains, representation, targetRecords = null) {
    if (!alive() || (!targetChains.length && !targetRecords?.length)) return;
    const instance = plugin;
    const hierarchy = instance.plugin.managers.structure.hierarchy;
    const originals = hierarchy.current.structures.flatMap((structure) => structure.components);
    const query = (chains) => {
      const tests = chains.map((chain) => ({ head: { name: "core.rel.eq" }, args: [
        { head: { name: "structure-query.atom-property.macromolecular.auth_asym_id" }, args: [] }, chain,
      ] }));
      return { head: { name: "structure-query.generator.atom-groups" }, args: {
        "chain-test": tests.length === 1 ? tests[0] : { head: { name: "core.logic.or" }, args: tests },
      } };
    };
    const components = [];
    const groups = [["target", targetChains], ...(partnerChains.length ? [["partner", partnerChains]] : [])];
    for (const [key, chains] of groups) {
      const component = key === "target" && targetRecords?.length
        ? await instance.plugin.builders.structure.tryCreateComponent(instance.assemblyRef, {
          type: { name: "bundle", params: instance.getBundle(targetRecords) },
          nullIfEmpty: true, label: "Analysis target",
        }, "ecd-target", ["structure-component-static-polymer"])
        : await instance.plugin.builders.structure.tryCreateComponentFromExpression(
        instance.assemblyRef, query(chains), `ecd-${key}`, {
          label: key === "target" ? "Analysis target" : "Protein partner",
          tags: ["structure-component-static-polymer"],
        });
      if (component) {
        await instance.plugin.builders.structure.representation.addRepresentation(component, { type: representation });
        components.push(component);
      }
    }
    if (instance !== plugin) return;
    if (components.length !== groups.length) throw new Error("Analysis target or partner chain is missing from the display structure");
    chainView = { originals, targetChains, targetRecords, partnerChains, partner: components[1], visible: false };
    applyChainVisibility();
  }

  function applyChainVisibility() {
    if (!alive() || !chainView) return;
    const hierarchy = plugin.plugin.managers.structure.hierarchy;
    // Mol* rebuilds representations as visible even under a hidden component.
    hierarchy.toggleVisibility(chainView.originals, "hide");
    if (chainView.partner) hierarchy.toggleVisibility([chainView.partner], chainView.visible ? "show" : "hide");
  }

  async function setPartnerVisible(visible) {
    if (!alive() || !chainView) return;
    chainView.visible = visible;
    applyChainVisibility();
    await reset();
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

  async function focus(records, { durationMs = 600 } = {}) {
    if (!plugin) return;
    if (!alive()) return;
    const reducedMotion = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches === true;
    const loci = plugin.getLociForParams(records);
    // PDBe-Mol* 3.10.1 visual.focus() calls focusLoci() without options, whose
    // pinned default is 250 ms. Pass an explicit transition to retain the
    // current camera orientation while smoothly moving the focus target.
    await plugin.plugin.managers.camera.focusLoci(loci, { durationMs: reducedMotion ? 0 : durationMs });
  }

  async function reset() {
    if (!plugin) return;
    if (!alive()) return;
    if (chainView) {
      const target = chainView.targetRecords || chainView.targetChains.map((chain) => ({ auth_asym_id: chain }));
      await focus([...target, ...(chainView.visible ? chainView.partnerChains.map((chain) => ({ auth_asym_id: chain })) : [])]);
      return;
    }
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
    try { dispose(plugin); } catch (e) { console.warn("Mol* dispose failed", e); }
    plugin = null;
    initOptions = null;
    chainView = null;
    ready = false;
    window.__ecdMolstar = null;
  }

  return { init, update, waitReady, setRepresentation, configureChains, setPartnerVisible, color, highlight, clearHighlight, tooltips, focus, reset, resize, spin, destroy, isReady: () => ready && alive() };
})();

const BeadsView = (() => {
  let viewer = null;
  let residues = [];
  let pickCb = null;
  let hoverCb = null;

  async function init(host, pdbUrl, residueList, inlinePdbText = null) {
    destroy();
    residues = residueList;
    host.innerHTML = '<div class="beads-canvas" style="position:absolute;inset:0"></div>';
    const instance = $3Dmol.createViewer(host.querySelector(".beads-canvas"), {
      backgroundColor: "#f1f5f7", antialias: true,
    });
    viewer = instance;
    const text = inlinePdbText ?? await fetch(pdbUrl).then((r) => {
      if (!r.ok) throw new Error(`structure fetch failed: HTTP ${r.status}`);
      return r.text();
    });
    if (instance !== viewer) return null;
    viewer.addModel(text, "pdb");
    viewer.zoomTo();
    viewer.zoom(1.35);
    apply(() => "#9CB89C", { selected: new Set(), sidechain: new Set() });
    return instance;
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
