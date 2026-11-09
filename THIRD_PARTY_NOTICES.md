# Third-party notices

The MIT license in `LICENSE` applies to this repository's own code. The
following components are either bundled or required separately.

## Bundled web viewers

- PDBe Mol* 3.10.1: Apache License 2.0. The complete notice is in
  `web/vendor/PDBe-Molstar-LICENSE`.
- 3Dmol.js: BSD-3-Clause, including its documented GLmol, Three.js and jQuery
  notices. The complete notice is in `web/vendor/3Dmol-LICENSE`.

## User-installed scientific tools

Users obtain these tools from their official sources and remain responsible
for the applicable terms. The SPURS inference patch in `docs/patches` is
derived from MIT-licensed SPURS code; its upstream license is included there.
Models, weights and executables are not bundled.

- SPURS: MIT-licensed code (<https://github.com/luo-group/SPURS>). The
  `cyclization9/SPURS` model card declares MIT.
- ESM-2 weights: obtained from Meta's official ESM download endpoint. Refer to
  the ESM repository and model terms before redistribution.
- DiscoTope-3.0: Creative Commons Attribution-NonCommercial 4.0 International.
  Commercial use requires separate permission from the upstream authors.
- ProCaliper: GPL-3.0-only according to its PyPI metadata. It is invoked as a
  separate optional worker and is not bundled here.
- FoldX 5.1: registration and acceptance of the FoldX license are required.
  Redistribution and commercial use are governed by that license.

Check the upstream repositories and license files again before publishing a
fork or offering a commercial service.
